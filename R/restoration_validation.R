# Numerical restoration contracts and receipts, separate from BIDS structure.

.lna_error_contract_schema_version <- "1.0.0"
.lna_numerical_receipt_schema_version <- "1.0.0"

.numerical_abort <- function(message, field, receipt = NULL, ...) {
  abort_lna(
    message,
    ...,
    field = field,
    receipt = receipt,
    .subclass = "lna_error_numerical_restoration",
    location = paste0("validate_lna_restoration:", field)
  )
}

.numerical_scalar_number <- function(x, field) {
  if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0) {
    .numerical_abort(
      paste0(field, " must be one finite non-negative number"), field
    )
  }
  as.double(x)
}

.numerical_choice <- function(x, choices, field) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !x %in% choices) {
    .numerical_abort(
      paste0(field, " must be one of: ", paste(choices, collapse = ", ")),
      field
    )
  }
  x
}

#' Declare a numerical error contract for LNA restoration
#'
#' @param loss_mode `"lossless"`, `"lossy"`, or `"unspecified"`.
#' @param output_dtype Declared reconstructed storage dtype.
#' @param max_abs_error Maximum permitted pointwise absolute error.
#' @param rmse_max Maximum permitted root mean squared error.
#' @param affine_tolerance Maximum permitted affine-element absolute error.
#' @return A validated JSON-safe `lna_error_contract`.
#' @export
lna_error_contract <- function(
    loss_mode = c("lossless", "lossy", "unspecified"),
    output_dtype = NULL,
    max_abs_error = NULL,
    rmse_max = NULL,
    affine_tolerance = 1e-6) {
  loss_mode <- match.arg(loss_mode)
  if (identical(loss_mode, "lossless")) {
    output_dtype <- output_dtype %||% "float64"
    max_abs_error <- max_abs_error %||% 0
    rmse_max <- rmse_max %||% 0
  } else if (identical(loss_mode, "unspecified")) {
    output_dtype <- output_dtype %||% "float32"
    max_abs_error <- max_abs_error %||% .Machine$double.xmax
    rmse_max <- rmse_max %||% .Machine$double.xmax
  }
  validate_lna_error_contract(list(
    schema_version = .lna_error_contract_schema_version,
    loss_mode = loss_mode,
    output_dtype = output_dtype,
    value_domain = "finite-reconstructed-values",
    max_abs_error = max_abs_error,
    rmse_max = rmse_max,
    affine_tolerance = affine_tolerance,
    bit_exact = identical(loss_mode, "lossless")
  ))
}

#' @rdname lna_error_contract
#' @param contract A contract or decoded JSON object.
#' @export
validate_lna_error_contract <- function(contract) {
  if (inherits(contract, "lna_error_contract")) contract <- unclass(contract)
  required <- c(
    "schema_version", "loss_mode", "output_dtype", "value_domain",
    "max_abs_error", "rmse_max", "affine_tolerance", "bit_exact"
  )
  if (!is.list(contract) || !setequal(names(contract), required)) {
    .numerical_abort(
      paste0("contract must contain exactly: ", paste(required, collapse = ", ")),
      "contract"
    )
  }
  if (!identical(contract$schema_version, .lna_error_contract_schema_version)) {
    .numerical_abort("unsupported error-contract schema version", "schema_version")
  }
  loss_mode <- .numerical_choice(
    contract$loss_mode, c("lossless", "lossy", "unspecified"),
    "contract.loss_mode"
  )
  output_dtype <- .numerical_choice(
    contract$output_dtype, c("float32", "float64", "integer32"),
    "contract.output_dtype"
  )
  if (!identical(contract$value_domain, "finite-reconstructed-values")) {
    .numerical_abort("unsupported error-contract value domain", "value_domain")
  }
  bit_exact <- contract$bit_exact
  if (!is.logical(bit_exact) || length(bit_exact) != 1L || is.na(bit_exact)) {
    .numerical_abort("contract.bit_exact must be one logical value", "bit_exact")
  }
  max_abs <- .numerical_scalar_number(contract$max_abs_error, "max_abs_error")
  rmse <- .numerical_scalar_number(contract$rmse_max, "rmse_max")
  affine <- .numerical_scalar_number(
    contract$affine_tolerance, "affine_tolerance"
  )
  if (identical(loss_mode, "lossless") &&
      (!isTRUE(bit_exact) || max_abs != 0 || rmse != 0)) {
    .numerical_abort(
      "lossless contracts require bit_exact and zero value-error bounds",
      "contract"
    )
  }
  if (identical(loss_mode, "lossy") && isTRUE(bit_exact)) {
    .numerical_abort("lossy contracts may not claim bit exactness", "bit_exact")
  }
  structure(list(
    schema_version = .lna_error_contract_schema_version,
    loss_mode = loss_mode,
    output_dtype = output_dtype,
    value_domain = "finite-reconstructed-values",
    max_abs_error = max_abs,
    rmse_max = rmse,
    affine_tolerance = affine,
    bit_exact = isTRUE(bit_exact)
  ), class = c("lna_error_contract", "list"))
}

.numerical_dtype_from_header <- function(header, values) {
  code <- header$datatype %||% NULL
  if (!is.null(code)) {
    code <- suppressWarnings(as.integer(code))
    mapped <- c(
      `8` = "integer32", `16` = "float32", `64` = "float64",
      `768` = "integer32", `1024` = "integer32"
    )
    if (as.character(code) %in% names(mapped)) return(unname(mapped[[as.character(code)]]))
  }
  if (is.integer(values)) "integer32" else "float64"
}

.numerical_object <- function(x, affine = NULL, field = "object") {
  original_path <- NULL
  if (is.character(x) && length(x) == 1L && !is.na(x)) {
    original_path <- x
    if (!file.exists(x)) .numerical_abort(paste0(field, " file does not exist"), field)
    x <- if (endsWith(x, ".lna.h5") || endsWith(x, ".h5")) {
      as_neurovec(x, output_dtype = "float64")
    } else {
      neuroim2::read_vec(x)
    }
  }
  is_neuro <- methods::is(x, "NeuroObj")
  values <- if (is_neuro) as.array(x) else x
  if (!is.numeric(values) || is.null(dim(values)) ||
      !length(values) || any(!is.finite(values))) {
    .numerical_abort(
      paste0(field, " must contain a finite numeric image array"), field
    )
  }
  if (!length(dim(values)) %in% c(3L, 4L)) {
    .numerical_abort(paste0(field, " must be three- or four-dimensional"), field)
  }
  if (length(dim(values)) == 3L) {
    values <- array(values, dim = c(dim(values), 1L))
  }
  object_affine <- affine
  header <- list()
  if (is_neuro) {
    object_affine <- neuroim2::trans(neuroim2::space(x))
    if (methods::.hasSlot(x, "header")) header <- methods::slot(x, "header")
  }
  object_affine <- if (is.null(object_affine)) NULL else .lna_normalize_affine(object_affine)
  list(
    values = values,
    dims = as.integer(dim(values)),
    affine = object_affine,
    dtype = .numerical_dtype_from_header(header, values),
    path = original_path
  )
}

.numerical_support <- function(support, dims, field) {
  if (is.null(support)) return(array(TRUE, dim = dims[1:3]))
  normalized <- .lna_normalize_support(support)
  if (length(normalized) != prod(dims[1:3])) {
    .numerical_abort(
      paste0(field, " support does not match spatial dimensions"), field
    )
  }
  array(normalized, dim = dims[1:3])
}

.numerical_values_on_support <- function(values, support) {
  matrix <- matrix(values, nrow = prod(dim(values)[1:3]))
  as.double(matrix[as.vector(support), , drop = FALSE])
}

.numerical_double_bytes <- function(values) {
  connection <- rawConnection(raw(), "wb")
  on.exit(close(connection), add = TRUE)
  writeBin(as.double(values), connection, size = 8L, endian = "little")
  rawConnectionValue(connection)
}

#' Validate decoded numerical fidelity independently of BIDS structure
#'
#' @param reference Reference array, neuroim2 image, NIfTI path, or LNA path.
#' @param restored Restored array, neuroim2 image, NIfTI path, or LNA path.
#' @param contract A declared [lna_error_contract()].
#' @param support Optional logical/binary spatial support defining the compared
#'   estimand.
#' @param reference_affine,restored_affine Explicit affines for plain arrays.
#' @param expected_support_id Optional canonical support identity to verify.
#' @param record_id Optional record identity included in the receipt.
#' @param fail Whether a failed check raises an error carrying the receipt.
#' @return An `lna_numerical_receipt` with pointwise and aggregate metrics.
#' @export
validate_lna_restoration <- function(
    reference,
    restored,
    contract,
    support = NULL,
    reference_affine = NULL,
    restored_affine = NULL,
    expected_support_id = NULL,
    record_id = NULL,
    fail = TRUE) {
  contract <- validate_lna_error_contract(contract)
  reference <- .numerical_object(reference, reference_affine, "reference")
  restored <- .numerical_object(restored, restored_affine, "restored")
  dims_equal <- identical(reference$dims, restored$dims)
  comparison_dims <- if (dims_equal) reference$dims else reference$dims
  support <- .numerical_support(support, comparison_dims, "support")
  reference_values <- .numerical_values_on_support(reference$values, support)
  restored_values <- if (dims_equal) {
    .numerical_values_on_support(restored$values, support)
  } else numeric()
  errors <- if (dims_equal) restored_values - reference_values else Inf
  max_abs <- if (length(errors)) max(abs(errors)) else Inf
  rmse <- if (length(errors)) sqrt(mean(errors^2)) else Inf
  mean_abs <- if (length(errors)) mean(abs(errors)) else Inf
  value_exact <- dims_equal && identical(
    .numerical_double_bytes(reference_values),
    .numerical_double_bytes(restored_values)
  )
  have_affines <- !is.null(reference$affine) && !is.null(restored$affine)
  affine_error <- if (have_affines) {
    max(abs(reference$affine - restored$affine))
  } else if (is.null(reference$affine) && is.null(restored$affine)) 0 else Inf
  grid_id <- if (!is.null(reference$affine)) {
    .lna_grid_id(
      reference$dims[1:3], reference$affine, "mm", "first-axis-fastest"
    )
  } else NULL
  support_id <- if (!is.null(grid_id)) .lna_support_id(grid_id, support) else NULL
  checks <- list(
    dimensions = dims_equal,
    temporal_extent = dims_equal && reference$dims[[4L]] == restored$dims[[4L]],
    dtype = identical(restored$dtype, contract$output_dtype),
    affine = affine_error <= contract$affine_tolerance,
    support = is.null(expected_support_id) || identical(support_id, expected_support_id),
    max_abs_error = max_abs <= contract$max_abs_error,
    rmse = rmse <= contract$rmse_max,
    bit_exact = !isTRUE(contract$bit_exact) || value_exact
  )
  passed <- all(unlist(checks, use.names = FALSE))
  receipt <- structure(list(
    schema_version = .lna_numerical_receipt_schema_version,
    record_id = record_id,
    passed = passed,
    contract = unclass(contract),
    reference = list(
      dims = reference$dims,
      dtype = reference$dtype,
      affine_id = if (is.null(reference$affine)) NULL else .lna_affine_id(reference$affine),
      path = reference$path
    ),
    restored = list(
      dims = restored$dims,
      dtype = restored$dtype,
      affine_id = if (is.null(restored$affine)) NULL else .lna_affine_id(restored$affine),
      path = restored$path
    ),
    support_id = support_id,
    n_compared_values = length(reference_values),
    metrics = list(
      max_abs_error = max_abs,
      mean_abs_error = mean_abs,
      rmse = rmse,
      affine_max_abs_error = affine_error,
      value_exact = value_exact
    ),
    checks = checks,
    serialization = list(
      source_and_nifti_bytes_compared = FALSE,
      reason = "NIfTI serialization bytes are outside the decoded-value estimand"
    )
  ), class = c("lna_numerical_receipt", "list"))
  if (!passed && isTRUE(fail)) {
    failed <- names(checks)[!unlist(checks, use.names = FALSE)]
    .numerical_abort(
      paste0("numerical restoration failed: ", paste(failed, collapse = ", ")),
      "receipt", receipt = receipt, failed_checks = failed
    )
  }
  receipt
}

#' @export
print.lna_numerical_receipt <- function(x, ...) {
  cat("<lna_numerical_receipt>", if (x$passed) "PASS" else "FAIL", "\n")
  cat("  values:", x$n_compared_values, "\n")
  cat("  max abs error:", format(x$metrics$max_abs_error), "\n")
  cat("  RMSE:", format(x$metrics$rmse), "\n")
  invisible(x)
}
