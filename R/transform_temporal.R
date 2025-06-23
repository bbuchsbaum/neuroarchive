#' Temporal Transform - Core Infrastructure
#'
#' This file contains the core temporal transform infrastructure including
#' forward_step and invert_step functions. Individual basis types, projection,
#' and reconstruction methods are now organized in separate files:
#' 
#' - temporal_basis.R: Basis generation for DCT, B-spline, DPSS, polynomial, wavelet
#' - temporal_modwt.R: All MODWT-specific functions and methods  
#' - temporal_project.R: Projection methods with thresholding
#' - temporal_reconstruct.R: Reconstruction methods
#' - temporal_utils.R: Utility functions like suggest_dpss_fmri()

#' Temporal Transform - Forward Step
#'
#' Projects data onto a temporal basis (DCT, B-spline, DPSS, polynomial, or wavelet).
#' Debug messages are controlled by the `lna.debug.temporal` option.
#' @keywords internal
#' @export
forward_step.temporal <- function(type, desc, handle) {
  dbg <- isTRUE(getOption("lna.debug.temporal", FALSE))
  p <- desc$params %||% list()
  # Extract temporal-specific parameters and remove them from p to avoid duplication
  kind <- p$kind %||% "dct"
  n_basis <- p$n_basis
  p$kind <- NULL
  p$n_basis <- NULL
  order <- p$order %||% 3
  p$order <- NULL

  if (!is.null(n_basis)) {
    if (!is.numeric(n_basis) || length(n_basis) != 1 ||
        n_basis <= 0 || n_basis %% 1 != 0) {
      abort_lna(
        "n_basis must be a positive integer",
        .subclass = "lna_error_validation",
        location = "forward_step.temporal:n_basis"
      )
    }
    n_basis <- as.integer(n_basis)
  }

  if (!is.null(order)) {
    if (!is.numeric(order) || length(order) != 1 ||
        order <= 0 || order %% 1 != 0) {
      abort_lna(
        "order must be a positive integer",
        .subclass = "lna_error_validation",
        location = "forward_step.temporal:order"
      )
    }
    order <- as.integer(order)
  }
  # Determine input key based on previous transform's output.
  # When temporal coefficients are already present we treat them as
  # the input so additional temporal steps operate on the projected
  # coefficients rather than reusing the raw matrix.
  if (handle$has_key("temporal_coefficients")) {
    input_key <- "temporal_coefficients"
  } else if (handle$has_key("delta_stream")) {
    input_key <- "delta_stream"
  } else if (handle$has_key("sparsepca_embedding")) {
    input_key <- "sparsepca_embedding"
  } else if (handle$has_key("aggregated_matrix")) {
    input_key <- "aggregated_matrix"
  } else {
    input_key <- if (!is.null(desc$inputs)) desc$inputs[[1]] else "input"
  }

  X <- handle$get_inputs(input_key)[[1]]
  X <- as_dense_mat(X)

  n_time <- nrow(X)
  if (is.null(n_basis)) n_basis <- n_time
  n_basis <- min(n_basis, n_time)

  # After resolving defaults, store parameters back in desc$params
  p_final <- c(list(kind = kind, n_basis = n_basis, order = order), p)
  desc$params <- p_final

  args <- c(list(kind = kind, n_time = n_time, n_basis = n_basis, order = order),
            p)
  
  basis <- do.call(temporal_basis, args)

  # Delegate projection logic to per-kind methods for extensibility
  coeff <- temporal_project(kind, basis, X)

  if (dbg) {
    # DEBUG: Check reconstruction locally
    if (is.matrix(basis) && is.matrix(coeff) && ncol(basis) == nrow(coeff)) {
      if (identical(kind, "polynomial")) {
        message("[forward_step.temporal POLY DEBUG] Checking orthogonality of basis (t(basis) %*% basis):")
        # Ensure it's a plain matrix for printing, and round for clarity
        t_basis_basis <- as.matrix(crossprod(basis))
        print(round(t_basis_basis, 5))
      }
      X_reconstructed_debug <- basis %*% coeff # Should be time x features
      if (!isTRUE(all.equal(X, X_reconstructed_debug, tolerance = 1e-7))) {
        message("[forward_step.temporal DEBUG] Local reconstruction MISMATCH.")
        if (identical(kind, "polynomial")) {
           message("Sum of squared differences: ", sum((X - X_reconstructed_debug)^2))
        }
      } else {
        message("[forward_step.temporal DEBUG] Local reconstruction MATCHES.")
      }
    } else {
      message("[forward_step.temporal DEBUG] Could not perform local reconstruction check due to matrix non-conformance.")
    }
    # END DEBUG
  }

  run_id <- handle$current_run_id %||% "run-01"
  run_id <- sanitize_run_id(run_id)
  plan <- handle$plan
  fname <- plan$get_next_filename(type)
  base_name <- tools::file_path_sans_ext(fname)
  temporal_root <- lna_options("paths.temporal_root")[[1]]
  scans_root <- lna_options("paths.scans_root")[[1]]
  basis_path <- paste0(temporal_root, base_name, "/basis")
  coef_path <- paste0(scans_root, run_id, "/", base_name, "/coefficients")
  knots_path <- paste0(temporal_root, base_name, "/knots")
  params_json <- as.character(jsonlite::toJSON(desc$params, auto_unbox = TRUE))

  desc$version <- "1.0"
  desc$inputs <- c(input_key)
  desc$outputs <- c("temporal_coefficients")
  datasets <- list(
    list(path = basis_path, role = "temporal_basis"),
    list(path = coef_path, role = "temporal_coefficients")
  )
  knots_data <- attr(basis, "knots")
  if (!is.null(knots_data)) {
    datasets[[length(datasets) + 1]] <- list(path = knots_path, role = "knots")
  }
  desc$datasets <- datasets

  plan$add_descriptor(fname, desc)
  
  basis_payload <- basis
  plan$add_payload(basis_path, basis_payload)
  


  plan$add_dataset_def(basis_path, "temporal_basis", as.character(type), run_id,
                       as.integer(plan$next_index), params_json,
                       basis_path, "eager", dtype = "float32")

  if (!is.null(knots_data)) {
    knots_payload <- knots_data
    plan$add_payload(knots_path, knots_payload)
    plan$add_dataset_def(knots_path, "knots", as.character(type), run_id,
                         as.integer(plan$next_index), params_json,
                         knots_path, "eager", dtype = "float32")
  }
  
  coeff_payload <- coeff
  plan$add_payload(coef_path, coeff_payload)

  plan$add_dataset_def(coef_path, "temporal_coefficients", as.character(type), run_id,
                       as.integer(plan$next_index), params_json,
                       coef_path, "eager", dtype = "float32")
  handle$plan <- plan

  handle$update_stash(keys = c(input_key),
                      new_values = list(temporal_coefficients = coeff))
}

#' Temporal Transform - Inverse Step
#'
#' Reconstructs data from stored temporal basis coefficients.
#' Debug messages are controlled by the `lna.debug.temporal` option.
#' @keywords internal
#' @export
invert_step.temporal <- function(type, desc, handle) {
  dbg <- isTRUE(getOption("lna.debug.temporal", FALSE))
  if (dbg) message(sprintf("[invert_step.temporal ENTRY] Incoming handle stash keys: %s. Is input NULL? %s", paste(names(handle$stash), collapse=", "), is.null(handle$stash$input)))
  basis_path <- NULL
  coeff_path <- NULL
  
  if (!is.null(desc$datasets)) {
    roles <- vapply(desc$datasets, function(d) d$role, character(1))
    idx_b <- which(roles == "temporal_basis")
    if (length(idx_b) > 0) basis_path <- desc$datasets[[idx_b[1]]]$path
    idx_c <- which(roles == "temporal_coefficients")
    if (length(idx_c) > 0) coeff_path <- desc$datasets[[idx_c[1]]]$path
  }

  if (is.null(basis_path)) {
    abort_lna(
      "temporal_basis path not found in descriptor",
      .subclass = "lna_error_descriptor",
      location = "invert_step.temporal:basis_path"
    )
  }
  if (is.null(coeff_path)) {
    abort_lna("temporal_coefficients path not found in descriptor datasets", .subclass = "lna_error_descriptor", location = "invert_step.temporal")
  }

  output_stash_key  <- desc$inputs[[1]] %||% "input"

  root <- handle$h5[["/"]]
  basis <- h5_read(root, basis_path)
  coeff <- h5_read(root, coeff_path)

  if (dbg) {
    message(sprintf("[invert_step.temporal] Basis dims: %s", paste(dim(basis), collapse = "x")))
    message(sprintf("[invert_step.temporal] Coeff dims: %s", paste(dim(coeff), collapse = "x")))
  }
  

  
  if (dbg) {
    message("--- Invert Step Pre-Dense Calculation Debug ---")
    if (nrow(basis) >= 2 && ncol(basis) >= 2) {
      message("basis_loaded[1:2, 1:2]:")
      print(basis[1:2, 1:2, drop = FALSE])
    }
    if (nrow(coeff) >= 2 && ncol(coeff) >= 2) {
      message("coeff_loaded[1:2, 1:2]:")
      print(coeff[1:2, 1:2, drop = FALSE])
    }
  }
  
  # Check for valid matrix dimensions before multiplication
  if (!is.matrix(basis) || !is.matrix(coeff)) {
    # Handle case where basis/coeff are stored dimension vectors from empty arrays
    if (length(basis) == 2 && length(coeff) == 2 && 
        all(basis >= 0) && all(coeff >= 0) && 
        all(basis == as.integer(basis)) && all(coeff == as.integer(coeff))) {
      # These look like stored dimensions - reconstruct the original empty matrices
      basis_dims <- as.integer(basis)
      coeff_dims <- as.integer(coeff)
      
      
      
      # Reconstruct the matrices
      basis <- array(numeric(0), dim = basis_dims)
      coeff <- array(numeric(0), dim = coeff_dims)
      
      # Now check if we can do matrix multiplication
      if (basis_dims[2] != coeff_dims[1]) {
        # Dimensions don't match for multiplication - create empty result
        dense <- matrix(numeric(0), nrow = basis_dims[1], ncol = coeff_dims[2])
        
      } else {
        # Dimensions match - do the multiplication (which will result in empty matrix)
        dense <- basis %*% coeff
        
      }
    } else if (length(basis) == 0 && length(coeff) == 0) {
      # Both are empty - create empty result matrix
      dense <- matrix(numeric(0), nrow = 0, ncol = 0)
      
    } else {
      abort_lna("Invalid matrix dimensions for multiplication", .subclass="lna_error_internal", location="invert_step.temporal")
    }
  } else {
    # Special case: if basis is 0x0 (n_time=0, n_basis=0), create empty result with correct dimensions
    if (nrow(basis) == 0 && ncol(basis) == 0) {
      # When n_time=0, we need to reconstruct to match the expected output dimensions
      # The output should have 0 rows (time) and the same number of columns as the original data
      # We can infer the number of columns from the coefficient matrix
      n_features <- if (is.matrix(coeff) && ncol(coeff) > 0) ncol(coeff) else 1
      dense <- matrix(numeric(0), nrow = 0, ncol = n_features)
      
    } else {
      # Check that matrix dimensions are compatible for multiplication
      if (ncol(basis) != nrow(coeff)) {
        abort_lna(
          sprintf("Matrix dimension mismatch: basis has %d columns but coeff has %d rows", 
                  ncol(basis), nrow(coeff)),
          .subclass = "lna_error_internal", 
          location = "invert_step.temporal"
        )
      }

      dense <- temporal_reconstruct(desc$params$kind %||% "dct", basis, coeff)
    }
  }
  
  if (dbg) message(sprintf("[invert_step.temporal] Dense dims after matmult: %s", paste(dim(dense), collapse="x")))
  if (dbg && nrow(dense) >= 2 && ncol(dense) >= 2) {
    message("dense[1:2, 1:2]:")
    print(dense[1:2, 1:2, drop = FALSE])
  }
  
  subset <- handle$subset
  if (!is.null(subset$roi_mask)) {
    roi <- as.logical(subset$roi_mask)
    if (length(roi) == ncol(dense)) { 
      dense <- dense[, roi, drop = FALSE]
    }
  }
  if (!is.null(subset$time_idx)) {
    idx <- as.integer(subset$time_idx)
    # Ensure that idx is not empty and all indices are within bounds
    if (length(idx) > 0 && nrow(dense) >= max(idx) && min(idx) > 0) { 
        dense <- dense[idx, , drop = FALSE]
    } else if (length(idx) > 0) {
        # Handle out-of-bounds or empty idx if necessary, or let it error if that's desired.
        warning("time_idx for temporal subsetting is invalid or out of bounds.")
    }
  }
  if (dbg) message(sprintf("[invert_step.temporal] Dense dims after subsetting: %s", paste(dim(dense), collapse="x")))
  
  if (is.null(dense)) {
    abort_lna("Reconstructed data (dense) is NULL before stashing", .subclass="lna_error_internal", location="invert_step.temporal")
  }
  if (dbg) message(sprintf("[invert_step.temporal] Stashing to key: '%s'. Is dense NULL? %s", output_stash_key, is.null(dense)))
  new_values_list <- setNames(list(dense), output_stash_key)

  handle <- handle$update_stash(keys = character(),
                                new_values = new_values_list)
  if (dbg) message(sprintf("[invert_step.temporal] invert_step.temporal IS RETURNING handle with Stash keys: %s. Is input NULL? %s", paste(names(handle$stash), collapse=", "), is.null(handle$stash$input)))
  return(handle)
} 