#' Delta Transform - Forward Step
#' Computes first-order differences along a chosen axis. When
#' `coding_method` is set to `'rle'`, the differences are compressed using
#' the **rle** package and stored as a two column matrix with `lengths` and
#' `values`.
#'
#' @importFrom rle as.rle compress.rle inverse.rle
#' @keywords internal
forward_step.delta <- function(type, desc, handle) {
  p <- desc$params %||% list()
  order <- p$order %||% 1
  axis <- p$axis %||% -1
  ref_store <- p$reference_value_storage %||% "first_value_verbatim"
  coding <- p$coding_method %||% "none"
  if (!coding %in% c("none", "rle")) {
    abort_lna(
      sprintf("Unsupported coding_method '%s'", coding),
      .subclass = "lna_error_validation",
      location = "forward_step.delta:coding_method"
    )
  }

  if (!identical(order, 1L)) {
    abort_lna("only order=1 supported", .subclass = "lna_error_validation",
              location = "forward_step.delta:order")
  }

  # Determine input key primarily from the descriptor
  input_key <- desc$inputs[[1]]
  if (is.null(input_key)) {
    warning(
      "desc$inputs missing for delta forward step; defaulting to 'input'",
      call. = FALSE
    )
    input_key <- "input"
  }
  x <- handle$get_inputs(input_key)[[1]]
  dims <- dim(x)
  if (is.null(dims)) dims <- length(x)
  if (axis == -1) axis <- length(dims)
  axis <- as.integer(axis)
  stopifnot(axis >= 1, axis <= length(dims))

  perm <- c(axis, setdiff(seq_along(dims), axis))
  xp <- aperm(x, perm)
  dim(xp) <- c(dims[axis], prod(dims[-axis]))

  first_vals <- xp[1, , drop = FALSE]
  dim(first_vals) <- c(1, prod(dims[-axis]))
  deltas <- xp[-1, , drop = FALSE] - xp[-nrow(xp), , drop = FALSE]

  if (identical(coding, "rle")) {
    # Flatten the delta matrix to a vector for run-length encoding,
    # compress it, and store the run lengths and values as a 2-column matrix.
    vec <- as.vector(deltas)
    r_obj <- rle::as.rle(vec)
    r_obj <- rle::compress.rle(r_obj)
    delta_stream <- cbind(lengths = r_obj$lengths, values = r_obj$values)
  } else {
    delta_stream <- deltas
  }

  run_id <- handle$current_run_id %||% "run-01"
  run_id <- sanitize_run_id(run_id)
  plan <- handle$plan
  fname <- plan$get_next_filename(type)
  base <- tools::file_path_sans_ext(fname)
  delta_path <- paste0("/scans/", run_id, "/deltas/", base, "/delta_stream")
  first_path <- paste0("/scans/", run_id, "/deltas/", base, "/first_values")
  step_index <- plan$next_index
  p$orig_dims <- dims
  params_json <- as.character(jsonlite::toJSON(p, auto_unbox = TRUE))

  desc$version <- "1.0"
  desc$inputs <- c(input_key)

  output_key <- NULL
  if (!is.null(desc$outputs)) {
    if (length(desc$outputs) > 0) {
      output_key <- desc$outputs[[1]]
      desc$outputs <- output_key
    } else {
      desc$outputs <- character()
    }
  } else {
    output_key <- "delta_stream"
    desc$outputs <- output_key
  }

  desc$params <- p
  ds <- list(list(path = delta_path, role = "delta_stream"))
  if (identical(ref_store, "first_value_verbatim")) {
    ds[[length(ds) + 1]] <- list(path = first_path, role = "first_values")
  }
  desc$datasets <- ds

  plan$add_descriptor(fname, desc)
  plan$add_payload(delta_path, delta_stream)
  plan$add_dataset_def(delta_path, "delta_stream", as.character(type), run_id,
                       as.integer(step_index), params_json,
                       delta_path, "eager", dtype = NA_character_)
  if (identical(ref_store, "first_value_verbatim")) {
    plan$add_payload(first_path, first_vals)
    plan$add_dataset_def(first_path, "first_values", as.character(type), run_id,
                         as.integer(step_index), params_json,
                         first_path, "eager", dtype = NA_character_)
  }

  handle$plan <- plan
  if (!is.null(output_key)) {
    handle$update_stash(
      keys = c(input_key),
      new_values = setNames(list(delta_stream), output_key)
    )
  } else {
    handle$update_stash(keys = c(input_key), new_values = list())
  }
}

#' Delta Transform - Inverse Step
#'
#' Decodes the delta representation written by
#' \code{\link{forward_step.delta}}. When `coding_method` was
#' `'rle'`, the stored matrix of run lengths and values is expanded back to
#' the full set of differences before reconstruction.
#' @keywords internal
invert_step.delta <- function(type, desc, handle) {
  p <- desc$params %||% list()
  axis <- p$axis %||% -1
  ref_store <- p$reference_value_storage %||% "first_value_verbatim"
  coding <- p$coding_method %||% "none"
  if (!coding %in% c("none", "rle")) {
    abort_lna(
      sprintf("Unsupported coding_method '%s'", coding),
      .subclass = "lna_error_validation",
      location = "invert_step.delta:coding_method"
    )
  }
  dims <- p$orig_dims
  if (is.null(dims)) {
    abort_lna("orig_dims missing in descriptor", .subclass = "lna_error_descriptor",
              location = "invert_step.delta")
  }
  if (axis == -1) axis <- length(dims)
  axis <- as.integer(axis)

  run_id <- handle$current_run_id %||% "run-01"
  delta_path <- desc$datasets[[1]]$path
  first_path <- NULL
  if (identical(ref_store, "first_value_verbatim")) {
    idx <- which(vapply(desc$datasets, function(d) d$role, character(1)) == "first_values")
    if (length(idx) > 0) first_path <- desc$datasets[[idx[1]]]$path
  }

  root <- handle$h5[["/"]]
  delta_stream <- h5_read(root, delta_path)
  expected_ncols <- if (length(dims) == 1) 1 else prod(dims[-axis])

  if (identical(ref_store, "first_value_verbatim")) {
    first_vals <- h5_read(root, first_path)
  } else {
    first_vals <- array(0, dim = expected_ncols)
  }

  if (length(first_vals) != expected_ncols) {
    abort_lna(
      sprintf(
        "first_vals length (%d) mismatch, expected %d",
        length(first_vals), expected_ncols
      ),
      .subclass = "lna_error_runtime",
      location = "invert_step.delta:first_vals_length"
    )
  }
  dim(first_vals) <- c(1, expected_ncols)

  expected_nrows_deltas <- max(0, dims[axis] - 1L)

  if (identical(coding, "rle")) {
    # Reconstruct the delta vector from run-length encoding stored in the
    # two-column matrix produced in the forward step.
    lengths <- delta_stream[, 1]
    values <- delta_stream[, 2]
    r_obj <- structure(list(lengths = lengths, values = values), class = "rle")
    delta_vec <- rle::inverse.rle(r_obj)
    if (expected_nrows_deltas > 0 && length(delta_vec) != expected_nrows_deltas * expected_ncols) {
      abort_lna(
        sprintf(
          "RLE decoded data length (%d) mismatch. Expected %d elements for %dx%d deltas matrix.",
          length(delta_vec), expected_nrows_deltas * expected_ncols,
          expected_nrows_deltas, expected_ncols
        ),
        .subclass = "lna_error_runtime",
        location = "invert_step.delta:rle_decode"
      )
    }
    deltas <- matrix(delta_vec, nrow = expected_nrows_deltas, ncol = expected_ncols)
  } else {
    deltas <- matrix(delta_stream, nrow = expected_nrows_deltas, ncol = expected_ncols)
  }

  cums <- matrix(apply(deltas, 2, cumsum), nrow = expected_nrows_deltas, ncol = expected_ncols)
  recon <- rbind(first_vals, sweep(cums, 2, first_vals, "+"))

  perm <- c(axis, setdiff(seq_along(dims), axis))
  recon_perm <- array(recon, dim = dims[perm])
  if (length(dims) > 1) {
    out <- aperm(recon_perm, order(perm))
  } else {
    out <- as.vector(recon_perm)
  }

  subset <- handle$subset
  if (!is.null(subset$roi_mask)) {
    vox_idx <- which(as.logical(subset$roi_mask))
    out <- out[vox_idx, , drop = FALSE]
  }
  if (!is.null(subset$time_idx)) {
    idx <- as.integer(subset$time_idx)
    if (axis == length(dims)) {
      out <- out[,, idx, drop = FALSE]
    } else {
      # for simplicity, convert to array and subset along axis
      ind <- vector("list", length(dims))
      for (i in seq_along(ind)) ind[[i]] <- seq_len(dim(out)[i])
      ind[[axis]] <- idx
      out <- do.call(`[`, c(list(out), ind, drop = FALSE))
    }
  }

  input_key <- desc$inputs[[1]] %||% "input"
  handle$update_stash(keys = character(), new_values = setNames(list(out), input_key))
}
