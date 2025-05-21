#' Delta Transform - Forward Step
#'
#' Computes first-order differences along a specified axis and optionally
#' run-length encodes the result.
#' @keywords internal
forward_step.delta <- function(type, desc, handle) {
  p <- desc$params %||% list()
  order <- p$order %||% 1
  axis <- p$axis %||% -1
  ref_store <- p$reference_value_storage %||% "first_value_verbatim"
  coding <- p$coding_method %||% "none"

  if (!identical(order, 1L)) {
    abort_lna("only order=1 supported", .subclass = "lna_error_validation",
              location = "forward_step.delta:order")
  }

  input_key <- if (!is.null(desc$inputs)) desc$inputs[[1]] else "input"
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
  deltas <- xp[-1, , drop = FALSE] - xp[-nrow(xp), , drop = FALSE]

  if (identical(coding, "rle")) {
    vec <- as.vector(deltas)
    r <- rle(vec)
    delta_stream <- cbind(lengths = r$lengths, values = r$values)
  } else {
    delta_stream <- deltas
  }

  run_id <- handle$current_run_id %||% "run-01"
  plan <- handle$plan
  fname <- plan$get_next_filename(type)
  base <- tools::file_path_sans_ext(fname)
  delta_path <- paste0("/scans/", run_id, "/deltas/", base, "/delta_stream")
  first_path <- paste0("/scans/", run_id, "/deltas/", base, "/first_values")
  step_index <- plan$next_index
  p$orig_dims <- dims
  params_json <- jsonlite::toJSON(p, auto_unbox = TRUE)

  desc$version <- "1.0"
  desc$inputs <- c(input_key)
  desc$outputs <- character()
  desc$params <- p
  ds <- list(list(path = delta_path, role = "delta_stream"))
  if (identical(ref_store, "first_value_verbatim")) {
    ds[[length(ds) + 1]] <- list(path = first_path, role = "first_values")
  }
  desc$datasets <- ds

  plan$add_descriptor(fname, desc)
  plan$add_payload(delta_path, delta_stream)
  plan$add_dataset_def(delta_path, "delta_stream", type, run_id,
                       as.integer(step_index), params_json,
                       delta_path, "eager")
  if (identical(ref_store, "first_value_verbatim")) {
    plan$add_payload(first_path, first_vals)
    plan$add_dataset_def(first_path, "first_values", type, run_id,
                         as.integer(step_index), params_json,
                         first_path, "eager")
  }

  handle$plan <- plan
  handle$update_stash(keys = input_key, new_values = list())
}

#' Delta Transform - Inverse Step
#'
#' Reconstructs data from delta representation stored in HDF5.
#' @keywords internal
invert_step.delta <- function(type, desc, handle) {
  p <- desc$params %||% list()
  axis <- p$axis %||% -1
  ref_store <- p$reference_value_storage %||% "first_value_verbatim"
  coding <- p$coding_method %||% "none"
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
  if (identical(ref_store, "first_value_verbatim")) {
    first_vals <- h5_read(root, first_path)
  } else {
    first_vals <- array(0, dim = prod(dims[-axis]))
  }
  if (is.null(dim(first_vals))) {
    dim(first_vals) <- c(1, length(first_vals))
  } else if (length(dim(first_vals)) == length(dims) - 1) {
    dim(first_vals) <- c(1, prod(dims[-axis]))
  }

  if (identical(coding, "rle")) {
    delta_vec <- rep(delta_stream[,2], delta_stream[,1])
    deltas <- matrix(delta_vec, nrow = dims[axis] - 1)
  } else {
    deltas <- delta_stream
  }

  if (ncol(deltas) != prod(dims[-axis])) {
    dim(deltas) <- c(dims[axis] - 1, prod(dims[-axis]))
  }

  cums <- apply(deltas, 2, cumsum)
  recon <- rbind(first_vals, sweep(cums, 2, first_vals, "+"))
  dim(recon) <- c(dims[axis], prod(dims[-axis]))
  aperm_recon <- array(recon, dim = c(dims[axis], dims[-axis]))
  inv_perm <- match(seq_along(dims), c(axis, setdiff(seq_along(dims), axis)))
  out <- aperm(aperm_recon, inv_perm)

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
