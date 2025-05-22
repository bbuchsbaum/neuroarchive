#' Delta Transform - Forward Step
#'
#' Computes first-order differences along a specified axis and optionally
#' run-length encodes the result.
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

  # Determine input key based on previous transform's output
  if (handle$has_key("sparsepca_embedding")) {
    input_key <- "sparsepca_embedding"
  } else if (handle$has_key("aggregated_matrix")) {
    input_key <- "aggregated_matrix"
  } else {
    input_key <- if (!is.null(desc$inputs)) desc$inputs[[1]] else "input"
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
  deltas <- xp[-1, , drop = FALSE] - xp[-nrow(xp), , drop = FALSE]

  if (identical(coding, "rle")) {
    vec <- as.vector(deltas)
    r <- rle(vec)
    delta_stream <- cbind(lengths = r$lengths, values = r$values)
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
  desc$outputs <- "delta_stream"
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
                       delta_path, "eager")
  if (identical(ref_store, "first_value_verbatim")) {
    plan$add_payload(first_path, first_vals)
    plan$add_dataset_def(first_path, "first_values", as.character(type), run_id,
                         as.integer(step_index), params_json,
                         first_path, "eager")
  }

  handle$plan <- plan
  handle$update_stash(keys = c(input_key), 
                      new_values = list(delta_stream = delta_stream))
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
  if (identical(ref_store, "first_value_verbatim")) {
    first_vals <- h5_read(root, first_path)
  } else {
    first_vals <- array(0, dim = prod(dims[-axis]))
  }
  if (is.null(dim(first_vals))) {
    dim(first_vals) <- c(1, length(first_vals))
  } else if (length(dims) > 1 && length(dim(first_vals)) == length(dims) - 1) {
    # If original was >1D and first_vals is N-1 D, make it 1 x prod(other_dims)
    dim(first_vals) <- c(1, prod(dims[-axis]))
  } else if (length(dims) == 1 && length(dim(first_vals)) == 1) {
    # Original was 1D, first_vals is 1D (should be length 1)
    dim(first_vals) <- c(1,1)
  }

  # Calculate expected number of columns for the 2D representation of deltas/cums
  expected_ncols <- if(length(dims) == 1) 1 else prod(dims[-axis])
  # Number of rows for deltas is one less than the original dimension along the axis
  expected_nrows_deltas <- max(0, dims[axis]-1)

  if (identical(coding, "rle")) {
    delta_vec <- rep(delta_stream[,2], delta_stream[,1])
    if (expected_nrows_deltas > 0 && length(delta_vec) != expected_nrows_deltas * expected_ncols) {
        abort_lna(sprintf("RLE decoded data length (%d) mismatch. Expected %d elements for %dx%d deltas matrix.", 
                        length(delta_vec), expected_nrows_deltas * expected_ncols, 
                        expected_nrows_deltas, expected_ncols),
                  .subclass = "lna_error_runtime", location="invert_step.delta:rle_decode")
    }
    deltas <- matrix(delta_vec, nrow = expected_nrows_deltas, ncol = expected_ncols)
  } else {
    deltas <- delta_stream
    # Ensure deltas is reshaped correctly. Total elements must match.
    if (length(as.vector(deltas)) != expected_nrows_deltas * expected_ncols) {
        abort_lna(sprintf("Delta stream data length (%d) mismatch. Expected %d elements for %dx%d deltas matrix.",
                        length(as.vector(deltas)), expected_nrows_deltas * expected_ncols,
                        expected_nrows_deltas, expected_ncols),
                  .subclass = "lna_error_runtime", location="invert_step.delta:delta_stream_shape")
    }
    dim(deltas) <- c(expected_nrows_deltas, expected_ncols)
  }

  # Ensure first_vals also has expected_ncols for sweep compatibility
  if (ncol(first_vals) != expected_ncols) {
     # This should ideally not happen if forward step saved it correctly (1 x prod(dims[-axis]))
     # and prod(dims[-axis]) is expected_ncols.
     dim(first_vals) <- c(1, expected_ncols) # Attempt to fix, assuming length(first_vals) == expected_ncols
  }

  cums <- apply(deltas, 2, cumsum)
  # apply on 0-row matrix results in 0-row matrix with same ncols
  # apply on 0-col matrix results in list() if simplify=T, or error if simplify=F 
  # expected_ncols cannot be 0 if dims is not empty.
  
  # Ensure cums has correct dimensions before sweep, esp. if deltas was empty
  if (nrow(deltas) == 0 && ncol(deltas) > 0) { # If deltas was 0xM
      dim(cums) <- c(0, ncol(deltas)) # apply might return list() or different shape
  } else if (nrow(deltas) > 0 && ncol(deltas) == 0) {
      # This case (Mx0 deltas) should not occur with current logic
      dim(cums) <- c(nrow(deltas),0)
  }
  # If deltas was 0x0, cums will be 0x0 (from matrix(...,0,0))

  recon <- rbind(first_vals, sweep(cums, 2, first_vals, "+"))

  # Reshape 'recon' back to original permuted dimensions
  if (length(dims) == 1) {
    # Original data was a vector. recon is N x 1 (or 1 x N if first_vals was row).
    # first_vals dim is c(1, prod(dims[-axis])), but if dims is 1D, prod(dims[-axis]) is prod(empty)=1.
    # So first_vals is 1x1. deltas is (N-1)x1. cums is (N-1)x1. recon is Nx1.
    out <- as.vector(recon)
  } else {
    # Original data was multi-dimensional array.
    # Target shape before inverse permutation is (dim_along_diff_axis, other_dim1, other_dim2, ...)
    # These are: retrieved_orig_dims[effective_axis], followed by retrieved_orig_dims[-effective_axis]
    
    current_shape_of_recon <- c(dims[axis], prod(dims[-axis]))
    # Ensure recon actually has this 2D shape. It should from rbind/sweep.
    dim(recon) <- current_shape_of_recon

    # The shape for aperm needs to match the permuted version from forward step
    # which was (dim_of_diff_axis, other_dims_collapsed_then_expanded...)
    # The actual dimensions for array() call prior to aperm should be:
    # retrieved_orig_dims but with effective_axis first.
    # e.g., if orig_dims = c(D1, D2, D3) and effective_axis = 2,
    # then permuted_dims_for_array = c(D2, D1, D3)
    
    permuted_dims_order <- c(axis, setdiff(seq_along(dims), axis))
    shape_before_aperm <- dims[permuted_dims_order]
    
    # Total elements in recon must match total elements in shape_before_aperm
    # numel_recon = prod(current_shape_of_recon)
    # numel_target = prod(shape_before_aperm)
    # These must be equal. They are, as both are prod(dims).
    
    dim(recon) <- shape_before_aperm # Reshape recon to this permuted multi-dimensional form

    inv_perm <- match(seq_along(dims), permuted_dims_order)
    out <- aperm(recon, inv_perm)
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
