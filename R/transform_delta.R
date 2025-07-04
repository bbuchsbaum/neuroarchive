#' Delta Transform - Forward Step
#'
#' Computes first-order differences along a specified axis and optionally
#' run-length encodes the result.
#' The name stored in `desc$outputs` (if supplied) controls the key used to
#' stash the resulting delta stream. If `desc$outputs` is `NULL`, the default
#' key `"delta_stream"` is used. An empty character vector results in no
#' stash update.
#' @keywords internal
forward_step.delta <- function(type, desc, handle) {
  # Get incoming params, or an empty list if NULL
  incoming_params <- desc$params %||% list()

  # Determine effective parameters, applying defaults
  # These are the parameters that will be used for processing AND stored in the descriptor
  p <- list()
  p$order <- incoming_params$order %||% 1L

  # Determine input key primarily from the descriptor and fetch the data
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
  if (is.null(dims)) dims <- c(length(x))

  # Use original dimensions if available, otherwise fall back to current dimensions
  orig_dims_attr <- attr(x, "lna.orig_dims")
  if (!is.null(orig_dims_attr)) {
    orig_dims <- orig_dims_attr
  } else {
    orig_dims <- dims
  }

  # Resolve actual axis for processing and storage in descriptor
  p$axis <- incoming_params$axis %||% -1L
  if (p$axis == -1L) {
    p$axis <- length(orig_dims)
  }
  p$axis <- as.integer(p$axis)

  p$reference_value_storage <- incoming_params$reference_value_storage %||% "first_value_verbatim"
  p$coding_method <- incoming_params$coding_method %||% "none"

  # Validate coding_method
  if (!p$coding_method %in% c("none", "rle")) {
    abort_lna(
      sprintf("Unsupported coding_method '%s'", p$coding_method),
      .subclass = "lna_error_validation",
      location = "forward_step.delta:coding_method"
    )
  }

  # Validate order
  if (!identical(p$order, 1L)) {
    abort_lna("only order=1 supported", .subclass = "lna_error_validation",
              location = "forward_step.delta:order")
  }

  # Store the true original dimensions in p
  p$orig_dims <- orig_dims

  # Ensure p$axis is valid for these original dims
  if (p$axis < 1 || p$axis > length(orig_dims)) {
    abort_lna("axis out of bounds",
              .subclass = "lna_error_validation",
              location = "forward_step.delta:axis")
  }

  # Processing: we need to handle axis based on original dimensions but work with current array
  # For now, work with current dimensions but store original dimensions for invert step
  perm <- c(p$axis, setdiff(seq_along(dims), p$axis))
  xp <- aperm(x, perm)
  
  # Calculate dimensions for processing: use original dims for axis calculation
  if (length(orig_dims) == 1) {
    # 1D input: axis must be 1, and we have only 1 "column" in reshaped view
    dim_xp_col <- 1L
  } else {
    # Multi-dimensional: calculate product of non-axis dimensions
    non_axis_dims <- orig_dims[-p$axis]
    dim_xp_col <- if (length(non_axis_dims) > 0) prod(non_axis_dims) else 1L
  }
  
  # Reshape based on current array dimensions but using original axis logic
  axis_len <- dims[p$axis]  # Use current array's axis length
  dim(xp) <- c(axis_len, dim_xp_col)

  if (axis_len == 0) {
    first_vals <- array(numeric(0), dim = c(0, dim_xp_col))
    deltas <- array(numeric(0), dim = c(0, dim_xp_col)) # For consistency, though deltas rows would be max(0, nrow-1)
  } else {
    first_vals <- xp[1, , drop = FALSE]
    # Ensure dim is 1xN even if only 1 col after prod(dims[-p$axis])
    dim(first_vals) <- c(1, dim_xp_col) 
    deltas <- xp[-1, , drop = FALSE] - xp[-nrow(xp), , drop = FALSE]
    if (nrow(xp) == 1) { # Special case: input has only 1 element along diff axis
        deltas <- array(numeric(0), dim = c(0, dim_xp_col))
    }
  }

  if (identical(p$coding_method, "rle")) {
    delta_stream <- .encode_rle(as.vector(deltas))
  } else {
    delta_stream <- deltas
  }

  run_id <- handle$current_run_id %||% "run-01"
  run_id <- sanitize_run_id(run_id)
  plan <- handle$plan
  fname <- plan$get_next_filename(type)
  base <- tools::file_path_sans_ext(fname)
  scans_root <- lna_options("paths.scans_root")[[1]]
  delta_path <- paste0(scans_root, run_id, "/deltas/", base, "/delta_stream")
  first_path <- paste0(scans_root, run_id, "/deltas/", base, "/first_values")
  step_index <- plan$next_index
  
  # params_json uses the fully populated p
  params_json <- as.character(jsonlite::toJSON(p, auto_unbox = TRUE))

  # Update desc with fully populated params and other fields
  desc$version <- "1.0"
  desc$inputs <- c(input_key)
  desc$params <- p # This ensures the descriptor stored has all resolved params

  output_key <- NULL
  if (!is.null(desc$outputs)) {
    if (length(desc$outputs) > 0) {
      output_key <- desc$outputs[[1]]
      desc$outputs <- output_key # Ensure it's a single string if multiple were passed
    } else {
      desc$outputs <- character() # Empty string means no output key to stash
    }
  } else {
    output_key <- "delta_stream" # Default output key
    desc$outputs <- output_key
  }
  
  ds_list <- list(list(path = delta_path, role = "delta_stream"))
  if (identical(p$reference_value_storage, "first_value_verbatim")) {
    ds_list[[length(ds_list) + 1]] <- list(path = first_path, role = "first_values")
  }
  desc$datasets <- ds_list

  plan$add_descriptor(fname, desc)
  
  plan$add_payload(delta_path, delta_stream)
  
  plan$add_dataset_def(delta_path, "delta_stream", as.character(type), run_id,
                       as.integer(step_index), params_json,
                       delta_path, "eager")
  if (identical(p$reference_value_storage, "first_value_verbatim")) {
    plan$add_payload(first_path, first_vals)
    plan$add_dataset_def(first_path, "first_values", as.character(type), run_id,
                         as.integer(step_index), params_json,
                         first_path, "eager", dtype = NA_character_)
  }

  handle$plan <- plan
  if (!is.null(output_key) && nzchar(output_key)) {
    handle$update_stash(
      keys = c(input_key),
      new_values = setNames(list(delta_stream), output_key)
    )
  } else {
    # If output_key is NULL or empty, effectively remove input_key from subsequent visibility
    # Or, ensure no accidental new value is stashed if output_key is empty string
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
  if (is.null(delta_stream)) {
    delta_stream <- numeric(0)
  }

  expected_ncols <- if (length(dims) == 1) 1 else prod(dims[-axis])


  expected_rows_first_vals <- if (dims[axis] == 0) 0L else 1L

  if (identical(ref_store, "first_value_verbatim")) {
    first_vals <- h5_read(root, first_path)
    if (is.null(first_vals)) {
      first_vals <- numeric(0)
    }
  } else {
    first_vals <- array(0, dim = c(expected_rows_first_vals, expected_ncols))
  }

  # Validate and ensure correct dimensions for first_vals
  if (!is.matrix(first_vals) || ncol(first_vals) != expected_ncols || nrow(first_vals) != expected_rows_first_vals) {
      if (length(first_vals) == (expected_rows_first_vals * expected_ncols)) {
          dim(first_vals) <- c(expected_rows_first_vals, expected_ncols)
      } else {
          abort_lna(
              sprintf(
                  "first_vals dimensions are incorrect. Expected %dx%d, got length %d or dims %s",
                  expected_rows_first_vals, expected_ncols, length(first_vals), paste(dim(first_vals), collapse="x")
              ),
              .subclass = "lna_error_runtime",
              location = "invert_step.delta:first_vals_dim_check"
          )
      }
  }

  expected_nrows_deltas <- max(0, dims[axis] - 1L)

  if (identical(coding, "rle")) {
    delta_vec <- .decode_rle(delta_stream,
                             expected_nrows_deltas * expected_ncols,
                             location = "invert_step.delta")
    deltas <- matrix(delta_vec, nrow = expected_nrows_deltas, ncol = expected_ncols)
  } else {
    deltas <- matrix(delta_stream, nrow = expected_nrows_deltas, ncol = expected_ncols)
  }

  subset <- handle$subset
  roi_mask <- subset$roi_mask
  time_idx <- subset$time_idx

  if (!is.null(roi_mask)) {
    vox_idx <- which(as.logical(roi_mask))
    first_vals <- first_vals[, vox_idx, drop = FALSE]
    deltas <- deltas[, vox_idx, drop = FALSE]
    dims[-axis] <- length(vox_idx)
  }

  if (!is.null(time_idx)) {
    original_idx <- as.integer(time_idx)
    valid_t_values <- original_idx[original_idx > 0 & original_idx <= (nrow(deltas) + 1)]

    if (length(valid_t_values) == 0) {
        recon <- matrix(numeric(0), nrow = 0, ncol = ncol(deltas))
        dims[axis] <- 0L
    } else {
        idx <- valid_t_values 
        max_idx <- max(idx) 

        if (max_idx > 1) {
          nrow_for_cums <- min(max_idx - 1, nrow(deltas))
          if (nrow_for_cums > 0) {
            cums <- .col_cumsums(deltas[seq_len(nrow_for_cums), , drop = FALSE])
          } else {
            cums <- matrix(numeric(0), nrow = 0, ncol = ncol(deltas))
          }
        } else { 
          cums <- matrix(numeric(0), nrow = 0, ncol = ncol(deltas))
        }

        recon <- matrix(0, nrow = length(idx), ncol = ncol(deltas))
        for (i in seq_along(idx)) {
          t <- idx[i]
          if (t == 1) {
            recon[i, ] <- first_vals
          } else {
            # The following check should ideally not be needed if valid_t_values logic is correct
            # but kept for safety during debugging this complex interaction.
            if ((t - 1) > nrow(cums) || (t-1) <= 0) {
                stop(sprintf("Internal error in delta invert: trying to access cums[%d,], but cums only has %d rows or t-1 is invalid. t=%d, valid_idx_max=%d, nrow_for_cums=%d, nrow_deltas=%d, p$orig_dims[axis]=%d",
                             t-1, nrow(cums), t, max_idx, 
                             if(exists("nrow_for_cums")) nrow_for_cums else NA, 
                             nrow(deltas), p$orig_dims[axis]))
            }
            recon[i, ] <- first_vals + cums[t - 1, ]
          }
        }
        dims[axis] <- length(idx) 
    }
  } else {
    cums <- .col_cumsums(deltas)
    recon <- rbind(first_vals, sweep(cums, 2, first_vals, "+"))
  }

  # Determine output format based on subsetting
  subset <- handle$subset
  roi_applied <- !is.null(subset$roi_mask)
  time_applied <- !is.null(subset$time_idx)
  
  if (roi_applied && time_applied) {
    # When both ROI and time subsetting are applied, return as matrix (voxels x time)
    out <- t(recon)
  } else if (roi_applied || time_applied) {
    # When only one type of subsetting is applied, maintain original structure but with reduced dimensions
    perm <- c(axis, setdiff(seq_along(dims), axis))
    recon_perm <- array(recon, dim = dims[perm])
    if (length(dims) > 1) {
      out <- aperm(recon_perm, order(perm))
    } else {
      out <- as.vector(recon_perm)
    }
  } else {
    # No subsetting applied, restore original array structure
    perm <- c(axis, setdiff(seq_along(dims), axis))
    recon_perm <- array(recon, dim = dims[perm])
    if (length(dims) > 1) {
      out <- aperm(recon_perm, order(perm))
    } else {
      out <- as.vector(recon_perm)
    }
  }

  input_key <- desc$inputs[[1]] %||% "input"
  handle$update_stash(keys = character(), new_values = setNames(list(out), input_key))
}
