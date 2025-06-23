#' Forward step for aggregate_runs transform
#'
#' This function aggregates multiple runs into a single matrix. It expects
#' the initial list of run objects to be available in the stash.
#'
#' @param type The transform type string.
#' @param desc The transform descriptor list.
#' @param handle The DataHandle object.
#' @return The updated DataHandle object with aggregated data in the stash.
#' @keywords internal
forward_step.aggregate_runs <- function(type, desc, handle) {
  # Get parameters with defaults
  params <- desc$params %||% list()
  agg_op <- params$agg_op %||% "sum"
  
  # Check if we have the initial input list in the stash
  # This should be set by core_write when aggregate_runs is the first transform
  input_key <- "initial_input_list"
  if (!input_key %in% names(handle$stash)) {
    # Try alternative: maybe we have individual run data
    input_key <- "input"
  }
  
  input_data <- handle$stash[[input_key]]
  
  if (is.null(input_data)) {
    abort_lna(
      "No input data found for aggregation. Expected 'initial_input_list' or 'input' in stash.",
      .subclass = "lna_error_missing_data",
      location = "forward_step.aggregate_runs"
    )
  }
  
  # Handle different input formats
  if (is.list(input_data) && !is.matrix(input_data) && !is.array(input_data)) {
    # List of runs - aggregate them
    matrices <- lapply(input_data, function(x) {
      if (is.matrix(x)) {
        x
      } else if (is.array(x)) {
        # Convert to matrix (time x voxels)
        dim_x <- dim(x)
        if (length(dim_x) == 4) {
          # 4D array: reshape to matrix
          matrix(x, nrow = dim_x[4], ncol = prod(dim_x[1:3]))
        } else if (length(dim_x) == 3) {
          # 3D array: single time point
          matrix(x, nrow = 1, ncol = prod(dim_x))
        } else {
          abort_lna(
            sprintf("Unsupported array dimensions: %s", paste(dim_x, collapse = "x")),
            .subclass = "lna_error_invalid_input",
            location = "forward_step.aggregate_runs"
          )
        }
      } else {
        abort_lna(
          sprintf("Unsupported data type for run: %s", class(x)[1]),
          .subclass = "lna_error_invalid_input", 
          location = "forward_step.aggregate_runs"
        )
      }
    })
    
    # Perform aggregation
    if (agg_op == "sum") {
      aggregated <- Reduce("+", matrices)
    } else if (agg_op == "mean") {
      aggregated <- Reduce("+", matrices) / length(matrices)
    } else if (agg_op == "concatenate") {
      aggregated <- do.call(rbind, matrices)
    } else {
      abort_lna(
        sprintf("Unsupported aggregation operation: %s", agg_op),
        .subclass = "lna_error_invalid_params",
        location = "forward_step.aggregate_runs"
      )
    }
    
    # Record which runs were included
    run_ids <- names(input_data)
    if (is.null(run_ids)) {
      run_ids <- paste0("run-", seq_along(input_data))
    }
    desc$params$runs_included <- run_ids
    
  } else if (is.matrix(input_data) || is.array(input_data)) {
    # Single run - just pass through
    aggregated <- if (is.matrix(input_data)) input_data else {
      dim_x <- dim(input_data)
      if (length(dim_x) == 4) {
        matrix(input_data, nrow = dim_x[4], ncol = prod(dim_x[1:3]))
      } else {
        matrix(input_data, nrow = 1, ncol = prod(dim_x))
      }
    }
    desc$params$runs_included <- "single_run"
  } else {
    abort_lna(
      sprintf("Unsupported input data type: %s", class(input_data)[1]),
      .subclass = "lna_error_invalid_input",
      location = "forward_step.aggregate_runs"
    )
  }
  
  # Store original dimensions for inverse step
  desc$params$orig_dims <- paste(dim(aggregated), collapse = "x")
  desc$params$agg_op <- agg_op
  
  # Update descriptor
  desc$inputs <- c(input_key)
  desc$outputs <- c("aggregated_matrix")
  desc$version <- "1.0"
  
  # Update stash
  handle <- handle$update_stash(
    keys = desc$outputs,
    new_values = list(aggregated_matrix = aggregated)
  )
  
  # Add descriptor to plan
  handle$plan$add_descriptor(
    transform_name = handle$plan$get_next_filename(type),
    desc_list = desc
  )
  
  handle
}

#' Inverse step for aggregate_runs transform
#'
#' This function attempts to invert the aggregation operation. For operations
#' like "sum" or "mean", it passes through the aggregated data. 
#' The dimensions of the input data X_agg are expected to match desc$params$orig_dims.
#'
#' @param type The transform type string.
#' @param desc The transform descriptor list.
#' @param handle The DataHandle object.
#' @return The updated DataHandle object with the inverted data in the stash.
#' @keywords internal
invert_step.aggregate_runs <- function(type, desc, handle) {
  # Invert: the forward step's output key is now our conceptual input key from stash
  input_key <- desc$outputs[[1]] 
  X_agg <- handle$stash[[input_key]]

  if (is.null(X_agg)) {
    abort_lna(
      sprintf("Input data for key '%s' is NULL in invert_step.aggregate_runs.", input_key),
      .subclass = "lna_error_missing_data",
      location = "invert_step.aggregate_runs:input_X_agg_null"
    )
  }

  orig_dims_str <- desc$params$orig_dims 
  agg_op <- desc$params$agg_op %||% "sum" # Default to sum if not specified 
  # agg_dim <- desc$params$agg_dim # May need for mean scaling factor

  if (is.null(orig_dims_str)) {
    abort_lna(
      "Required parameter 'orig_dims' missing from descriptor.",
      .subclass = "lna_error_descriptor",
      location = "invert_step.aggregate_runs:params_orig_dims_missing"
    )
  }
  orig_dims_vec <- as.integer(strsplit(orig_dims_str, "x")[[1]])

  # Ensure X_agg has the dimensions recorded by the forward step's aggregator
  if (!identical(as.integer(dim(X_agg)), orig_dims_vec)) {
    # If X_agg is a vector but its length matches prod(orig_dims_vec), try to reshape it.
    if (is.null(dim(X_agg)) && length(X_agg) == prod(orig_dims_vec)) {
      dim(X_agg) <- orig_dims_vec
    } else {
      # Check if this might be a time subsetting case
      # If the first dimension is different but other dimensions match, this might be subsetting
      if (length(dim(X_agg)) == length(orig_dims_vec) && 
          length(orig_dims_vec) >= 1 && 
          all(dim(X_agg)[-1] == orig_dims_vec[-1])) {
        # This is likely time subsetting, proceed without error
      } else {
        abort_lna(
          sprintf("Dimension mismatch for X_agg in invert_step.aggregate_runs. Expected %s (from desc$params$orig_dims), got %s. X_agg class: %s, input_key: %s",
                  orig_dims_str, paste(dim(X_agg), collapse="x"), class(X_agg)[1], input_key),
          .subclass = "lna_error_dimension_mismatch",
          location = "invert_step.aggregate_runs"
        )
      }
    }
  }
  
  X_inverted_orig_shape <- X_agg
  
  if (agg_op == "mean") {
    # This is a placeholder for mean. True inverse of mean requires knowing N_items.
    # For now, it doesn't change X_agg if it was a mean.
    # A proper inverse of mean would be X_agg * N_items_aggregated.
    # N_items_aggregated would be related to desc$params$agg_dim and original run structure.
    # This needs to be properly implemented if mean aggregation is used and needs inversion.
  }

  # The output of this inverse step should be what the forward step originally took as input key name
  output_key <- desc$inputs[[1]] %||% "input" 
  handle <- handle$update_stash(keys = names(handle$stash), # Clear all previous keys from stash
                               new_values = setNames(list(X_inverted_orig_shape), output_key))
  handle
}

#' Default parameters for aggregate_runs
#' @export
#' @keywords internal
lna_default.aggregate_runs <- function() {
  list(
    agg_op = "sum", 
    agg_dim = NULL, 
    orig_dims = NULL 
  )
} 