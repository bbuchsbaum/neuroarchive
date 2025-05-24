#' Inverse step for myorg.aggregate_runs transform
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
invert_step.myorg.aggregate_runs <- function(type, desc, handle) {
  # Invert: the forward step's output key is now our conceptual input key from stash
  input_key <- desc$outputs[[1]] 
  X_agg <- handle$stash[[input_key]]

  if (is.null(X_agg)) {
    abort_lna(
      sprintf("Input data for key '%s' is NULL in invert_step.myorg.aggregate_runs.", input_key),
      .subclass = "lna_error_missing_data",
      location = "invert_step.myorg.aggregate_runs:input_X_agg_null"
    )
  }

  orig_dims_str <- desc$params$orig_dims 
  agg_op <- desc$params$agg_op %||% "sum" # Default to sum if not specified 
  # agg_dim <- desc$params$agg_dim # May need for mean scaling factor

  if (is.null(orig_dims_str)) {
    abort_lna(
      "Required parameter 'orig_dims' missing from descriptor.",
      .subclass = "lna_error_descriptor",
      location = "invert_step.myorg.aggregate_runs:params_orig_dims_missing"
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
        stop(sprintf("Dimension mismatch for X_agg in invert_step.aggregate_runs. Expected %s (from desc$params$orig_dims), got %s. X_agg class: %s, input_key: %s",
                     orig_dims_str, paste(dim(X_agg), collapse="x"), class(X_agg)[1], input_key))
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

#' Default parameters for myorg.aggregate_runs
#' @export
#' @keywords internal
lna_default.myorg.aggregate_runs <- function() {
  list(
    agg_op = "sum", 
    agg_dim = NULL, 
    orig_dims = NULL 
  )
} 