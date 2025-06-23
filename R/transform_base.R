#' Transform Base Utilities
#'
#' @description
#' Internal utilities for transform implementations to reduce code duplication
#' and ensure consistency across transform types.
#'
#' @keywords internal

#' Extract and validate transform parameters
#'
#' @param desc Transform descriptor containing params
#' @param param_spec Named list defining parameter specifications
#' @return Named list of validated parameters
#' @keywords internal
extract_transform_params <- function(desc, param_spec) {
  p <- desc$params %||% list()
  result <- list()
  
  for (param_name in names(param_spec)) {
    spec <- param_spec[[param_name]]
    value <- p[[param_name]] %||% spec$default
    
    # Type validation
    if (!is.null(spec$type)) {
      type_check <- switch(spec$type,
        "integer" = is.numeric(value) && all(value == as.integer(value)),
        "numeric" = is.numeric(value),
        "character" = is.character(value),
        "logical" = is.logical(value),
        TRUE
      )
      
      if (!type_check) {
        abort_lna(
          sprintf("Parameter '%s' must be of type %s", param_name, spec$type),
          .subclass = "lna_error_validation",
          location = sprintf("extract_transform_params:%s", param_name)
        )
      }
    }
    
    # Length validation
    if (!is.null(spec$length)) {
      if (length(value) != spec$length) {
        abort_lna(
          sprintf("Parameter '%s' must have length %d", param_name, spec$length),
          .subclass = "lna_error_validation",
          location = sprintf("extract_transform_params:%s", param_name)
        )
      }
    }
    
    # Range validation
    if (!is.null(spec$min) || !is.null(spec$max)) {
      if (!is.null(spec$min) && any(value < spec$min)) {
        abort_lna(
          sprintf("Parameter '%s' must be >= %s", param_name, spec$min),
          .subclass = "lna_error_validation",
          location = sprintf("extract_transform_params:%s", param_name)
        )
      }
      if (!is.null(spec$max) && any(value > spec$max)) {
        abort_lna(
          sprintf("Parameter '%s' must be <= %s", param_name, spec$max),
          .subclass = "lna_error_validation",
          location = sprintf("extract_transform_params:%s", param_name)
        )
      }
    }
    
    # Enum validation
    if (!is.null(spec$enum)) {
      if (!all(value %in% spec$enum)) {
        abort_lna(
          sprintf("Parameter '%s' must be one of: %s", 
                  param_name, paste(spec$enum, collapse = ", ")),
          .subclass = "lna_error_validation",
          location = sprintf("extract_transform_params:%s", param_name)
        )
      }
    }
    
    # Custom validator
    if (!is.null(spec$validator) && is.function(spec$validator)) {
      validation_result <- spec$validator(value, param_name)
      if (!isTRUE(validation_result)) {
        error_msg <- if (is.character(validation_result)) {
          validation_result
        } else {
          sprintf("Parameter '%s' failed custom validation", param_name)
        }
        abort_lna(
          error_msg,
          .subclass = "lna_error_validation",
          location = sprintf("extract_transform_params:%s", param_name)
        )
      }
    }
    
    result[[param_name]] <- value
  }
  
  result
}

#' Create validator functions for common patterns
#'
#' @name validators
#' @keywords internal
NULL

#' @rdname validators
validate_integer_range <- function(min, max) {
  function(value, param_name) {
    if (!is.numeric(value) || length(value) != 1) {
      return(sprintf("'%s' must be a single numeric value", param_name))
    }
    if (value != as.integer(value)) {
      return(sprintf("'%s' must be an integer", param_name))
    }
    if (value < min || value > max) {
      return(sprintf("'%s' must be between %d and %d", param_name, min, max))
    }
    TRUE
  }
}

#' @rdname validators
validate_enum <- function(allowed_values) {
  function(value, param_name) {
    if (length(value) != 1) {
      return(sprintf("'%s' must be a single value", param_name))
    }
    if (!value %in% allowed_values) {
      return(sprintf("'%s' must be one of: %s", 
                     param_name, paste(allowed_values, collapse = ", ")))
    }
    TRUE
  }
}

#' @rdname validators
validate_probability <- function(value, param_name) {
  if (!is.numeric(value) || length(value) != 1) {
    return(sprintf("'%s' must be a single numeric value", param_name))
  }
  if (value < 0 || value > 1) {
    return(sprintf("'%s' must be between 0 and 1", param_name))
  }
  TRUE
}

#' Generate standard paths for transform outputs
#'
#' @param transform_type Transform type name
#' @param run_id Run identifier
#' @param base_name Base name for the transform
#' @param dataset_type Type of dataset (e.g., "quantized", "coefficients")
#' @return Character path string
#' @keywords internal
generate_transform_path <- function(transform_type, run_id, base_name = NULL, 
                                    dataset_type = NULL) {
  
  run_id <- sanitize_run_id(run_id)
  
  # Get root paths from configuration
  scans_root <- lna_options("paths.scans_root")[[1]]
  transforms_root <- lna_options("paths.transforms_root")[[1]]
  temporal_root <- lna_options("paths.temporal_root")[[1]]
  
  # Determine which root to use based on transform type and dataset
  if (transform_type == "temporal" && dataset_type == "basis") {
    root <- temporal_root
    path <- paste0(root, base_name, "/", dataset_type)
  } else if (dataset_type %in% c("report", "metadata")) {
    root <- transforms_root
    path <- paste0(root, base_name, "_", dataset_type, ".json")
  } else {
    root <- scans_root
    if (!is.null(base_name)) {
      path <- paste0(root, run_id, "/", base_name, "/", dataset_type)
    } else {
      path <- paste0(root, run_id, "/", dataset_type)
    }
  }
  
  path
}

#' Resolve input key for transform
#'
#' @param handle DataHandle object
#' @param desc Transform descriptor
#' @param preferred_keys Character vector of preferred keys to check first
#' @return Character string of resolved input key
#' @keywords internal
resolve_input_key <- function(handle, desc, preferred_keys = character()) {
  # Check preferred keys first
  for (key in preferred_keys) {
    if (handle$has_key(key)) {
      return(key)
    }
  }
  
  # Check descriptor inputs
  if (!is.null(desc$inputs) && length(desc$inputs) > 0) {
    # If descriptor specifies inputs, use the first available one
    for (input in desc$inputs) {
      if (handle$has_key(input)) {
        return(input)
      }
    }
    # If none exist, still use the first specified
    return(desc$inputs[[1]])
  }
  
  # Ultimate fallback
  lna_options("defaults.input_key")[[1]] %||% "input"
}

#' Apply subset to data
#'
#' @param data Data array or matrix
#' @param subset Subset specification from handle
#' @param data_type Type of data structure ("matrix" or "array")
#' @return Subsetted data
#' @keywords internal
apply_subset <- function(data, subset, data_type = c("matrix", "array")) {
  if (is.null(subset)) {
    return(data)
  }
  
  data_type <- match.arg(data_type)
  
  # Apply ROI mask
  if (!is.null(subset$roi_mask)) {
    roi <- as.logical(subset$roi_mask)
    
    if (data_type == "matrix") {
      data <- data[, roi, drop = FALSE]
    } else if (data_type == "array") {
      # For 4D arrays, ROI applies to spatial dimensions
      dims <- dim(data)
      if (length(dims) == 4) {
        # Reshape to matrix, apply ROI, reshape back
        mat <- matrix(data, prod(dims[1:3]), dims[4])
        mat <- mat[roi, , drop = FALSE]
        # Can't reshape back to original spatial dims, keep as matrix
        data <- mat
      }
    }
  }
  
  # Apply time subset
  if (!is.null(subset$time_idx)) {
    idx <- as.integer(subset$time_idx)
    
    if (data_type == "matrix") {
      data <- data[idx, , drop = FALSE]
    } else if (data_type == "array") {
      dims <- dim(data)
      if (length(dims) == 4) {
        data <- data[, , , idx, drop = FALSE]
      } else if (length(dims) == 2) {
        # Already reshaped by ROI
        data <- data[, idx, drop = FALSE]
      }
    }
  }
  
  data
}

#' Execute transform with consistent error handling
#'
#' @param expr Expression to execute
#' @param type Transform type
#' @param step Step name (e.g., "forward", "invert")
#' @param handle DataHandle object
#' @return Result of expression
#' @keywords internal
with_transform_error_handling <- function(expr, type, step, handle = NULL) {
  tryCatch(
    expr,
    lna_error_validation = function(e) {
      # Add transform context to validation errors
      e$transform_type <- type
      e$transform_step <- step
      
      # Log validation errors if handle available
      if (!is.null(handle) && !is.null(handle$meta)) {
        handle$meta$validation_errors <- c(
          handle$meta$validation_errors,
          list(list(
            type = type,
            step = step,
            error = conditionMessage(e),
            timestamp = Sys.time()
          ))
        )
      }
      rlang::cnd_signal(e)  # Re-signal
    },
    error = function(e) {
      # Wrap generic errors with transform context
      if (!inherits(e, "lna_error")) {
        abort_lna(
          sprintf("Transform '%s' failed at %s: %s", type, step, conditionMessage(e)),
          .subclass = "lna_error_transform",
          parent = e,
          location = sprintf("%s:%s", type, step),
          transform_type = type,
          transform_step = step
        )
      } else {
        # Already an lna_error, just re-signal
        rlang::cnd_signal(e)
      }
    }
  )
}