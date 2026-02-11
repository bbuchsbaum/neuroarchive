#' Simple Transform Builder API
#'
#' @description
#' Simplified API for creating custom transforms without needing to understand
#' the internal DataHandle/Plan/HDF5 machinery. Provides high-level functions
#' that handle all the boilerplate.
#'
#' @name simple_transform
NULL

#' Create a simple matrix-to-matrix transform
#'
#' @description
#' The easiest way to create a custom compression transform. You only need to
#' provide two functions: one that compresses a matrix, and one that decompresses it.
#' All the HDF5, schema validation, and pipeline integration is handled automatically.
#'
#' @param name Character string, name of your transform (e.g., "robustpca")
#' @param compress_fn Function that takes a matrix and returns a list with:
#'   \itemize{
#'     \item \code{compressed}: The compressed representation (matrix or list of matrices)
#'     \item \code{metadata}: Any metadata needed for decompression (list)
#'   }
#' @param decompress_fn Function that takes compressed data and metadata and returns
#'   the reconstructed matrix
#' @param schema List defining parameters (simplified format, will be converted to JSON schema)
#' @param default_params Named list of default parameter values
#' @param description Character string describing what this transform does
#'
#' @return Invisible NULL. Registers the transform methods with the package.
#'
#' @examples
#' \dontrun{
#' # Create a simple "top-k" compressor that keeps only the k largest values
#' simple_matrix_transform(
#'   name = "topk",
#'
#'   compress_fn = function(X, params) {
#'     k <- params$k %||% 100
#'
#'     # Find top-k values per column
#'     top_indices <- apply(abs(X), 2, function(col) {
#'       order(abs(col), decreasing = TRUE)[1:min(k, length(col))]
#'     })
#'
#'     # Store as sparse representation
#'     compressed <- list(
#'       indices = top_indices,
#'       values = sapply(1:ncol(X), function(i) X[top_indices[,i], i])
#'     )
#'
#'     metadata <- list(
#'       original_dims = dim(X),
#'       k_actual = k
#'     )
#'
#'     list(compressed = compressed, metadata = metadata)
#'   },
#'
#'   decompress_fn = function(compressed, metadata, params) {
#'     dims <- metadata$original_dims
#'     X_reconstructed <- matrix(0, nrow = dims[1], ncol = dims[2])
#'
#'     for (i in 1:ncol(X_reconstructed)) {
#'       indices <- compressed$indices[, i]
#'       values <- compressed$values[, i]
#'       X_reconstructed[indices, i] <- values
#'     }
#'
#'     X_reconstructed
#'   },
#'
#'   schema = list(
#'     k = list(type = "integer", default = 100, minimum = 1,
#'              description = "Number of top values to keep per column")
#'   ),
#'
#'   description = "Keep only the k largest magnitude values per column"
#' )
#'
#' # Now use it!
#' write_lna(my_data, "compressed.h5", transforms = "topk",
#'           transform_params = list(topk = list(k = 50)))
#' }
#'
#' @export
simple_matrix_transform <- function(name,
                                   compress_fn,
                                   decompress_fn,
                                   schema = list(),
                                   default_params = list(),
                                   description = "") {

  # Validate inputs
  stopifnot(
    is.character(name), length(name) == 1,
    is.function(compress_fn),
    is.function(decompress_fn),
    is.list(schema),
    is.list(default_params)
  )

  # Create JSON schema
  schema_path <- system.file("schemas", sprintf("%s.schema.json", name),
                            package = "neuroarchive")
  if (!nzchar(schema_path)) {
    # Generate schema from simplified format
    json_schema <- create_simple_schema(name, schema, default_params)

    # Try to write to inst/schemas if in development
    schema_dir <- file.path("inst", "schemas")
    if (dir.exists(schema_dir)) {
      schema_file <- file.path(schema_dir, sprintf("%s.schema.json", name))
      writeLines(json_schema, schema_file)
    }

    # Register schema in runtime cache for validation/defaults in installed environments
    try({
      schema_obj <- jsonlite::fromJSON(json_schema, simplifyVector = FALSE)
      assign(name, schema_obj, envir = .schema_cache)
    }, silent = TRUE)
  }

  # Create forward_step method
  forward_fn <- function(type, desc, handle) {
    # Get parameters
    params <- desc$params %||% default_params

    # Get input matrix using helper
    X <- .get_input_matrix(handle, desc)

    # Apply compression
    result <- compress_fn(X, params)
    compressed <- result$compressed
    metadata <- result$metadata

    # Store using helper
    .store_compressed_simple(
      handle = handle,
      desc = desc,
      compressed = compressed,
      metadata = metadata,
      transform_name = name
    )
  }

  # Create invert_step method
  inverse_fn <- function(type, desc, handle) {
    # Load compressed data using helper
    data_list <- .load_compressed_simple(handle, desc, name)

    # Apply decompression
    params <- desc$params %||% default_params
    reconstructed <- decompress_fn(
      compressed = data_list$compressed,
      metadata = data_list$metadata,
      params = params
    )

    # Handle subsetting
    reconstructed <- .apply_subsetting(reconstructed, handle$subset)

    # Store result
    input_key <- desc$inputs[[1]]
    handle$update_stash(keys = character(), new_values = setNames(list(reconstructed), input_key))
  }

  # Register methods without writing into the (locked) package namespace.
  .register_runtime_transform_method("forward_step", name, forward_fn)
  .register_runtime_transform_method("invert_step", name, inverse_fn)
  .register_runtime_transform_method("lna_default", name, function() default_params)

  # Register with transform registry
  register_transform(
    name,
    metadata = list(
      description = description,
      category = "custom",
      simplified_api = TRUE
    ),
    validate = FALSE
  )

  invisible(NULL)
}

#' Create a simple function-based transform
#'
#' @description
#' Even simpler than \code{simple_matrix_transform} - just provide a compression
#' function and its inverse. Perfect for simple operations like centering, scaling,
#' or filtering that don't need to store extra metadata.
#'
#' @param name Character string, transform name
#' @param forward_fn Function(X, params) that returns transformed X
#' @param inverse_fn Function(X, params) that returns original X
#' @param params List of parameter definitions (name = default_value)
#' @param description Character string describing the transform
#'
#' @examples
#' \dontrun{
#' # Create a simple robust centering transform
#' simple_function_transform(
#'   name = "robust_center",
#'
#'   forward_fn = function(X, params) {
#'     method <- params$method %||% "median"
#'     centers <- switch(method,
#'       median = apply(X, 2, median, na.rm = TRUE),
#'       mean = colMeans(X, na.rm = TRUE)
#'     )
#'     sweep(X, 2, centers, "-")
#'   },
#'
#'   inverse_fn = function(X, params) {
#'     # Centering is not reversible without storing centers!
#'     # This is why simple_matrix_transform is usually better
#'     warning("Cannot reverse robust_center without stored centers")
#'     X
#'   },
#'
#'   params = list(
#'     method = "median"
#'   ),
#'
#'   description = "Center each column using robust statistics"
#' )
#' }
#'
#' @export
simple_function_transform <- function(name,
                                     forward_fn,
                                     inverse_fn,
                                     params = list(),
                                     description = "") {

  # Wrapper that stores/retrieves minimal metadata
  compress_fn <- function(X, params) {
    transformed <- forward_fn(X, params)
    list(
      compressed = transformed,
      metadata = list(original_dims = dim(X))
    )
  }

  decompress_fn <- function(compressed, metadata, params) {
    inverse_fn(compressed, params)
  }

  # Convert params to schema format
  schema <- lapply(names(params), function(nm) {
    val <- params[[nm]]
    list(
      type = if (is.numeric(val)) "number" else "string",
      default = val
    )
  })
  names(schema) <- names(params)

  simple_matrix_transform(
    name = name,
    compress_fn = compress_fn,
    decompress_fn = decompress_fn,
    schema = schema,
    default_params = params,
    description = description
  )
}

#' Quick transform for dimensionality reduction
#'
#' @description
#' Specialized helper for PCA-like transforms. Handles the common pattern
#' of computing a basis and projecting data onto it.
#'
#' @param name Transform name
#' @param basis_fn Function(X, params) that returns basis matrix (features x components)
#' @param params Parameter definitions
#' @param description Transform description
#'
#' @examples
#' \dontrun{
#' simple_basis_transform(
#'   name = "robust_pca",
#'
#'   basis_fn = function(X, params) {
#'     k <- params$k %||% 50
#'
#'     # Median center
#'     centers <- apply(X, 2, median)
#'     X_centered <- sweep(X, 2, centers, "-")
#'
#'     # Robust weights using MAD
#'     mads <- apply(X_centered, 2, mad)
#'     weights <- 1 / (1 + abs(X_centered) / mads)
#'     X_weighted <- X_centered * sqrt(weights)
#'
#'     # SVD
#'     svd_result <- svd(X_weighted, nu = 0, nv = k)
#'
#'     # Return basis and metadata
#'     list(
#'       basis = svd_result$v,  # features x k
#'       centers = centers,
#'       mads = mads
#'     )
#'   },
#'
#'   params = list(k = 50),
#'   description = "Robust PCA using weighted covariance"
#' )
#' }
#'
#' @export
simple_basis_transform <- function(name, basis_fn, params = list(), description = "") {

  compress_fn <- function(X, params) {
    # Compute basis
    basis_result <- basis_fn(X, params)
    basis <- basis_result$basis

    # Extract metadata (everything except basis)
    metadata <- basis_result[names(basis_result) != "basis"]

    # Project onto basis
    coefficients <- X %*% basis  # time x k

    list(
      compressed = list(
        coefficients = coefficients,
        basis = basis
      ),
      metadata = c(metadata, list(original_dims = dim(X)))
    )
  }

  decompress_fn <- function(compressed, metadata, params) {
    # Reconstruct: X = coefficients * basis'
    coefficients <- compressed$coefficients
    basis <- compressed$basis

    reconstructed <- coefficients %*% t(basis)

    # Add back any stored offsets (like centers)
    if (!is.null(metadata$centers)) {
      reconstructed <- sweep(reconstructed, 2, metadata$centers, "+")
    }

    reconstructed
  }

  # Add k parameter to schema if not present
  if (!"k" %in% names(params)) {
    params$k <- 50
  }

  schema <- list(
    k = list(
      type = "integer",
      default = params$k,
      minimum = 1,
      description = "Number of basis components"
    )
  )

  simple_matrix_transform(
    name = name,
    compress_fn = compress_fn,
    decompress_fn = decompress_fn,
    schema = schema,
    default_params = params,
    description = description
  )
}

# Internal helpers --------------------------------------------------------

.get_input_matrix <- function(handle, desc) {
  input_key <- if (!is.null(desc$inputs) && length(desc$inputs) > 0) {
    desc$inputs[[1]]
  } else {
    inp <- handle$pull_first(c("input", "dense_mat", "aggregated_matrix"))
    inp$key
  }

  X <- handle$stash[[input_key]]
  if (is.null(X)) {
    abort_lna(
      sprintf("Input key '%s' not found", input_key),
      .subclass = "lna_error_contract"
    )
  }

  # Convert to 2D matrix
  X <- as_dense_mat(X)

  # Reshape if 4D
  if (length(dim(X)) == 4) {
    orig_dims <- dim(X)
    n_time <- orig_dims[4]
    n_voxels <- prod(orig_dims[1:3])
    X <- matrix(X, nrow = n_time, ncol = n_voxels)
  }

  X
}

.store_compressed_simple <- function(handle, desc, compressed, metadata, transform_name) {
  run_label <- handle$current_run_id
  step_idx <- handle$plan$next_index
  output_key <- desc$outputs[[1]]

  # Store compressed data in stash for next transform
  handle <- handle$update_stash(keys = character(), new_values = setNames(list(compressed), output_key))

  # Store metadata for reconstruction
  metadata_key <- sprintf("%s_metadata", transform_name)
  handle <- handle$update_stash(keys = character(), new_values = setNames(list(metadata), metadata_key))

  # Add compressed data to plan
  datasets_for_desc <- list()
  # Allow optional dtype hints from metadata$dtypes (named list)
  dtypes <- metadata$dtypes %||% list()
  if (is.list(compressed)) {
    for (nm in names(compressed)) {
      payload_key <- sprintf("%s_%s_%s_s%d", transform_name, nm, run_label, step_idx)
      handle$plan$add_payload(payload_key, compressed[[nm]])
      ds_path <- file.path("/scans", run_label, "transforms", sprintf("%s_s%d", transform_name, step_idx), nm)
      handle$plan$add_dataset_def(
        path = ds_path,
        role = "transform_data",
        producer = transform_name,
        origin = run_label,
        step_index = step_idx,
        params_json = jsonlite::toJSON(desc$params %||% list(), auto_unbox = TRUE),
        payload_key = payload_key,
        write_mode = "eager",
        dtype = as.character(dtypes[[nm]] %||% "float32")
      )
      datasets_for_desc[[length(datasets_for_desc) + 1L]] <- list(path = ds_path, role = "transform_data")
    }
  } else {
    # Single matrix/array
    payload_key <- sprintf("%s_data_%s_s%d", transform_name, run_label, step_idx)
    handle$plan$add_payload(payload_key, compressed)
    ds_path <- file.path("/scans", run_label, "transforms", sprintf("%s_s%d", transform_name, step_idx), "data")
    handle$plan$add_dataset_def(
      path = ds_path,
      role = "transform_data",
      producer = transform_name,
      origin = run_label,
      step_index = step_idx,
      params_json = jsonlite::toJSON(desc$params %||% list(), auto_unbox = TRUE),
      payload_key = payload_key,
      write_mode = "eager",
      dtype = as.character(dtypes[["data"]] %||% "float32")
    )
    datasets_for_desc[[length(datasets_for_desc) + 1L]] <- list(path = ds_path, role = "transform_data")
  }

  # Store metadata
  for (nm in names(metadata)) {
    if (!identical(nm, "dtypes")) { # do not persist dtype hints
      metadata_payload_key <- sprintf("%s_%s_%s_s%d", transform_name, nm, run_label, step_idx)
      handle$plan$add_payload(metadata_payload_key, metadata[[nm]])
      md_path <- file.path("/scans", run_label, "transforms", sprintf("%s_s%d", transform_name, step_idx), nm)
      handle$plan$add_dataset_def(
        path = md_path,
        role = "transform_metadata",
        producer = transform_name,
        origin = run_label,
        step_index = step_idx,
        params_json = jsonlite::toJSON(desc$params %||% list(), auto_unbox = TRUE),
        payload_key = metadata_payload_key,
        write_mode = "eager",
        dtype = NA_character_
      )
      datasets_for_desc[[length(datasets_for_desc) + 1L]] <- list(path = md_path, role = "transform_metadata")
    }
  }

  # Emit a transform descriptor so readers can discover this step
  fname <- handle$plan$get_next_filename(transform_name)
  desc_for_plan <- desc
  desc_for_plan$type <- transform_name
  desc_for_plan$datasets <- datasets_for_desc
  if (is.null(desc_for_plan$version)) desc_for_plan$version <- "1.0"
  if (is.null(desc_for_plan$outputs)) desc_for_plan$outputs <- list(output_key)
  handle$plan$add_descriptor(fname, desc_for_plan)

  handle
}

.load_compressed_simple <- function(handle, desc, transform_name) {
  h5 <- handle$h5
  datasets <- desc$datasets

  if (is.null(datasets) || length(datasets) == 0) {
    abort_lna(
      sprintf("No datasets found for %s", transform_name),
      .subclass = "lna_error_contract"
    )
  }

  # Separate compressed data from metadata
  compressed_data <- list()
  metadata <- list()

  for (ds in datasets) {
    path <- ds$path
    data <- h5[[path]]$read()

    # Determine if this is compressed data or metadata
    basename <- basename(path)
    if (ds$role == "transform_data") {
      compressed_data[[basename]] <- data
    } else if (ds$role == "transform_metadata") {
      metadata[[basename]] <- data
    }
  }

  # If single compressed item, unwrap from list
  if (length(compressed_data) == 1) {
    compressed_data <- compressed_data[[1]]
  }

  list(compressed = compressed_data, metadata = metadata)
}

.apply_subsetting <- function(X, subset) {
  if (is.null(subset)) return(X)

  if (!is.null(subset$time_idx)) {
    X <- X[subset$time_idx, , drop = FALSE]
  }
  if (!is.null(subset$roi_mask)) {
    roi <- subset$roi_mask
    if (is.array(roi)) roi <- as.vector(roi)
    roi <- as.logical(roi)
    if (length(roi) == ncol(X)) {
      X <- X[, roi, drop = FALSE]
    } else {
      warn_lna("ROI mask length does not match number of columns; skipping ROI subsetting",
               .subclass = "lna_warning_validation",
               location = "simple_transform:.apply_subsetting")
    }
  }

  X
}

create_simple_schema <- function(name, schema_def, defaults) {
  # Convert simplified schema to full JSON schema
  properties <- list()

  for (param_name in names(schema_def)) {
    param_spec <- schema_def[[param_name]]

    # Build property definition
    prop <- list(type = param_spec$type %||% "number")

    if (!is.null(param_spec$minimum)) prop$minimum <- param_spec$minimum
    if (!is.null(param_spec$maximum)) prop$maximum <- param_spec$maximum
    if (!is.null(param_spec$enum)) prop$enum <- param_spec$enum
    if (!is.null(param_spec$description)) prop$description <- param_spec$description
    if (!is.null(param_spec$default)) {
      prop$default <- param_spec$default
    } else if (param_name %in% names(defaults)) {
      prop$default <- defaults[[param_name]]
    }

    properties[[param_name]] <- prop
  }

  schema <- list(
    `$schema` = "http://json-schema.org/draft-07/schema#",
    `$id` = sprintf("https://neuroarchive.org/schemas/lna/2.0/%s.schema.json", name),
    title = sprintf("Parameters for '%s' transform", name),
    type = "object",
    properties = list(
      type = list(type = "string"),
      version = list(type = "string"),
      inputs = list(oneOf = list(
        list(type = "array"),
        list(type = "string")
      )),
      outputs = list(oneOf = list(
        list(type = "array"),
        list(type = "string")
      )),
      datasets = list(type = "array"),
      params = list(
        type = "object",
        properties = properties,
        additionalProperties = FALSE
      )
    ),
    additionalProperties = TRUE
  )

  jsonlite::toJSON(schema, auto_unbox = TRUE, pretty = TRUE)
}
