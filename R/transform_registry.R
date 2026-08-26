#' Transform Registry System
#'
#' @description
#' Internal system for registering and discovering transform implementations.
#' The registry provides metadata about available transforms, their capabilities,
#' and requirements.
#'
#' @keywords internal

# Create registry environment
.transform_registry <- new.env(parent = emptyenv())

#' Register a transform
#'
#' @param type Character transform type name
#' @param metadata List of metadata about the transform
#' @param validate Whether to validate the transform has required methods
#' @return Invisible NULL
#' @keywords internal
register_transform <- function(type, metadata = list(), validate = TRUE) {
  if (!is.character(type) || length(type) != 1) {
    abort_lna(
      "Transform type must be a single character string",
      .subclass = "lna_error_validation",
      location = "register_transform"
    )
  }
  
  # Check if already registered
  if (type %in% names(.transform_registry)) {
    warn_lna(
      sprintf("Transform '%s' already registered, overwriting", type),
      .subclass = "lna_warning_overwrite",
      location = "register_transform"
    )
  }
  
  # Build transform info
  info <- list(
    type = type,
    registered_at = Sys.time()
  )
  
  # Check for S3 methods
  if (validate) {
    info$has_forward <- !is.null(utils::getS3method("forward_step", type, optional = TRUE))
    info$has_invert <- !is.null(utils::getS3method("invert_step", type, optional = TRUE))
    
    if (!info$has_forward) {
      warn_lna(
        sprintf("Transform '%s' has no forward_step method", type),
        .subclass = "lna_warning_missing_method",
        location = "register_transform"
      )
    }
  } else {
    # Skip validation, assume methods exist
    info$has_forward <- TRUE
    info$has_invert <- TRUE
  }
  
  # Check for schema
  schema_file <- system.file("schemas", paste0(type, ".schema.json"), 
                            package = "neuroarchive")
  info$has_schema <- nzchar(schema_file) && file.exists(schema_file)
  
  # Check for additional methods
  info$has_validate_params <- !is.null(utils::getS3method("validate_params", type, optional = TRUE))
  info$has_prepare <- !is.null(utils::getS3method("prepare_transform", type, optional = TRUE))
  info$has_cleanup <- !is.null(utils::getS3method("cleanup_transform", type, optional = TRUE))
  
  # Get minimum dimensions
  info$min_dims <- tryCatch(
    transform_min_dims(type),
    error = function(e) 3L  # Default
  )
  
  # Check for special capabilities based on type patterns
  info$capabilities <- list()
  
  # Temporal transforms have special methods
  if (type == "temporal") {
    info$capabilities$temporal_basis <- !is.null(utils::getS3method("temporal_basis", "default", optional = TRUE))
    info$capabilities$temporal_project <- !is.null(utils::getS3method("temporal_project", "default", optional = TRUE))
    info$capabilities$temporal_reconstruct <- !is.null(utils::getS3method("temporal_reconstruct", "default", optional = TRUE))
  }
  
  # Basis transforms can have subtypes
  if (type == "basis" || grepl("basis", type)) {
    info$capabilities$supports_pca <- TRUE
    info$capabilities$supports_empirical = grepl("empirical", type)
    info$capabilities$supports_hrbf = grepl("hrbf", type)
  }
  
  # Compression transforms
  if (type %in% c("quant", "delta") || grepl("quant|delta", type)) {
    info$capabilities$compression = TRUE
    info$capabilities$lossy = (type == "quant" || grepl("quant", type))
    info$capabilities$lossless = (type == "delta" || grepl("delta", type))
  }
  
  # Spatial transforms
  if (grepl("spat|hrbf", type)) {
    info$capabilities$spatial = TRUE
  }
  
  # Add user metadata
  info$metadata <- metadata

  subset_caps <- .transform_supports_subset(type)
  info$lifecycle <- .transform_lifecycle_status(type, metadata)
  info$roundtrip <- isTRUE(info$has_forward) && isTRUE(info$has_invert)
  info$supports_roi_mask <- isTRUE(subset_caps$roi_mask)
  info$supports_time_idx <- isTRUE(subset_caps$time_idx)
  info$supports_parcel_idx <- isTRUE(subset_caps$parcel_idx)
  info$physical_roi_mask <- isTRUE(subset_caps$physical$roi_mask)
  info$physical_time_idx <- isTRUE(subset_caps$physical$time_idx)
  info$physical_parcel_idx <- isTRUE(subset_caps$physical$parcel_idx)
  info$temporal_granularity <- subset_caps$temporal_granularity
  info$parcel_granularity <- subset_caps$parcel_granularity
  info$roi_strategy <- subset_caps$roi_strategy
  info$capabilities$roundtrip <- info$roundtrip
  info$capabilities$supports_roi_mask <- info$supports_roi_mask
  info$capabilities$supports_time_idx <- info$supports_time_idx
  info$capabilities$supports_parcel_idx <- info$supports_parcel_idx
  info$capabilities$physical_roi_mask <- info$physical_roi_mask
  info$capabilities$physical_time_idx <- info$physical_time_idx
  info$capabilities$physical_parcel_idx <- info$physical_parcel_idx
  info$capabilities$temporal_granularity <- info$temporal_granularity
  info$capabilities$parcel_granularity <- info$parcel_granularity
  info$capabilities$roi_strategy <- info$roi_strategy
  
  # Store in registry
  .transform_registry[[type]] <- info
  
  invisible(NULL)
}

#' List registered transforms
#'
#' @param category Optional category filter
#' @param capability Optional capability filter
#' @param with_details Include detailed information
#' @return List or data frame of transforms
#' @export
list_transforms <- function(category = NULL, capability = NULL, with_details = FALSE) {
  transforms <- as.list(.transform_registry)
  
  # Filter by category if provided
  if (!is.null(category)) {
    transforms <- Filter(function(t) {
      isTRUE(t$metadata$category == category)
    }, transforms)
  }
  
  # Filter by capability if provided
  if (!is.null(capability)) {
    transforms <- Filter(function(t) {
      capability %in% names(t$capabilities) && isTRUE(t$capabilities[[capability]])
    }, transforms)
  }
  
  if (length(transforms) == 0) {
    return(if (with_details) list() else character())
  }
  
  if (with_details) {
    return(transforms)
  }
  
  # Simple list of names
  names(transforms)
}

#' Get transform information
#'
#' @param type Transform type
#' @return Transform info list or NULL if not found
#' @export
get_transform_info <- function(type) {
  if (!type %in% names(.transform_registry)) {
    return(NULL)
  }
  .transform_registry[[type]]
}

#' Check if transform is registered
#'
#' @param type Transform type
#' @return Logical
#' @export
is_transform_registered <- function(type) {
  type %in% names(.transform_registry)
}

#' Get transform capabilities
#'
#' @param type Transform type
#' @return Named list of capabilities
#' @export
get_transform_capabilities <- function(type) {
  info <- get_transform_info(type)
  if (is.null(info)) {
    return(list())
  }
  info$capabilities
}

#' Discover and register transforms
#'
#' Scans for transform implementations and registers them.
#' Called automatically on package load.
#'
#' @param pattern Regex pattern for finding transform files
#' @param validate Whether to validate each transform
#' @return Number of transforms registered
#' @keywords internal
discover_and_register_transforms <- function(pattern = "^transform_.*\\.R$", 
                                           validate = TRUE) {
  types <- character()

  # Prefer installed schema metadata because it is available both from source
  # and from an installed package during R CMD check.
  schema_dir <- system.file("schemas", package = "neuroarchive")
  if (nzchar(schema_dir) && dir.exists(schema_dir)) {
    schema_files <- list.files(schema_dir, pattern = "\\.schema\\.json$", full.names = FALSE)
    schema_types <- sub("\\.schema\\.json$", "", schema_files)
    types <- union(types, schema_types)
  }

  # Fall back to source-file discovery when running from a checkout.
  pkg_dir <- system.file("R", package = "neuroarchive")
  if (!nzchar(pkg_dir)) {
    pkg_dir <- "R"
  }
  if (dir.exists(pkg_dir)) {
    transform_files <- list.files(pkg_dir, pattern = pattern, full.names = FALSE)
    file_types <- gsub("^transform_", "", gsub("\\.R$", "", transform_files))
    types <- union(types, file_types)
  }

  if (length(types) == 0L) {
    return(0L)
  }
  
  # Filter out some special cases
  types <- setdiff(types, c("base", "builder", "registry", "meta"))
  
  # Register each transform
  count <- 0L
  for (type in types) {
    tryCatch({
      # Determine category based on type
      category <- if (type %in% c("quant", "delta")) {
        "compression"
      } else if (type %in% c("basis", "embed", "sparsepca") || grepl("basis", type)) {
        "dimensionality"
      } else if (type == "temporal" || grepl("temporal", type)) {
        "temporal"
      } else if (grepl("spat|hrbf", type)) {
        "spatial"
      } else if (type == "aggregate_runs") {
        "utility"
      } else {
        "other"
      }
      
      register_transform(type, list(category = category), validate = validate)
      count <- count + 1L
    }, error = function(e) {
      # Silently skip transforms that fail to register during discovery
      NULL
    })
  }
  
  count
}

#' Clear transform registry
#'
#' @return Invisible NULL
#' @keywords internal
clear_transform_registry <- function() {
  rm(list = ls(.transform_registry), envir = .transform_registry)
  invisible(NULL)
}

#' Get transform summary
#'
#' @return Data frame summarizing registered transforms
#' @export
transform_summary <- function() {
  transforms <- as.list(.transform_registry)
  
  if (length(transforms) == 0) {
    return(data.frame(
      type = character(),
      category = character(),
      lifecycle = character(),
      has_forward = logical(),
      has_invert = logical(),
      has_schema = logical(),
      roundtrip = logical(),
      supports_roi_mask = logical(),
      supports_time_idx = logical(),
      min_dims = integer(),
      stringsAsFactors = FALSE
    ))
  }
  
  # Build summary data frame
  do.call(rbind, lapply(names(transforms), function(type) {
    info <- transforms[[type]]
    data.frame(
      type = type,
      category = info$metadata$category %||% "other",
      lifecycle = info$lifecycle %||% "custom",
      has_forward = info$has_forward,
      has_invert = info$has_invert,
      has_schema = info$has_schema,
      roundtrip = isTRUE(info$roundtrip),
      supports_roi_mask = isTRUE(info$supports_roi_mask),
      supports_time_idx = isTRUE(info$supports_time_idx),
      min_dims = info$min_dims,
      compression = isTRUE(info$capabilities$compression),
      spatial = isTRUE(info$capabilities$spatial),
      temporal = !is.null(info$capabilities$temporal_basis),
      stringsAsFactors = FALSE
    )
  }))
}
