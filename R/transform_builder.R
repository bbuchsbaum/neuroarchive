#' Transform Builder
#'
#' @description
#' R6 class that provides a builder pattern for transform implementations,
#' reducing boilerplate code and ensuring consistency.
#'
#' @import R6
#' @keywords internal

#' @export
#' @importFrom R6 R6Class
TransformBuilder <- R6::R6Class("TransformBuilder",
  public = list(
    #' @description
    #' Create a new TransformBuilder
    #' @param handle DataHandle object
    #' @param type Transform type string
    #' @param desc Transform descriptor
    initialize = function(handle, type, desc) {
      private$handle <- handle
      private$type <- type
      private$desc <- desc
      private$plan <- handle$plan
      private$run_id <- sanitize_run_id(handle$current_run_id %||% "run-01")
      private$fname <- private$plan$get_next_filename(type)
      private$base_name <- tools::file_path_sans_ext(private$fname)
      private$datasets <- list()
      private$step_index <- private$plan$next_index
    },
    
    #' @description
    #' Add a dataset to the transform
    #' @param path Dataset path
    #' @param data Data to store (or NULL if written directly to HDF5)
    #' @param role Dataset role
    #' @param dtype Data type string
    #' @return Self for chaining
    add_dataset = function(path, data, role, dtype = NA_character_) {
      private$datasets[[length(private$datasets) + 1]] <- list(
        path = path,
        data = data,
        role = role,
        dtype = dtype,
        payload_key = if (is.null(data)) "" else path
      )
      invisible(self)
    },
    
    #' @description
    #' Add a standard dataset with automatic path generation
    #' @param dataset_type Type of dataset (e.g., "quantized", "coefficients")
    #' @param data Data to store
    #' @param role Dataset role (defaults to dataset_type)
    #' @param dtype Data type string
    #' @return Self for chaining
    add_standard_dataset = function(dataset_type, data, role = NULL, dtype = NA_character_) {
      if (is.null(role)) {
        role <- dataset_type
      }
      
      path <- generate_transform_path(
        private$type, 
        private$run_id, 
        private$base_name,
        dataset_type
      )
      
      self$add_dataset(path, data, role, dtype)
    },
    
    #' @description
    #' Add a report dataset
    #' @param report_data Report data (will be JSON encoded and gzipped)
    #' @param report_type Type of report (default: "report")
    #' @return Self for chaining
    add_report = function(report_data, report_type = "report") {
      report_path <- generate_transform_path(
        private$type,
        private$run_id,
        private$base_name,
        report_type
      )
      
      # Convert to JSON and gzip
      json_str <- jsonlite::toJSON(report_data, auto_unbox = TRUE, pretty = TRUE)
      gzipped <- memCompress(charToRaw(as.character(json_str)), type = "gzip")
      
      # Store in params for reference
      if (is.null(private$desc$params)) {
        private$desc$params <- list()
      }
      private$desc$params$report_path <- report_path
      
      # Add as dataset
      self$add_dataset(report_path, gzipped, paste0(private$type, "_report"), "uint8")
    },
    
    #' @description
    #' Set transform version
    #' @param version Version string
    #' @return Self for chaining
    set_version = function(version = "1.0") {
      private$desc$version <- version
      invisible(self)
    },
    
    #' @description
    #' Set input/output specifications
    #' @param inputs Character vector of input keys
    #' @param outputs Character vector of output keys
    #' @return Self for chaining
    set_io = function(inputs = NULL, outputs = NULL) {
      if (!is.null(inputs)) {
        private$desc$inputs <- inputs
      }
      if (!is.null(outputs)) {
        private$desc$outputs <- outputs
      }
      invisible(self)
    },
    
    #' @description
    #' Build the transform and update handle
    #' @param stash_keys Keys to remove from stash
    #' @param stash_values Named list of values to add to stash
    #' @return Updated handle
    build = function(stash_keys = character(), stash_values = list()) {
      # Ensure descriptor has required fields
      private$desc$type <- private$type
      if (is.null(private$desc$version)) {
        private$desc$version <- "1.0"
      }
      
      # Add descriptor
      private$plan$add_descriptor(private$fname, private$desc)
      
      # Generate params JSON once
      params_json <- as.character(jsonlite::toJSON(
        private$desc$params %||% list(), 
        auto_unbox = TRUE
      ))
      
      # Add all datasets
      for (dataset in private$datasets) {
        # Add payload if data provided
        if (!is.null(dataset$data)) {
          private$plan$add_payload(dataset$path, dataset$data)
        }
        
        # Add dataset definition
        private$plan$add_dataset_def(
          dataset$path,
          dataset$role,
          private$type,
          private$run_id,
          as.integer(private$step_index),
          params_json,
          dataset$payload_key,
          "eager",
          dtype = dataset$dtype
        )
      }
      
      # Update handle
      private$handle$plan <- private$plan
      
      # Update stash
      if (length(stash_keys) > 0 || length(stash_values) > 0) {
        private$handle <- private$handle$update_stash(stash_keys, stash_values)
      }
      
      private$handle
    },
    
    #' @description
    #' Get the base name for this transform
    #' @return Character string
    get_base_name = function() {
      private$base_name
    },
    
    #' @description
    #' Get the run ID
    #' @return Character string
    get_run_id = function() {
      private$run_id
    },
    
    #' @description
    #' Get or set transform parameters
    #' @param params Named list of parameters (if setting)
    #' @return Current parameters (if getting)
    params = function(params = NULL) {
      if (!is.null(params)) {
        private$desc$params <- params
        invisible(self)
      } else {
        private$desc$params
      }
    }
  ),
  
  private = list(
    handle = NULL,
    type = NULL,
    desc = NULL,
    plan = NULL,
    run_id = NULL,
    fname = NULL,
    base_name = NULL,
    datasets = NULL,
    step_index = NULL
  )
)