#' Plan Class for LNA Write Operations
#'
#' @description Defines the structure and operations for planning the write
#'   process of an LNA file, including dataset definitions, transform descriptors,
#'   and payload management.
#' @importFrom R6 R6Class
#' @import tibble
#' @import jsonlite
#' @keywords internal
Plan <- R6::R6Class("Plan",
  public = list(
    #' @field datasets A tibble storing definitions for HDF5 datasets to be created.
    datasets = NULL,
    #' @field descriptors A list storing transform descriptor lists.
    descriptors = NULL,
    #' @field payloads A list storing data payloads to be written.
    payloads = NULL,
    #' @field next_index An integer counter for naming transforms sequentially.
    next_index = NULL,
    #' @field origin_label A string label identifying the source (e.g., run ID).
    origin_label = NULL, # Added for clarity based on spec usage

    #' @description
    #' Initialize a new Plan object.
    #' @param origin_label A string label for the origin (e.g., run ID).
    initialize = function(origin_label = "global") {
      if (!is.character(origin_label) || length(origin_label) != 1) {
        abort_lna(
          "origin_label must be a single character string",
          .subclass = "lna_error_validation",
          location = "Plan$initialize"
        )
      }
      self$datasets <- tibble::tibble(
        path = character(),
        role = character(),
        producer = character(),
        origin = character(),
        step_index = integer(),
        params_json = character(),
        payload_key = character(),
        write_mode = character(),
        write_mode_effective = character(), # Added based on Spec v1.4
        dtype = character()
      )
      self$descriptors <- list()
      self$payloads <- list()
      self$next_index <- 0L
      self$origin_label <- origin_label
    },

    #' @description
    #' Add a data payload to be written later.
    #' @param key Character string identifier (often HDF5 path).
    #' @param value The R object to be written.
    #' @param overwrite Logical flag; if `TRUE`, an existing payload with the
    #'   same key will be replaced. Defaults to `FALSE` which raises an error on
    #'   duplicates.
    add_payload = function(key, value, overwrite = FALSE) {
      if (!is.character(key) || length(key) != 1) {
        abort_lna(
          "key must be a single character string",
          .subclass = "lna_error_validation",
          location = "Plan$add_payload"
        )
      }
      if (!is.logical(overwrite) || length(overwrite) != 1) {
        abort_lna(
          "overwrite must be a single logical value",
          .subclass = "lna_error_validation",
          location = "Plan$add_payload"
        )
      }
      if (key %in% names(self$payloads) && !overwrite) {
        abort_lna(
          sprintf("Payload key '%s' already exists in plan", key),
          .subclass = "lna_error_duplicate_key",
          location = "Plan$add_payload"
        )
      }
      self$payloads[[key]] <- value
      invisible(self)
    },

    #' @description
    #' Add a definition for an HDF5 dataset.
    #' @param path Character string, HDF5 path for the dataset.
    #' @param role Character string, semantic role of the dataset.
    #' @param producer Character string, type of the transform producing this.
    #' @param origin Character string, label of the originating run/source.
    #' @param step_index Integer, index of the transform step.
    #' @param params_json Character string, JSON representation of transform params.
    #' @param payload_key Character string, key linking to the entry in `self$payloads`.
    #' @param write_mode Character string, requested write mode ("eager"/"stream").
    #' @param dtype Optional character string naming the storage datatype (e.g.,
    #'   "uint8", "uint16").
    add_dataset_def = function(path, role, producer, origin, step_index, params_json, payload_key, write_mode, dtype = NA_character_) {
      # Basic type checks with additional validation
      if (!is.character(path) || length(path) != 1) {
        abort_lna("path must be a single character string", .subclass = "lna_error_validation", location = "Plan$add_dataset_def")
      }
      if (!is.character(role) || length(role) != 1) {
        abort_lna("role must be a single character string", .subclass = "lna_error_validation", location = "Plan$add_dataset_def")
      }
      if (!is.character(producer) || length(producer) != 1) {
        abort_lna("producer must be a single character string", .subclass = "lna_error_validation", location = "Plan$add_dataset_def")
      }
      if (!is.character(origin) || length(origin) != 1) {
        abort_lna("origin must be a single character string", .subclass = "lna_error_validation", location = "Plan$add_dataset_def")
      }
      if (!is.numeric(step_index) || length(step_index) != 1 || is.na(step_index) || step_index %% 1 != 0) {
        abort_lna("step_index must be a single integer value", .subclass = "lna_error_validation", location = "Plan$add_dataset_def")
      }
      if (!is.character(params_json) || length(params_json) != 1) {
        abort_lna("params_json must be a single character string", .subclass = "lna_error_validation", location = "Plan$add_dataset_def")
      }
      if (!is.character(payload_key) || length(payload_key) != 1) {
        abort_lna("payload_key must be a single character string", .subclass = "lna_error_validation", location = "Plan$add_dataset_def")
      }
      if (!is.character(write_mode) || length(write_mode) != 1) {
        abort_lna("write_mode must be a single character string", .subclass = "lna_error_validation", location = "Plan$add_dataset_def")
      }
      if (!is.character(dtype) || length(dtype) != 1) {
        abort_lna("dtype must be a single character string", .subclass = "lna_error_validation", location = "Plan$add_dataset_def")
      }

      # Validate write_mode values
      if (!write_mode %in% c("eager", "stream")) {
        abort_lna(
          "write_mode must be either 'eager' or 'stream'",
          .subclass = "lna_error_validation",
          location = "Plan$add_dataset_def"
        )
      }

      # Validate JSON
      valid_json <- jsonlite::validate(params_json)
      if (!isTRUE(valid_json)) {
        abort_lna(
          sprintf("Invalid params_json: %s", as.character(valid_json)),
          .subclass = "lna_error_validation",
          location = "Plan$add_dataset_def"
        )
      }

      self$datasets <- tibble::add_row(
        self$datasets,
        path = path,
        role = role,
        producer = producer,
        origin = origin,
        step_index = as.integer(step_index),
        params_json = params_json,
        payload_key = payload_key,
        write_mode = write_mode,
        write_mode_effective = NA_character_, # To be filled during materialization
        dtype = as.character(dtype)
      )
      invisible(self)
    },

    #' @description
    #' Add a transform descriptor to the plan.
    #' @param transform_name Character string, name for the descriptor (e.g., "00_type.json").
    #' @param desc_list List, the descriptor content.
    add_descriptor = function(transform_name, desc_list) {
      if (!is.character(transform_name) || length(transform_name) != 1) {
        abort_lna(
          "transform_name must be a single character string",
          .subclass = "lna_error_validation",
          location = "Plan$add_descriptor"
        )
      }
      if (!is.list(desc_list)) {
        abort_lna(
          "desc_list must be a list",
          .subclass = "lna_error_validation",
          location = "Plan$add_descriptor"
        )
      }
      if (transform_name %in% names(self$descriptors)) {
        abort_lna(
          sprintf("Descriptor name '%s' already exists in plan", transform_name),
          .subclass = "lna_error_duplicate_key",
          location = "Plan$add_descriptor"
        )
      }
      self$descriptors[[transform_name]] <- desc_list
      self$next_index <- self$next_index + 1L
      invisible(self)
    },

    #' @description
    #' Get the next sequential filename prefix for a transform descriptor.
    #' @param type Character string, the transform type.
    #' @return Character string (e.g., "00_type.json").
    get_next_filename = function(type) {
      if (!is.character(type) || length(type) != 1) {
        abort_lna(
          "type must be a single character string",
          .subclass = "lna_error_validation",
          location = "Plan$get_next_filename"
        )
      }

      if (grepl("..", type, fixed = TRUE) || grepl("/", type, fixed = TRUE) || grepl("\\", type, fixed = TRUE)) {
        abort_lna(
          sprintf("Invalid characters found in type '%s'", type),
          .subclass = "lna_error_validation",
          location = "Plan$get_next_filename"
        )
      }

      safe_pat <- "^[A-Za-z][A-Za-z0-9_.]*$"
      if (!grepl(safe_pat, type)) {
        abort_lna(
          sprintf("Invalid transform type '%s'. Must match %s", type, safe_pat),
          .subclass = "lna_error_validation",
          location = "Plan$get_next_filename"
        )
      }

      index_str <- sprintf("%02d", self$next_index)
      filename <- paste0(index_str, "_", type, ".json")
      return(filename)
    },

    #' @description
    #' Return the first run identifier appearing in the plan. If no dataset
    #' definitions exist, fall back to `origin_label` when it matches the run
    #' pattern. Returns `NULL` when no run information is available.
    first_run_id = function() {
      if (nrow(self$datasets) > 0) {
        ids <- self$datasets$origin
        run_like <- grep("^run-[0-9]+$", ids, value = TRUE)
        if (length(run_like) > 0) return(run_like[1])
      }
      if (grepl("^run-[0-9]+$", self$origin_label)) {
        return(self$origin_label)
      }
      NULL
    },

    #' @description
    #' Convenience helper to add an array as the initial payload for a run.
    #' This is used by `core_write` when no transforms are specified.
    #' @param x Array to add.
    #' @param run_id Optional run identifier. Defaults to "run-01" when neither
    #'   `run_id` nor `origin_label` specifies a run pattern.
    import_from_array = function(x, run_id = NULL) {
      if (!is.array(x)) {
        abort_lna(
          "x must be an array",
          .subclass = "lna_error_validation",
          location = "Plan$import_from_array"
        )
      }
      rid <- if (!is.null(run_id)) {
        run_id
      } else if (grepl("^run-[0-9]+$", self$origin_label)) {
        self$origin_label
      } else {
        "run-01"
      }
      key <- paste0(rid, "_initial")
      self$add_payload(key, x, overwrite = TRUE)
      self$add_dataset_def(
        path = file.path("/scans", rid, "data", "values"),
        role = "raw_data",
        producer = "core_write_initial_input",
        origin = rid,
        step_index = 0L,
        params_json = "{}",
        payload_key = key,
        write_mode = "eager",
        dtype = NA_character_
      )
      invisible(self)
    },

    #' @description
    #' Mark a payload as written (e.g., by setting its value to NULL).
    #' @param key Character string, the key of the payload to mark.
    mark_payload_written = function(key) {
      if (!is.character(key) || length(key) != 1) {
        abort_lna(
          "key must be a single character string",
          .subclass = "lna_error_validation",
          location = "Plan$mark_payload_written"
        )
      }
      if (!key %in% names(self$payloads)) {
        warn_lna(
          sprintf("Payload key '%s' not found in plan when trying to mark as written", key),
          .subclass = "lna_warning_missing_key",
          location = "Plan$mark_payload_written"
        )
      } else {
        self$payloads[[key]] <- NULL
      }
      invisible(self)
    }
  )
) 