#' Plan Class for LNA Write Operations
#'
#' @description Defines the structure and operations for planning the write
#'   process of an LNA file, including dataset definitions, transform descriptors,
#'   and payload management.
#' @importFrom R6 R6Class
#' @import tibble
#' @keywords internal
Plan <- R6::R6Class("Plan",
  public = list(
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
      private$.datasets_cols <- list(
        path = list(),
        role = list(),
        producer = list(),
        origin = list(),
        step_index = list(),
        params_json = list(),
        payload_key = list(),
        write_mode = list(),
        write_mode_effective = list(),
        dtype = list()
      )
      private$.datasets_cache <- private$.materialize_datasets()
      private$.datasets_cache_valid <- TRUE
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
          sprintf("Payload key '%s' already exists in plan.", key),
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

      private$.append_dataset_row(
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
          sprintf("Descriptor name '%s' already exists in plan.", transform_name),
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
      datasets <- self$datasets
      if (nrow(datasets) > 0) {
        ids <- datasets$origin
        run_like <- grep("^run-[0-9]+$", ids, value = TRUE)
        if (length(run_like) > 0) return(run_like[1])
      }
      if (grepl("^run-[0-9]+$", self$origin_label)) {
        return(self$origin_label)
      }
      NULL
    },

    #' @description
    #' Set the effective write mode for a dataset definition.
    #' @param i Integer row index in `datasets`.
    #' @param mode Character scalar, e.g. `"eager"`.
    set_write_mode_effective = function(i, mode) {
      if (!is.numeric(i) || length(i) != 1 || is.na(i) || i %% 1 != 0) {
        abort_lna(
          "i must be a single integer index",
          .subclass = "lna_error_validation",
          location = "Plan$set_write_mode_effective"
        )
      }
      if (!is.character(mode) || length(mode) != 1 || !nzchar(mode)) {
        abort_lna(
          "mode must be a non-empty character scalar",
          .subclass = "lna_error_validation",
          location = "Plan$set_write_mode_effective"
        )
      }

      idx <- as.integer(i)
      n <- length(private$.datasets_cols$path)
      if (idx < 1L || idx > n) {
        abort_lna(
          sprintf("Dataset index %d out of bounds (n=%d)", idx, n),
          .subclass = "lna_error_validation",
          location = "Plan$set_write_mode_effective"
        )
      }

      private$.datasets_cols$write_mode_effective[[idx]] <- mode
      if (isTRUE(private$.datasets_cache_valid)) {
        private$.datasets_cache$write_mode_effective[[idx]] <- mode
      }
      invisible(self)
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
          sprintf("Payload key '%s' not found in plan when trying to mark as written.", key),
          .subclass = "lna_warning_missing_key",
          location = "Plan$mark_payload_written"
        )
      } else {
        self$payloads[key] <- list(NULL)
      }
      invisible(self)
    }
  ),
  active = list(
    #' @field datasets A tibble storing definitions for HDF5 datasets to be created.
    datasets = function(value) {
      if (missing(value)) {
        if (!isTRUE(private$.datasets_cache_valid)) {
          private$.datasets_cache <- private$.materialize_datasets()
          private$.datasets_cache_valid <- TRUE
        }
        return(private$.datasets_cache)
      }

      if (!tibble::is_tibble(value)) {
        abort_lna(
          "datasets must be a tibble",
          .subclass = "lna_error_validation",
          location = "Plan$datasets<-"
        )
      }
      required_cols <- c(
        "path", "role", "producer", "origin", "step_index", "params_json",
        "payload_key", "write_mode", "write_mode_effective", "dtype"
      )
      missing_cols <- setdiff(required_cols, names(value))
      if (length(missing_cols) > 0) {
        abort_lna(
          paste0("datasets missing required column(s): ", paste(missing_cols, collapse = ", ")),
          .subclass = "lna_error_validation",
          location = "Plan$datasets<-"
        )
      }

      private$.datasets_cols <- lapply(required_cols, function(nm) {
        as.list(value[[nm]])
      })
      names(private$.datasets_cols) <- required_cols
      private$.datasets_cache <- value[, required_cols, drop = FALSE]
      private$.datasets_cache_valid <- TRUE
      invisible(value)
    }
  ),
  private = list(
    .datasets_cols = NULL,
    .datasets_cache = NULL,
    .datasets_cache_valid = FALSE,

    .append_dataset_row = function(path, role, producer, origin, step_index, params_json,
                                   payload_key, write_mode, write_mode_effective, dtype) {
      cols <- private$.datasets_cols
      cols$path[[length(cols$path) + 1L]] <- path
      cols$role[[length(cols$role) + 1L]] <- role
      cols$producer[[length(cols$producer) + 1L]] <- producer
      cols$origin[[length(cols$origin) + 1L]] <- origin
      cols$step_index[[length(cols$step_index) + 1L]] <- step_index
      cols$params_json[[length(cols$params_json) + 1L]] <- params_json
      cols$payload_key[[length(cols$payload_key) + 1L]] <- payload_key
      cols$write_mode[[length(cols$write_mode) + 1L]] <- write_mode
      cols$write_mode_effective[[length(cols$write_mode_effective) + 1L]] <- write_mode_effective
      cols$dtype[[length(cols$dtype) + 1L]] <- dtype
      private$.datasets_cols <- cols
      private$.datasets_cache_valid <- FALSE
      invisible(NULL)
    },

    .materialize_datasets = function() {
      cols <- private$.datasets_cols
      if (is.null(cols) || length(cols$path) == 0) {
        return(tibble::tibble(
          path = character(),
          role = character(),
          producer = character(),
          origin = character(),
          step_index = integer(),
          params_json = character(),
          payload_key = character(),
          write_mode = character(),
          write_mode_effective = character(),
          dtype = character()
        ))
      }

      tibble::tibble(
        path = as.character(unlist(cols$path, use.names = FALSE)),
        role = as.character(unlist(cols$role, use.names = FALSE)),
        producer = as.character(unlist(cols$producer, use.names = FALSE)),
        origin = as.character(unlist(cols$origin, use.names = FALSE)),
        step_index = as.integer(unlist(cols$step_index, use.names = FALSE)),
        params_json = as.character(unlist(cols$params_json, use.names = FALSE)),
        payload_key = as.character(unlist(cols$payload_key, use.names = FALSE)),
        write_mode = as.character(unlist(cols$write_mode, use.names = FALSE)),
        write_mode_effective = as.character(unlist(cols$write_mode_effective, use.names = FALSE)),
        dtype = as.character(unlist(cols$dtype, use.names = FALSE))
      )
    }
  )
)
