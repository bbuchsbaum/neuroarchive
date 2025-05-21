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
      stopifnot(is.character(origin_label), length(origin_label) == 1)
      self$datasets <- tibble::tibble(
        path = character(),
        role = character(),
        producer = character(),
        origin = character(),
        step_index = integer(),
        params_json = character(),
        payload_key = character(),
        write_mode = character(),
        write_mode_effective = character() # Added based on Spec v1.4
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
      stopifnot(is.character(key), length(key) == 1)
      stopifnot(is.logical(overwrite), length(overwrite) == 1)
      if (key %in% names(self$payloads) && !overwrite) {
        stop(paste("Payload key '", key, "' already exists in plan.", sep = ""))
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
    add_dataset_def = function(path, role, producer, origin, step_index, params_json, payload_key, write_mode) {
      # Basic type checks with additional validation
      stopifnot(
        is.character(path), length(path) == 1,
        is.character(role), length(role) == 1,
        is.character(producer), length(producer) == 1,
        is.character(origin), length(origin) == 1,
        is.numeric(step_index), length(step_index) == 1, !is.na(step_index), step_index %% 1 == 0,
        is.character(params_json), length(params_json) == 1,
        is.character(payload_key), length(payload_key) == 1,
        is.character(write_mode), length(write_mode) == 1
      )

      # Validate write_mode values
      if (!write_mode %in% c("eager", "stream")) {
        stop("write_mode must be either 'eager' or 'stream'")
      }

      # Validate JSON
      valid_json <- jsonlite::validate(params_json)
      if (!isTRUE(valid_json)) {
        stop(paste("Invalid params_json:", valid_json))
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
        write_mode_effective = NA_character_ # To be filled during materialization
      )
      invisible(self)
    },

    #' @description
    #' Add a transform descriptor to the plan.
    #' @param transform_name Character string, name for the descriptor (e.g., "00_type.json").
    #' @param desc_list List, the descriptor content.
    add_descriptor = function(transform_name, desc_list) {
      stopifnot(
        is.character(transform_name), length(transform_name) == 1,
        is.list(desc_list)
      )
      if (transform_name %in% names(self$descriptors)) {
        stop(paste("Descriptor name '", transform_name, "' already exists in plan.", sep = ""))
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
      stopifnot(is.character(type), length(type) == 1)

      if (grepl("\.\.", type) || grepl("/", type) || grepl("\\\\", type)) {
        stop(sprintf(
          "Invalid characters found in type '%s'", type
        ), call. = FALSE)
      }

      safe_pat <- "^[A-Za-z][A-Za-z0-9_.]*$"
      if (!grepl(safe_pat, type)) {
        stop(sprintf(
          "Invalid transform type '%s'. Must match %s", type, safe_pat
        ), call. = FALSE)
      }

      index_str <- sprintf("%02d", self$next_index)
      filename <- paste0(index_str, "_", type, ".json")
      return(filename)
    },

    #' @description
    #' Mark a payload as written (e.g., by setting its value to NULL).
    #' @param key Character string, the key of the payload to mark.
    mark_payload_written = function(key) {
      stopifnot(is.character(key), length(key) == 1)
      if (!key %in% names(self$payloads)) {
        warning(paste("Payload key '", key, "' not found in plan when trying to mark as written.", sep = ""))
      } else {
        self$payloads[[key]] <- NULL
      }
      invisible(self)
    }
  )
) 