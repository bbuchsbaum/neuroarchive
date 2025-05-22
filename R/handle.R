#' DataHandle Class for LNA Operations
#'
#' @description Represents the state during LNA read/write operations, holding
#'   data being transformed, metadata, the write plan (if applicable),
#'   HDF5 file access, and subsetting information.
#' @importFrom R6 R6Class
#' @import rlang
#' @keywords internal
DataHandle <- R6::R6Class("DataHandle",
  public = list(
    #' @field stash A list holding temporary data objects during transform chain.
    stash = NULL,
    #' @field meta A list holding metadata associated with the data.
    meta = NULL,
    #' @field plan A Plan object (from R/plan.R) used during write operations.
    plan = NULL,
    #' @field h5 An H5File object (from hdf5r) providing access to the LNA file.
    h5 = NULL,
    #' @field subset A list specifying subsetting parameters (e.g., ROI, time indices).
    subset = NULL,
    #' @field run_ids Character vector of run identifiers for multi-run data.
    run_ids = NULL,
    #' @field current_run_id The run identifier currently being processed.
    current_run_id = NULL,
    #' @field mask_info List with mask array and active voxel count
    mask_info = NULL,

    #' @description
    #' Initialize a new DataHandle object.
    #' @param initial_stash Initial list of objects for the stash.
    #' @param initial_meta Initial list for metadata.
    #' @param plan A Plan object (optional, for writing).
    #' @param h5 An H5File object (optional, for reading/writing).
    #' @param subset A list specifying subsetting (optional, for reading).
    initialize = function(initial_stash = list(), initial_meta = list(), plan = NULL,
                          h5 = NULL, subset = list(), run_ids = character(),
                          current_run_id = NULL, mask_info = NULL) {
      # Basic input validation
      stopifnot(is.list(initial_stash))
      stopifnot(is.list(initial_meta))
      stopifnot(is.list(subset))
      # Placeholder validation for R6 objects - refine later if needed
      if (!is.null(plan) && !inherits(plan, "Plan")) {
        stop("'plan' must be a Plan R6 object or NULL")
      }
      if (!is.null(h5) && !inherits(h5, "H5File")) {
        # Assuming hdf5r class is H5File - verify this
        stop("'h5' must be an H5File object from hdf5r or NULL")
      }

      stopifnot(is.character(run_ids))
      if (!is.null(current_run_id)) {
        stopifnot(is.character(current_run_id), length(current_run_id) == 1)
      }

      self$stash <- initial_stash
      self$meta <- initial_meta
      self$plan <- plan
      self$h5 <- h5
      self$subset <- subset
      self$run_ids <- run_ids
      self$current_run_id <- current_run_id
      self$mask_info <- mask_info
    },

    #' @description
    #' Retrieve specified input objects from the stash.
    #' @param keys Character vector of keys to retrieve from the stash.
    #' @return A named list containing the requested objects.
    #' @details Raises an lna_error_contract if any key is not found.
    get_inputs = function(keys) {
      stopifnot(is.character(keys), length(keys) > 0)
      stash_names <- names(self$stash)
      missing_keys <- setdiff(keys, stash_names)

      if (length(missing_keys) > 0) {
        abort_lna(
          paste(
            "Required key(s) not found in stash:",
            paste(missing_keys, collapse = ", ")
          ),
          .subclass = "lna_error_contract",
          missing_keys = missing_keys,
          location = "DataHandle$get_inputs"
        )
      }
      return(self$stash[keys])
    },

    #' @description
    #' Update the stash with new or modified objects (immutable update).
    #' @param keys Character vector of keys to remove from the current stash.
    #' @param new_values Named list of new objects to add to the stash.
    #' @return A *new* DataHandle object with the updated stash.
    update_stash = function(keys, new_values) {
      message(sprintf("[DataHandle$update_stash ENTRY] Self stash keys: %s. ", paste(names(self$stash), collapse=", ")))
      message(sprintf("[DataHandle$update_stash ENTRY] new_values keys: %s. ", paste(names(new_values), collapse=", ")))
  
      stopifnot(is.character(keys))
      stopifnot(is.list(new_values))

      # Calculate the new stash based on current stash, keys to remove, and new values
      current_stash <- self$stash
      keys_to_remove <- intersect(keys, names(current_stash))
      if (length(keys_to_remove) > 0) {
         current_stash[keys_to_remove] <- NULL
      }

      # Warn if new_values will overwrite existing stash entries that were not removed
      if (length(new_values) > 0) {
          overlap <- intersect(names(new_values), names(current_stash))
          if (length(overlap) > 0) {
              warning(
                paste(
                  "Overwriting existing stash entries:",
                  paste(overlap, collapse = ", ")
                )
              )
          }
          # Use modifyList for safe merging/overwriting
          current_stash <- utils::modifyList(current_stash, new_values)
      }
      
      # Return a new DataHandle with the updated stash using the 'with' method
      return(self$with(stash = current_stash))
    },

    #' @description
    #' Create a new DataHandle with modified fields (immutable update).
    #' @param ... Named arguments corresponding to fields to update (e.g., meta = new_meta).
    #' @return A *new* DataHandle object with updated fields.
    with = function(...) {
      new_obj <- self$clone() # Use shallow clone; deep clone was causing issues with H5File objects
      updates <- list(...)
      # Get public field names from the R6 class generator
      class_generator <- get(class(self)[1])
      allowed_fields <- names(class_generator$public_fields)

      for (field_name in names(updates)) {
        if (!field_name %in% allowed_fields) {
          warning(paste("Field '", field_name, "' not found in DataHandle, skipping update.", sep = ""))
          next
        }
        # TODO: Add validation specific to field types? (e.g., plan must be Plan)
        new_obj[[field_name]] <- updates[[field_name]]
      }
      message(sprintf("[DataHandle$with] Returning new_obj. Stash keys: %s. Is input in stash NULL? %s", paste(names(new_obj$stash), collapse=", "), is.null(new_obj$stash$input)))
      return(new_obj)
    },

    #' @description
    #' Check if a key exists in the stash.
    #' @param key Character string, the key to check.
    #' @return Logical, TRUE if the key exists in the stash, FALSE otherwise.
    has_key = function(key) {
      stopifnot(is.character(key), length(key) == 1)
      return(key %in% names(self$stash))
    }
  )
)
