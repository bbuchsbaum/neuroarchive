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
    #' @field stash A list environment holding temporary data objects during transform chain.
    stash = NULL,
    #' @field meta A list holding metadata associated with the data.
    meta = NULL,
    #' @field plan A Plan object (from R/plan.R) used during write operations.
    plan = NULL,
    #' @field h5 An H5File object (from hdf5r) providing access to the LNA file.
    h5 = NULL,
    #' @field subset A list specifying subsetting parameters (e.g., ROI, time indices).
    subset = NULL,

    #' @description
    #' Initialize a new DataHandle object.
    #' @param initial_stash Initial list of objects for the stash.
    #' @param initial_meta Initial list for metadata.
    #' @param plan A Plan object (optional, for writing).
    #' @param h5 An H5File object (optional, for reading/writing).
    #' @param subset A list specifying subsetting (optional, for reading).
    initialize = function(initial_stash = list(), initial_meta = list(), plan = NULL, h5 = NULL, subset = list()) {
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

      self$stash <- initial_stash
      self$meta <- initial_meta
      self$plan <- plan
      self$h5 <- h5
      self$subset <- subset
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
          missing_keys = missing_keys
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
      stopifnot(is.character(keys))
      stopifnot(is.list(new_values))

      # Calculate the new stash based on current stash, keys to remove, and new values
      current_stash <- self$stash
      keys_to_remove <- intersect(keys, names(current_stash))
      if (length(keys_to_remove) > 0) {
         current_stash[keys_to_remove] <- NULL
      }

      # Add new values (potentially overwriting existing keys if not removed)
      # Consider if overwriting should be prevented/warned if key wasn't in `keys`?
      # For now, allow overwriting as list assignment does.
      if (length(new_values) > 0) {
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
      new_obj <- self$clone(deep = TRUE) # Use deep clone for safety with lists
      updates <- list(...)
      allowed_fields <- names(get(class(self)[1])$public_fields)

      for (field_name in names(updates)) {
        if (!field_name %in% allowed_fields) {
          warning(paste("Field '", field_name, "' not found in DataHandle, skipping update.", sep = ""))
          next
        }
        # TODO: Add validation specific to field types? (e.g., plan must be Plan)
        new_obj[[field_name]] <- updates[[field_name]]
      }
      return(new_obj)
    },

    #' @description
    #' Check if a key exists in the stash.
    #' @param key Character string, the key to check.
    #' @return Logical, TRUE if the key exists in the stash, FALSE otherwise.
    exists = function(key) {
      stopifnot(is.character(key), length(key) == 1)
      return(key %in% names(self$stash))
    }
  )
) 