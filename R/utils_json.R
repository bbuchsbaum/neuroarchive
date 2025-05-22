#' JSON Read/Write Helpers for LNA Descriptors
#'
#' @description Provides internal functions for reading and writing transform
#'   descriptors stored as JSON strings within HDF5 datasets.
#'
#' @import jsonlite
#' @import hdf5r
#' @importFrom hdf5r H5File H5Group H5S H5T H5T_CLASS H5D h5types H5T_C_S1 H5T_CSET_UTF8
#' @keywords internal

#' Read a JSON descriptor from an HDF5 group.
#'
#' @param h5_group An H5Group object from hdf5r.
#' @param name The name of the HDF5 dataset containing the JSON string.
#' @return A list object parsed from the JSON string.
#' @details Assumes the dataset stores a single UTF-8 string (potentially variable length). Numeric
#'   values in the JSON are coerced with `as.numeric()` so that whole-number values are not
#'   returned as integers.
read_json_descriptor <- function(h5_group, name) {
  stopifnot(inherits(h5_group, "H5Group")) # Basic type check
  stopifnot(is.character(name), length(name) == 1)


  assert_h5_path(h5_group, name)

  json_string <- NULL
  parsed_list <- NULL

  loc <- sprintf("read_json_descriptor:%s", name)

  tryCatch({
    json_string <- h5_read(h5_group, name)

    if (length(json_string) != 1 || !is.character(json_string)) {
      abort_lna(
        sprintf("Dataset '%s' did not contain a single string.", name),
        .subclass = "lna_error_invalid_descriptor",
        location = loc
      )
    }

    parsed_list <- jsonlite::fromJSON(
      json_string,
      simplifyVector = TRUE,
      simplifyDataFrame = FALSE,
      simplifyMatrix = FALSE
    )

    # Convert any integer values to base numeric to avoid integer
    # coercion when numbers appear as whole values in the JSON
    convert_numeric <- function(x) {
      if (is.list(x)) {
        lapply(x, convert_numeric)
      } else if (is.integer(x)) {
        as.numeric(x)
      } else if (is.numeric(x) && all(!is.na(x) & x == floor(x))) {
        as.numeric(x)
      } else {
        x
      }
    }

    parsed_list <- convert_numeric(parsed_list)
  }, error = function(e) {
    detailed_error <- tryCatch(
      conditionMessage(e),
      error = function(e2) paste("Failed to get message:", e2$message)
    )
    abort_lna(
      sprintf(
        "Error reading/parsing JSON descriptor '%s': %s",
        name,
        detailed_error
      ),
      .subclass = "lna_error_json_parse",
      location = loc,
      parent = e
    )
  })

  return(parsed_list)
}

#' Write a JSON descriptor to an HDF5 group.
#'
#' @param h5_group An H5Group object from hdf5r.
#' @param name The name of the HDF5 dataset to create or overwrite.
#' @param desc_list A list object to be converted to JSON.
#' @return Invisibly returns NULL.
#' @details Writes the list as a JSON string to a scalar HDF5 dataset with
#'   a variable-length string datatype (UTF-8). Overwrites existing dataset
#'   with the same name.
write_json_descriptor <- function(h5_group, name, desc_list) {
  stopifnot(inherits(h5_group, "H5Group"))
  stopifnot(is.character(name) && length(name) == 1)
  stopifnot(is.list(desc_list))

  json_string <- jsonlite::toJSON(desc_list, auto_unbox = TRUE, pretty = TRUE)

  if (h5_group$exists(name)) {
    h5_group$link_delete(name)     # overwrite semantics
  }

  # Define resources, ensure cleanup with on.exit
  str_type <- NULL
  space <- NULL
  dset <- NULL
  on.exit({
    # Close resources if they were successfully created
    if (!is.null(str_type) && inherits(str_type, "H5T")) str_type$close()
    if (!is.null(space) && inherits(space, "H5S")) space$close()
    if (!is.null(dset) && inherits(dset, "H5D")) dset$close()
  }, add = TRUE)

  # --- Define scalar, variable-length, UTF-8 string dataset ----
  # Use C-style string datatype and set to variable length
  str_type <- hdf5r::H5T_C_S1$new()
  str_type$set_size(Inf)
  str_type$set_cset("UTF-8")
  space <- hdf5r::H5S$new("scalar")                   # scalar (0-dim) dataspace

  # Create the dataset skeleton
  dset <- h5_group$create_dataset(name,
                                  dtype = str_type,
                                  space = space,
                                  chunk_dims = NULL)

  # Write data using slice assignment
  dset[] <- json_string

  invisible(NULL)
} 
#' Schema Cache Environment
#'
#' Internal environment used to store compiled JSON schema objects for
#' transform validation.  It is not intended for direct use but can be
#' emptied via [schema_cache_clear()] when needed (e.g. during unit
#' testing).
#' @keywords internal
.schema_cache <- new.env(parent = emptyenv())

#' Clear the schema cache
#'
#' Removes all entries from the internal \code{.schema_cache} environment.
#' Intended primarily for unit tests or to avoid stale compiled objects.
#'
#' @return invisible(NULL)
#' @keywords internal
schema_cache_clear <- function() {
  rm(list = ls(envir = .schema_cache, all.names = TRUE), envir = .schema_cache)
  invisible(NULL)
}

