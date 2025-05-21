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
#' @details Assumes the dataset stores a single UTF-8 string (potentially variable length).
read_json_descriptor <- function(h5_group, name) {
  stopifnot(inherits(h5_group, "H5Group")) # Basic type check
  stopifnot(is.character(name), length(name) == 1)

  
  if (!h5_group$exists(name)) {
    # Or should this error? Returning NULL allows checking existence implicitly.
    # Let's error for now, as descriptors are usually expected.
    # TODO: Revisit error handling based on usage context.
    stop(paste("JSON descriptor dataset '", name, "' not found in HDF5 group.", sep = ""))
    # return(NULL)
  }


  dset <- NULL
  json_string <- NULL
  parsed_list <- NULL

  tryCatch({
    dset <- h5_group[[name]]
    # Simple read call
    json_string <- dset$read()

    if (length(json_string) != 1 || !is.character(json_string)) {
        stop(paste("Dataset '", name, "' did not contain a single string.", sep=""))
    }

    # Parse the JSON, simplifying vectors
    parsed_list <- jsonlite::fromJSON(json_string, simplifyVector = TRUE, simplifyDataFrame = FALSE, simplifyMatrix = FALSE)

  # Restore robust error handler
  }, error = function(e) {
      # Use conditionMessage for robustness
      detailed_error <- tryCatch(conditionMessage(e), error = function(e2) paste("Failed to get message:", e2$message))
      stop(paste("Error reading/parsing JSON descriptor '", name, "': ", detailed_error, sep = ""))
  }, finally = {
      if (!is.null(dset) && inherits(dset, "H5D")) dset$close()
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
  str_type <- H5T_STRING$new(size = Inf)
  str_type$set_cset("UTF-8")
  space <- H5S$new("scalar")                   # scalar (0-dim) dataspace

  # Create the dataset skeleton
  dset <- h5_group$create_dataset(name,
                                  dtype = str_type,
                                  space = space,
                                  chunk_dim = NULL)

  # Write data using slice assignment
  dset[] <- json_string

  invisible(NULL)
} 