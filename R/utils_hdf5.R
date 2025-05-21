#' HDF5 Attribute Read/Write Helpers
#'
#' @description Provides internal functions for reading, writing, checking existence,
#'   and deleting attributes associated with HDF5 objects (groups or datasets).
#'
#' @import hdf5r
#' @keywords internal

# Check if the object is a valid hdf5r object that can hold attributes
.is_valid_h5_object <- function(h5_obj) {
  inherits(h5_obj, "H5Group") || inherits(h5_obj, "H5D")
}

#' Write an attribute to an HDF5 object.
#'
#' @param h5_obj An H5Group or H5D object from hdf5r.
#' @param name The name of the attribute.
#' @param value The value to write (basic R types and simple vectors supported).
#' @return Invisibly returns NULL.
#' @details Overwrites the attribute if it already exists.
h5_attr_write <- function(h5_obj, name, value) {
  stopifnot("h5_obj must be an H5Group or H5D object" = .is_valid_h5_object(h5_obj))
  stopifnot(is.character(name), length(name) == 1)

  # Use hdf5r's assignment function, which handles overwriting
  tryCatch({
    hdf5r::h5attr(h5_obj, name) <- value
  }, error = function(e) {
    stop(paste("Error writing attribute '", name, "': ", conditionMessage(e)), call. = FALSE)
  })

  invisible(NULL)
}

#' Read an attribute from an HDF5 object.
#'
#' @param h5_obj An H5Group or H5D object from hdf5r.
#' @param name The name of the attribute.
#' @return The value of the attribute.
#' @details Throws an error if the attribute does not exist.
h5_attr_read <- function(h5_obj, name) {
  stopifnot("h5_obj must be an H5Group or H5D object" = .is_valid_h5_object(h5_obj))
  stopifnot(is.character(name), length(name) == 1)

  if (!h5_attr_exists(h5_obj, name)) {
    stop(paste("Attribute '", name, "' not found."), call. = FALSE)
  }

  # Use hdf5r's read function
  tryCatch({
    hdf5r::h5attr(h5_obj, name)
  }, error = function(e) {
    stop(paste("Error reading attribute '", name, "': ", conditionMessage(e)), call. = FALSE)
  })
}

#' Check if an attribute exists on an HDF5 object.
#'
#' @param h5_obj An H5Group or H5D object from hdf5r.
#' @param name The name of the attribute.
#' @return Logical TRUE if the attribute exists, FALSE otherwise.
h5_attr_exists <- function(h5_obj, name) {
  stopifnot("h5_obj must be an H5Group or H5D object" = .is_valid_h5_object(h5_obj))
  stopifnot(is.character(name), length(name) == 1)

  tryCatch({
    h5_obj$attr_exists(name)
  }, error = function(e) {
    # Should generally not error, but catch just in case
    stop(paste("Error checking existence of attribute '", name, "': ", conditionMessage(e)), call. = FALSE)
  })
}

#' Delete an attribute from an HDF5 object.
#'
#' @param h5_obj An H5Group or H5D object from hdf5r.
#' @param name The name of the attribute to delete.
#' @return Invisibly returns NULL.
#' @details Does nothing if the attribute does not exist.
h5_attr_delete <- function(h5_obj, name) {
  stopifnot("h5_obj must be an H5Group or H5D object" = .is_valid_h5_object(h5_obj))
  stopifnot(is.character(name), length(name) == 1)

  # Check existence first to avoid potential error in attr_delete if it doesn't exist
  if (h5_attr_exists(h5_obj, name)) {
      tryCatch({
        h5_obj$attr_delete(name)
      }, error = function(e) {
        stop(paste("Error deleting attribute '", name, "': ", conditionMessage(e)), call. = FALSE)
      })
  }

  invisible(NULL)
} 