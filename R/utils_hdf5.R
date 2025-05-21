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

#' Write a dataset to an HDF5 group
#'
#' @description Creates or overwrites a dataset at `path`, optionally using
#'   chunking and gzip compression. Intermediate groups in `path` are created as
#'   needed. If `chunk_dims` is `NULL`, a basic heuristic limits chunk dimensions
#'   to at most 128 along each axis.
#'
#' @param h5_group An `H5Group` object used as the starting location for `path`.
#' @param path Character string giving the dataset path relative to `h5_group`.
#' @param data Numeric matrix/array to write.
#' @param chunk_dims Optional integer vector specifying HDF5 chunk dimensions.
#' @param compression_level Integer 0–9 giving gzip compression level.
#' @return The created `H5D` dataset object (invisibly).
h5_write_dataset <- function(h5_group, path, data,
                             chunk_dims = NULL, compression_level = 0) {
  stopifnot(inherits(h5_group, "H5Group"))
  stopifnot(is.character(path), length(path) == 1)
  stopifnot(is.numeric(compression_level), length(compression_level) == 1)

  if (!is.array(data)) {
    if (is.vector(data)) {
      dim(data) <- length(data)
    } else {
      stop("`data` must be a matrix or array")
    }
  }

  parts <- strsplit(path, "/")[[1]]
  parts <- parts[nzchar(parts)]
  stopifnot(length(parts) > 0)
  ds_name <- tail(parts, 1)

  grp <- h5_group
  if (length(parts) > 1) {
    for (g in parts[-length(parts)]) {
      grp <- if (!grp$exists(g)) grp$create_group(g) else grp[[g]]
    }
  }

  if (is.null(chunk_dims)) {
    chunk_dims <- pmin(dim(data), rep(128L, length(dim(data))))
  } else {
    chunk_dims <- as.integer(chunk_dims)
  }

  create_fun <- function(level) {
    grp$create_dataset(ds_name,
                       data = data,
                       chunk_dims = chunk_dims,
                       gzip_level = level)
  }

  dset <- NULL
  if (!is.null(compression_level) && compression_level > 0) {
    dset <- tryCatch(create_fun(compression_level), error = function(e) {
      warning("Compression filter unavailable, writing without compression")
      NULL
    })
    if (is.null(dset)) {
      dset <- create_fun(NULL)
    }
  } else {
    dset <- create_fun(NULL)
  }

  invisible(dset)
}
