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

#' Guess reasonable HDF5 chunk dimensions
#'
#' @description Heuristic used when `chunk_dims` is `NULL` in
#'   `h5_write_dataset`. Chunks are targeted to ~1 MiB. For datasets
#'   larger than 4 GiB, the first dimension is halved until the estimated
#'   chunk size falls below 1 GiB. If the chunk would still exceed about
#'   256 MiB, an additional reduction is applied with a warning.
#' @param dims Integer vector of dataset dimensions.
#' @param dtype_size Size in bytes of a single data element.
#' @return Integer vector of chunk dimensions.
guess_chunk_dims <- function(dims, dtype_size) {
  target <- 1024^2  # 1 MiB
  chunk <- hdf5r::guess_chunks(space_maxdims = dims,
                               dtype_size = dtype_size,
                               chunk_size = target)

  data_bytes <- prod(dims) * dtype_size
  chunk_bytes <- prod(chunk) * dtype_size

  if (data_bytes > 4 * 1024^3 && chunk_bytes > 1024^3) {
    while (chunk_bytes > 1024^3 && chunk[1] > 1) {
      chunk[1] <- ceiling(chunk[1] / 2)
      chunk_bytes <- prod(chunk) * dtype_size
    }
  }

  if (chunk_bytes > 256 * 1024^2) {
    warning("Auto-reducing chunk size to meet HDF5 limits")
    while (chunk_bytes > 256 * 1024^2 && chunk[1] > 1) {
      chunk[1] <- ceiling(chunk[1] / 2)
      chunk_bytes <- prod(chunk) * dtype_size
    }
  }

  chunk <- pmin(as.integer(chunk), dims)
  return(chunk)
}

#' Write a dataset to an HDF5 group
#'
#' @description Creates or overwrites a dataset at `path`, optionally using
#'   chunking and gzip compression. Intermediate groups in `path` are created as
#'   needed. If `chunk_dims` is `NULL`, a heuristic attempts to keep chunks
#'   around 1 MiB. For datasets larger than 4 GiB, the fastest changing axis is
#'   halved until the estimated chunk size is below 1 GiB. If the resulting chunk
#'   would still exceed roughly 256 MiB (HDF5 practical limit), an additional
#'   reduction is performed with a warning.
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
    element_size <- if (is.integer(data)) 4L else 8L
    chunk_dims <- guess_chunk_dims(dim(data), element_size)
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


#' Open an HDF5 file with basic error handling
#'
#' Wrapper around `hdf5r::H5File$new` that throws a clearer error message on
#' failure.
#'
#' @param path Path to the HDF5 file.
#' @param mode File mode passed to `H5File$new`.
#' @param ... Additional arguments forwarded to `H5File$new`.
#' @return An `H5File` object.
#' @keywords internal
open_h5 <- function(path, mode = "r", ...) {
  stopifnot(is.character(path), length(path) == 1)
  stopifnot(is.character(mode), length(mode) == 1)

  tryCatch(
    hdf5r::H5File$new(path, mode = mode, ...),
    error = function(e) {
      stop(
        sprintf("Failed to open HDF5 file '%s': %s", path, conditionMessage(e)),
        call. = FALSE
      )
    }
  )
}

#' Close an HDF5 file handle if valid
#'
#' Silently attempts to close an `H5File` handle, ignoring objects that are not
#' valid file handles.
#'
#' @param h5 Object returned by `open_h5`.
#' @return Invisible `NULL`.
#' @keywords internal
close_h5_safely <- function(h5) {
  if (inherits(h5, "H5File") && h5$is_valid()) {
    tryCatch(h5$close_all(), error = function(e) {
      warning(paste("Error closing HDF5 handle:", conditionMessage(e)))
    })
  }
  invisible(NULL)


#' Assert that an HDF5 path exists
#'
#' Convenience helper to verify that a dataset or group is present
#' at the given path relative to `h5`. Throws an `lna_error_missing_path`
#' error if the path does not exist.
#'
#' @param h5 An `H5File` or `H5Group` object.
#' @param path Character path to check.
#' @return Invisibly returns `NULL` when the path exists.
#' @keywords internal
assert_h5_path <- function(h5, path) {
  stopifnot(inherits(h5, c("H5File", "H5Group")))
  stopifnot(is.character(path), length(path) == 1)

  if (!h5$exists(path)) {
    abort_lna(
      sprintf("HDF5 path '%s' not found", path),
      .subclass = "lna_error_missing_path"
    )
  }
  invisible(NULL)
}

#' Map a datatype name to an HDF5 type
#'
#' Provides a small lookup used when creating datasets. The mapping is
#' intentionally simple but can be extended to support NIfTI conversions.
#'
#' @param dtype Character scalar naming the datatype.
#' @return An `H5T` object.
#' @keywords internal
map_dtype <- function(dtype) {
  if (inherits(dtype, "H5T")) {
    return(dtype)
  }

  stopifnot(is.character(dtype), length(dtype) == 1)

  switch(dtype,
    float32 = hdf5r::h5types$H5T_IEEE_F32LE,
    float64 = hdf5r::h5types$H5T_IEEE_F64LE,
    int8    = hdf5r::h5types$H5T_STD_I8LE,
    uint8   = hdf5r::h5types$H5T_STD_U8LE,
    int16   = hdf5r::h5types$H5T_STD_I16LE,
    uint16  = hdf5r::h5types$H5T_STD_U16LE,
    int32   = hdf5r::h5types$H5T_STD_I32LE,
    uint32  = hdf5r::h5types$H5T_STD_U32LE,
    int64   = hdf5r::h5types$H5T_STD_I64LE,
    uint64  = hdf5r::h5types$H5T_STD_U64LE,
    abort_lna(
      sprintf("Unknown dtype '%s'", dtype),
      .subclass = "lna_error_validation"
    )
  )
}

#' Guess an HDF5 datatype for an R object
#'
#' @param x R object to inspect.
#' @return An `H5T` datatype object.
#' @keywords internal
guess_h5_type <- function(x) {
  if (is.integer(x)) {
    return(map_dtype("int32"))
  } else if (is.double(x)) {
    return(map_dtype("float64"))
  } else if (is.logical(x) || is.raw(x)) {
    return(map_dtype("uint8"))
  } else if (is.character(x)) {
    t <- hdf5r::H5T_STRING$new(size = Inf)
    t$set_cset("UTF-8")
    return(t)
  }

  abort_lna(
    "Unsupported object type for HDF5 storage",
    .subclass = "lna_error_validation"
  )
}

#' Read a dataset from an HDF5 group
#'
#' @param h5_group An `H5Group` object used as the starting location for `path`.
#' @param path Character string giving the dataset path relative to `h5_group`.
#' @return The contents of the dataset.
#' @details Throws an error if the dataset does not exist or reading fails.
h5_read <- function(h5_group, path) {
  stopifnot(inherits(h5_group, "H5Group"))
  stopifnot(is.character(path), length(path) == 1)

  if (!h5_group$exists(path)) {
    stop(paste0("Dataset '", path, "' not found."), call. = FALSE)
  }

  dset <- NULL
  result <- NULL
  tryCatch({
    dset <- h5_group[[path]]
    result <- dset$read()
  }, error = function(e) {
    stop(paste0("Error reading dataset '", path, "': ", conditionMessage(e)), call. = FALSE)
  }, finally = {
    if (!is.null(dset) && inherits(dset, "H5D")) dset$close()
  })

  result
}

#' Read a subset of a dataset from an HDF5 group
#'
#' @param h5_group An `H5Group` object used as the starting location for `path`.
#' @param path Character string giving the dataset path relative to `h5_group`.
#' @param index List of indices for each dimension as accepted by `hdf5r`.
#' @return The selected subset of the dataset.
#' @details Throws an error if the dataset does not exist or reading fails.
h5_read_subset <- function(h5_group, path, index) {
  stopifnot(inherits(h5_group, "H5Group"))
  stopifnot(is.character(path), length(path) == 1)
  stopifnot(is.list(index))

  if (!h5_group$exists(path)) {
    stop(paste0("Dataset '", path, "' not found."), call. = FALSE)
  }

  dset <- NULL
  result <- NULL
  tryCatch({
    dset <- h5_group[[path]]
    result <- dset$read(args = list(index = index))
  }, error = function(e) {
    stop(paste0("Error reading subset from dataset '", path, "': ", conditionMessage(e)), call. = FALSE)
  }, finally = {
    if (!is.null(dset) && inherits(dset, "H5D")) dset$close()
  })

  result

}
