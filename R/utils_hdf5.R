#' HDF5 Attribute Read/Write Helpers
#'
#' @description Provides internal functions for reading, writing, checking existence,
#'   and deleting attributes associated with HDF5 objects (groups or datasets).
#'
#' @import hdf5r
#' @importFrom hdf5r H5P_DATASET_XFER H5P_FILE_CREATE H5P_DEFAULT
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
  assert_scalar_character(name, "name")

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
  assert_scalar_character(name, "name")

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
  assert_scalar_character(name, "name")

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
  assert_scalar_character(name, "name")

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
  target_mib <- lna_options("write.chunk_target_mib")[[1]]
  if (is.null(target_mib)) target_mib <- 1
  target <- as.numeric(target_mib) * 1024^2
  
  # Ensure dims is a valid integer vector
  dims <- as.integer(dims)
  if (length(dims) == 0) {
    stop("Invalid dimensions for chunk calculation")
  }
  
  # Handle empty dimensions - any dimension can be 0 for empty arrays
  if (any(dims == 0)) {
    # For empty arrays, return safe chunk dimensions (1s with same length)
    return(rep(1L, length(dims)))
  }
  
  # Now check for invalid negative dimensions
  if (any(dims < 1)) {
    stop("Invalid dimensions for chunk calculation")
  }
  
  chunk <- tryCatch({
    hdf5r::guess_chunks(space_maxdims = dims,
                        dtype_size = dtype_size,
                        chunk_size = target)
  }, error = function(e) {
    # Fallback: use the dims themselves, but capped to reasonable values
    pmin(dims, 1024L)
  })
  
  # Ensure chunk is valid and has the right length
  if (length(chunk) == 0 || length(chunk) != length(dims) || any(is.na(chunk)) || any(chunk < 1)) {
    # Fallback: use the dims themselves, but capped to reasonable values
    chunk <- pmin(dims, 1024L)
  }

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
  chunk <- pmax(chunk, 1L)
  
  # Final validation
  if (length(chunk) != length(dims)) {
    chunk <- pmin(dims, 1024L)
  }
  
  return(chunk)
}

#' Reduce chunk dimensions toward a byte target
#'
#' Helper used when retrying dataset writes. Starting from an existing
#' chunk dimension vector, halves the first dimension until the
#' estimated chunk size is below `target_bytes` or the dimension would
#' drop below 1. Returns the adjusted chunk vector.
#'
#' @param chunk Integer vector of current chunk dimensions.
#' @param dtype_size Size in bytes of the datatype being stored.
#' @param target_bytes Target maximum chunk size in bytes.
#' @return Integer vector of reduced chunk dimensions.
#' @keywords internal
reduce_chunk_dims <- function(chunk, dtype_size, target_bytes) {
  stopifnot(is.numeric(chunk))
  chunk <- as.integer(chunk)
  chunk_bytes <- prod(chunk) * dtype_size
  while (chunk_bytes > target_bytes && chunk[1] > 1) {
    chunk[1] <- ceiling(chunk[1] / 2)
    chunk_bytes <- prod(chunk) * dtype_size
  }
  chunk
}

#' Create an empty HDF5 dataset
#'
#' Helper used when block-wise algorithms need a dataset skeleton to write
#' slabs into. This mirrors \code{h5_write_dataset} but allocates the dataset
#' without supplying data.
#'
#' @param h5_group An `H5Group` object used as the starting location for `path`.
#' @param path Character string giving the dataset path relative to `h5_group`.
#' @param dims Integer vector of dataset dimensions.
#' @param dtype Character string naming the datatype (e.g. "uint8", "float32").
#' @param chunk_dims Optional integer vector specifying chunk layout. When
#'   `NULL`, [guess_chunk_dims()] is used.
#' @return Invisibly returns `TRUE` on success.
#' @keywords internal
h5_create_empty_dataset <- function(h5_group, path, dims, dtype,
                                    chunk_dims = NULL) {
  assert_h5group(h5_group)
  assert_scalar_character(path, "path")
  stopifnot(is.numeric(dims))
  stopifnot(is.character(dtype) || inherits(dtype, "H5T"))

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
    size <- map_dtype(dtype)$get_size(variable_as_inf = FALSE)
    chunk_dims <- guess_chunk_dims(as.integer(dims), size)
  } else {
    chunk_dims <- as.integer(chunk_dims)
  }

  dty <- map_dtype(dtype)
  ds <- grp$create_dataset(ds_name,
                           dims = as.integer(dims),
                           dtype = dty,
                           chunk_dims = chunk_dims)
  # if (inherits(ds, "H5D")) ds$close() # Do not close the dataset here
  invisible(ds) # Return the created dataset object
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
#' @return Invisibly returns `TRUE` on success.
h5_write_dataset <- function(h5_group, path, data,
                             chunk_dims = NULL, compression_level = 0,
                             dtype = NULL) {
  assert_h5group(h5_group)
  assert_scalar_character(path, "path")
  stopifnot(is.numeric(compression_level), length(compression_level) == 1)

  if (!is.array(data)) {
    if (is.vector(data)) {
      dim(data) <- length(data)
    } else {
      stop("`data` must be a matrix or array")
    }
  }
  
  data_dims <- dim(data) # These are the dimensions of the incoming data
  if (is.null(data_dims) || length(data_dims) == 0) {
    stop("Unable to determine valid dimensions for data")
  }
  
  is_empty_array <- any(data_dims == 0)

  if (is_empty_array) {
    placeholder_data <- array(as.integer(data_dims), dim = c(length(data_dims), 1))
    
    data_to_write_final <- placeholder_data
    dtype_final <- map_dtype("int32") 
  } else {
    data_to_write_final <- data 
    dtype_final <- if (!is.null(dtype)) map_dtype(dtype) else guess_h5_type(data)
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
    element_size <- dtype_final$get_size(variable_as_inf = FALSE)
    # If writing a placeholder, chunks are based on placeholder dims.
    # Otherwise, based on original data dims.
    effective_dims_for_chunking <- if(is_empty_array) dim(data_to_write_final) else data_dims
    chunk_dims <- guess_chunk_dims(effective_dims_for_chunking, element_size)
  } else {
    chunk_dims <- as.integer(chunk_dims)
  }

  create_fun <- function(level) {
    on.exit(if (inherits(dtype_final, "H5T") && is.null(dtype) && !is_empty_array) dtype_final$close(), add = TRUE)
    
    current_data_to_write <- data_to_write_final
    if (!is_empty_array && is.raw(current_data_to_write) && !is.null(dtype) && dtype == "uint8") {
      current_data_to_write <- as.integer(current_data_to_write)
      dim(current_data_to_write) <- dim(data_to_write_final) 
    }
    
    if (grp$exists(ds_name)) {
      existing_obj <- grp[[ds_name]]
      if (inherits(existing_obj, "H5D")) {
        existing_obj$close()
        grp$link_delete(ds_name)
      } else {
        existing_obj$close()
        stop(sprintf("Object '%s' already exists and is not a dataset", ds_name))
      }
    }
    
    created_dset <- grp$create_dataset(ds_name,
                                       robj = current_data_to_write, 
                                       chunk_dims = chunk_dims,
                                       dtype = dtype_final, 
                                       gzip_level = level)
    
    if (is_empty_array) {
      h5_attr_write(created_dset, "lna_empty_array_placeholder", TRUE)
    }
    created_dset 
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

  if (inherits(dset, "H5D")) dset$close()
  invisible(TRUE)
}

#' @importFrom hdf5r H5File
#' Open an HDF5 file with basic error handling
#'
#' Wrapper around `hdf5r::H5File$new` that throws a clearer error message on
#' failure.
#'
#' @param path Path to the HDF5 file.
#' @param mode File mode passed to `H5File$new`.
#' @return An `H5File` object.
#' @details When `mode` is `"w"` the file is truncated if it already
#'   exists. Use a unique temporary file and `file.rename()` when
#'   writing in parallel.
#' @keywords internal
open_h5 <- function(path, mode = "a") {
  assert_scalar_character(path, "path")
  assert_scalar_character(mode, "mode")

  tryCatch(
    hdf5r::H5File$new(path, mode = mode),
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
  if (inherits(h5, "H5File") && h5$is_valid) {
    tryCatch(h5$close_all(), error = function(e) {
      warning(paste("Error closing HDF5 handle:", conditionMessage(e)))
    })
  }
  invisible(NULL)
}

#' Safely check for the existence of an HDF5 path
#'
#' Wrapper around `H5Group$exists` that treats any error as the path not
#' existing. This is useful when paths may contain characters that `exists()`
#' cannot handle cleanly.
#'
#' @param group An `H5File` or `H5Group` to check.
#' @param path_name Character path to test.
#' @return Logical `TRUE` if the path exists, otherwise `FALSE`.
#' @keywords internal
path_exists_safely <- function(group, path_name) {
  if (is.null(path_name) || !nzchar(path_name)) return(FALSE)
  tryCatch({
    group$exists(path_name)
  }, error = function(e) {
    FALSE
  })
}

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
  assert_scalar_character(path, "path")

  if (!h5$exists(path)) {
    abort_lna(
      sprintf("HDF5 path '%s' not found", path),
      .subclass = "lna_error_missing_path",
      location = sprintf("assert_h5_path:%s", path)
    )
  }
  invisible(NULL)
}


#' Safely check if an HDF5 path exists
#'
#' Wrapper around `$exists` that catches errors (e.g., invalid paths)
#' and returns `FALSE` instead of propagating the error.
#'
#' @param h5 An `H5File` or `H5Group` object.
#' @param path Character scalar dataset or group path.
#' @return Logical scalar, `TRUE` if the path exists, `FALSE` otherwise.
#' @keywords internal
path_exists_safely <- function(h5, path) {
  if (is.null(path) || !nzchar(path)) return(FALSE)
  stopifnot(inherits(h5, c("H5File", "H5Group")))
  assert_scalar_character(path, "path")


  tryCatch({
    h5$exists(path)
  }, error = function(e) FALSE)
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

  assert_scalar_character(dtype, "dtype")

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
      .subclass = "lna_error_validation",
      location = "map_dtype"
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
    .subclass = "lna_error_validation",
    location = "guess_h5_type"
  )
}

#' Read a dataset from an HDF5 group
#'
#' @param h5_group An `H5Group` object used as the starting location for `path`.
#' @param path Character string giving the dataset path relative to `h5_group`.
#' @return The contents of the dataset.
#' @details Throws an error if the dataset does not exist or reading fails.
h5_read <- function(h5_group, path) {
  assert_h5group(h5_group)
  assert_scalar_character(path, "path")

  if (!h5_group$exists(path)) {
    stop(paste0("Dataset '", path, "' not found."), call. = FALSE)
  }

  dset <- NULL
  result <- NULL
  tryCatch({
    dset <- h5_group[[path]]
    # If not a placeholder, read the dataset normally
    result <- dset$read(drop = FALSE) 
    
    # Check for the empty array placeholder attribute
    is_placeholder <- FALSE
    if (h5_attr_exists(dset, "lna_empty_array_placeholder")) {
      is_placeholder <- isTRUE(h5_attr_read(dset, "lna_empty_array_placeholder"))
    }

    # If it's a placeholder, reconstruct the original empty structure
    if (is_placeholder) {
      # Check if placeholder is a vector (hdf5r might convert Nx1 matrix to vector)
      if (is.vector(result) && !is.list(result) && is.numeric(result) && 
          all(result >= 0) && all(result == as.integer(result))) {
        orig_dims <- as.integer(result)
        result <- array(numeric(0), dim = orig_dims)
      } else if (is.matrix(result) && ncol(result) == 1 && is.numeric(result) && 
          all(result >= 0) && all(result == as.integer(result))) {
        orig_dims <- as.integer(result[,1])
        result <- array(numeric(0), dim = orig_dims)
      } else {
        # This case should ideally not happen if writing was correct
        warning(sprintf("Dataset '%s' marked as empty placeholder but data format is unexpected. Result class: %s, dim: %s, values: %s. Returning as is.", 
                        path, class(result)[1], paste(dim(result), collapse="x"), paste(result, collapse=",")))
      }
    }
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
  assert_h5group(h5_group)
  assert_scalar_character(path, "path")
  stopifnot(is.list(index))

  if (!h5_group$exists(path)) {
    stop(paste0("Dataset '", path, "' not found."), call. = FALSE)
  }

  dset <- NULL
  result <- NULL
  tryCatch({
    dset <- h5_group[[path]]
    result <- dset$read(args = index)
  }, error = function(e) {
    stop(paste0("Error reading subset from dataset '", path, "': ", conditionMessage(e)), call. = FALSE)
  }, finally = {
    if (!is.null(dset) && inherits(dset, "H5D")) dset$close()
  })

  result
}

#' Discover run identifiers in an LNA file
#'
#' Lists available run groups under `/scans` that match the `run-XX` pattern.
#'
#' @param h5 An open `H5File` object.
#' @return Character vector of run identifiers sorted alphabetically.
#' @keywords internal
discover_run_ids <- function(h5) {
  stopifnot(inherits(h5, "H5File"))
  if (!h5$exists("scans")) {
    return(character())
  }
  grp <- h5[["scans"]]
  nms <- grp$ls()$name
  runs <- grep("^run-", nms, value = TRUE)
  sort(runs)
}

#' Resolve run_id patterns against available runs
#'
#' @param patterns Character vector of run_id patterns or names. `NULL` selects the first available run.
#' @param available Character vector of available run identifiers.
#' @return Character vector of matched run identifiers.
#' @keywords internal
resolve_run_ids <- function(patterns, available) {
  if (is.null(patterns)) {
    return(if (length(available) > 0) available[1] else character())
  }
  patterns <- as.character(patterns)
  out <- character()
  for (p in patterns) {
    if (grepl("[*?]", p)) {
      rx <- utils::glob2rx(p)
      out <- union(out, available[grepl(rx, available)])
    } else {
      out <- union(out, intersect(available, p))
    }
  }
  unique(out)
}

#' Validate and sanitize run identifiers
#'
#' Ensures that \code{run_id} matches the expected ``"run-XX"`` pattern and
#' does not contain path separators.  Returns the sanitized identifier or
#' throws an error on invalid input.
#'
#' @param run_id Character scalar run identifier.
#' @return The validated \code{run_id} string.
#' @keywords internal
sanitize_run_id <- function(run_id) {
  assert_scalar_character(run_id, "run_id")
  if (grepl("/|\\\\", run_id)) {
    abort_lna(
      "run_id must not contain path separators",
      .subclass = "lna_error_validation",
      location = "sanitize_run_id"
    )
  }
  if (!grepl("^run-[0-9]{2}$", run_id)) {
    abort_lna(
      "run_id must match 'run-XX' pattern",
      .subclass = "lna_error_validation",
      location = "sanitize_run_id"
    )
  }
  run_id
}
