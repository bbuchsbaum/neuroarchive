#' HDF5 N-bit utilities
#'
#' Internal helpers for writing integer code datasets using the HDF5 N-bit
#' filter, plus lightweight inspection helpers.
#'
#' @import hdf5r
#' @keywords internal
NULL

.safe_h5_close_any <- function(obj) {
  if (is.null(obj)) return(invisible(NULL))
  if (inherits(obj, "H5File")) {
    return(safe_h5_close(obj))
  }
  if (inherits(obj, c("H5Group", "H5D"))) {
    return(safe_h5_close(obj))
  }
  if (is.function(obj$close)) {
    tryCatch(obj$close(), error = function(e) invisible(NULL))
  }
  invisible(NULL)
}

.h5_nbit_dtype <- function(bits) {
  bits <- as.integer(bits)
  if (!is.numeric(bits) || length(bits) != 1 || is.na(bits)) {
    abort_lna(
      "bits must be a single integer",
      .subclass = "lna_error_validation",
      location = ".h5_nbit_dtype"
    )
  }
  if (bits < 1L || bits > 16L) {
    abort_lna(
      "bits must be between 1 and 16",
      .subclass = "lna_error_validation",
      location = ".h5_nbit_dtype"
    )
  }

  bytes <- if (bits <= 8L) 1L else 2L

  dt <- h5types$H5T_NATIVE_UINT32$copy()
  dt$set_size(bytes)$set_precision(bits)$set_sign(h5const$H5T_SGN_NONE)
  dt
}

.nbit_write_options <- function() {
  shuffle <- lna_options("write.quant.nbit.shuffle")[[1]]
  if (is.null(shuffle)) shuffle <- FALSE
  deflate_level <- lna_options("write.quant.nbit.deflate_level")[[1]]
  if (is.null(deflate_level)) deflate_level <- 0L
  list(
    shuffle = isTRUE(shuffle),
    deflate_level = as.integer(deflate_level)
  )
}

.h5_nbit_dcpl <- function(chunk_dims, shuffle = FALSE, deflate_level = 0L) {
  chunk_dims <- as.integer(chunk_dims)
  if (length(chunk_dims) == 0 || any(is.na(chunk_dims)) || any(chunk_dims < 1L)) {
    abort_lna(
      "chunk_dims must be positive integers",
      .subclass = "lna_error_validation",
      location = ".h5_nbit_dcpl"
    )
  }

  dcpl <- H5P_DATASET_CREATE$new()
  dcpl$set_chunk(chunk_dims)$set_nbit()
  if (isTRUE(shuffle)) dcpl$set_shuffle()
  if (!is.null(deflate_level) && is.numeric(deflate_level) &&
      length(deflate_level) == 1 && !is.na(deflate_level) &&
      deflate_level > 0) {
    dcpl$set_deflate(as.integer(deflate_level))
  }
  dcpl
}

.ensure_h5_group <- function(h5_group, location) {
  if (!inherits(h5_group, "H5Group")) {
    abort_lna(
      "h5_group must be an H5Group object",
      .subclass = "lna_error_validation",
      location = location
    )
  }
  invisible(TRUE)
}

.split_h5_path <- function(path, location) {
  if (!is.character(path) || length(path) != 1 || !nzchar(path)) {
    abort_lna(
      "path must be a non-empty string",
      .subclass = "lna_error_validation",
      location = location
    )
  }
  parts <- strsplit(path, "/", fixed = TRUE)[[1]]
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) {
    abort_lna(
      "path must contain at least one component",
      .subclass = "lna_error_validation",
      location = location
    )
  }
  parts
}

.ensure_parent_groups <- function(h5_group, parts) {
  grp <- h5_group
  if (length(parts) > 1) {
    for (g in parts[-length(parts)]) {
      grp <- if (!grp$exists(g)) grp$create_group(g) else grp[[g]]
    }
  }
  grp
}

.overwrite_dataset_if_present <- function(grp, ds_name, location) {
  if (!grp$exists(ds_name)) return(invisible(NULL))
  existing_obj <- grp[[ds_name]]
  if (inherits(existing_obj, "H5D")) {
    existing_obj$close()
    grp$link_delete(ds_name)
    return(invisible(NULL))
  }
  existing_obj$close()
  abort_lna(
    sprintf("Object '%s' already exists and is not a dataset", ds_name),
    .subclass = "lna_error_file_exists",
    location = location
  )
}

#' List HDF5 dataset creation filters
#'
#' Returns the names of dataset filters configured on a dataset's creation
#' property list (DCPL). This includes filters such as `"nbit"`, `"shuffle"`,
#' and `"deflate"`.
#'
#' @param dset An `hdf5r::H5D` dataset.
#' @return Character vector of filter names (lowercase).
#' @export
list_h5_filters <- function(dset) {
  if (!inherits(dset, "H5D")) {
    abort_lna(
      "dset must be an hdf5r::H5D dataset",
      .subclass = "lna_error_validation",
      location = "list_h5_filters"
    )
  }

  pl <- dset$get_create_plist()
  on.exit(.safe_h5_close_any(pl), add = TRUE)

  n <- tryCatch(pl$get_nfilters(), error = function(e) 0L)
  n <- as.integer(n %||% 0L)
  if (is.na(n) || n <= 0L) return(character())

  out <- character(n)
  for (i in seq_len(n)) {
    idx0 <- as.integer(i - 1L) # hdf5r uses 0-based indices for get_filter()
    info <- tryCatch(pl$get_filter(idx0), error = function(e) NULL)
    nm <- if (is.null(info) || is.null(info$name)) NA_character_ else tolower(as.character(info$name)[1])
    out[[i]] <- nm
  }

  out <- out[!is.na(out) & nzchar(out)]
  unique(out)
}

#' Check if a dataset uses the HDF5 N-bit filter
#'
#' @param dset An `hdf5r::H5D` dataset.
#' @return Logical scalar.
#' @export
is_nbit_dataset <- function(dset) {
  if (!inherits(dset, "H5D")) return(FALSE)
  "nbit" %in% list_h5_filters(dset)
}

#' Create an empty N-bit dataset for integer codes
#'
#' Creates a chunked integer dataset with an N-bit filter configured so that
#' only `bits` significant bits are stored on disk.
#'
#' @param h5_group Starting `H5Group`.
#' @param path Dataset path relative to `h5_group` (may start with `/`).
#' @param dims Integer vector of dataset dimensions.
#' @param bits Integer number of significant bits to store (1-16).
#' @param chunk_dims Optional chunk dimensions. Defaults to an internal heuristic.
#' @return An open `H5D` dataset handle.
#' @keywords internal
h5_create_empty_dataset_nbit <- function(h5_group, path, dims, bits, chunk_dims = NULL) {
  .ensure_h5_group(h5_group, "h5_create_empty_dataset_nbit")
  parts <- .split_h5_path(path, "h5_create_empty_dataset_nbit:path")
  ds_name <- tail(parts, 1)
  grp <- .ensure_parent_groups(h5_group, parts)
  .overwrite_dataset_if_present(grp, ds_name, "h5_create_empty_dataset_nbit:overwrite")

  dims <- as.integer(dims)
  if (!is.numeric(dims) || length(dims) == 0 || any(is.na(dims)) || any(dims < 0L)) {
    abort_lna(
      "dims must be a non-empty integer vector with non-negative values",
      .subclass = "lna_error_validation",
      location = "h5_create_empty_dataset_nbit:dims"
    )
  }

  opts <- .nbit_write_options()
  dt <- .h5_nbit_dtype(bits)
  on.exit(.safe_h5_close_any(dt), add = TRUE)

  if (is.null(chunk_dims)) {
    dtype_size <- if (as.integer(bits) <= 8L) 1L else 2L
    chunk_dims <- guess_chunk_dims(dims, dtype_size)
  }
  if (length(chunk_dims) != length(dims)) {
    abort_lna(
      "chunk_dims must have the same length as dims",
      .subclass = "lna_error_validation",
      location = "h5_create_empty_dataset_nbit:chunk_dims"
    )
  }
  chunk_dims <- pmax(1L, pmin(as.integer(chunk_dims), pmax(1L, dims)))

  dcpl <- .h5_nbit_dcpl(chunk_dims, shuffle = opts$shuffle, deflate_level = opts$deflate_level)
  on.exit(.safe_h5_close_any(dcpl), add = TRUE)

  dset <- grp$create_dataset(
    name = ds_name,
    dims = dims,
    dtype = dt,
    chunk_dims = NULL,
    gzip_level = NULL,
    dataset_create_pl = dcpl
  )

  # Lightweight metadata for downstream consumers.
  h5_attr_write(dset, "quant_bits", as.integer(bits))
  dset
}

#' Write integer codes to an N-bit dataset
#'
#' Materialises a dataset with N-bit storage enabled and writes `codes`.
#'
#' @param h5_group Starting `H5Group`.
#' @param path Dataset path relative to `h5_group` (may start with `/`).
#' @param codes Integer vector/matrix/array of codes.
#' @param bits Integer number of significant bits to store (1-16).
#' @param chunk_dims Optional chunk dimensions. Defaults to an internal heuristic.
#' @return An open `H5D` dataset handle.
#' @keywords internal
h5_write_qbits_nbit <- function(h5_group, path, codes, bits, chunk_dims = NULL) {
  .ensure_h5_group(h5_group, "h5_write_qbits_nbit")
  parts <- .split_h5_path(path, "h5_write_qbits_nbit:path")
  ds_name <- tail(parts, 1)
  grp <- .ensure_parent_groups(h5_group, parts)
  .overwrite_dataset_if_present(grp, ds_name, "h5_write_qbits_nbit:overwrite")

  if (!is.array(codes)) {
    if (is.vector(codes)) {
      dim(codes) <- length(codes)
    } else {
      abort_lna(
        "codes must be an array/matrix or atomic vector",
        .subclass = "lna_error_validation",
        location = "h5_write_qbits_nbit:codes"
      )
    }
  }

  dims <- dim(codes)
  if (is.null(dims) || length(dims) == 0) {
    abort_lna(
      "Unable to determine dataset dimensions from codes",
      .subclass = "lna_error_validation",
      location = "h5_write_qbits_nbit:dims"
    )
  }
  codes <- array(as.integer(codes), dim = dims)

  maxv <- bitwShiftL(1L, as.integer(bits)) - 1L
  if (any(codes < 0L, na.rm = TRUE) || any(codes > maxv, na.rm = TRUE)) {
    abort_lna(
      sprintf("codes must be in [0, %d] for bits=%d", maxv, as.integer(bits)),
      .subclass = "lna_error_validation",
      location = "h5_write_qbits_nbit:range"
    )
  }

  opts <- .nbit_write_options()
  dt <- .h5_nbit_dtype(bits)
  on.exit(.safe_h5_close_any(dt), add = TRUE)

  if (is.null(chunk_dims)) {
    dtype_size <- if (as.integer(bits) <= 8L) 1L else 2L
    chunk_dims <- guess_chunk_dims(as.integer(dims), dtype_size)
  }
  if (length(chunk_dims) != length(dims)) {
    abort_lna(
      "chunk_dims must have the same length as dims",
      .subclass = "lna_error_validation",
      location = "h5_write_qbits_nbit:chunk_dims"
    )
  }
  chunk_dims <- pmax(1L, pmin(as.integer(chunk_dims), pmax(1L, as.integer(dims))))

  dcpl <- .h5_nbit_dcpl(chunk_dims, shuffle = opts$shuffle, deflate_level = opts$deflate_level)
  on.exit(.safe_h5_close_any(dcpl), add = TRUE)

  dset <- grp$create_dataset(
    name = ds_name,
    robj = codes,
    dtype = dt,
    dims = as.integer(dims),
    chunk_dims = NULL,
    gzip_level = NULL,
    dataset_create_pl = dcpl
  )
  h5_attr_write(dset, "quant_bits", as.integer(bits))
  dset
}
