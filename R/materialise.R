#' Materialise Plan to HDF5
#'
#' @description Writes transform descriptors and payload datasets to an open
#'   HDF5 file according to the provided `plan`. Implements basic retry logic
#'   for common HDF5 errors.
#' @param h5 An open `H5File` object.
#' @param plan A `Plan` R6 object produced by `core_write`.
#' @param checksum Character string indicating checksum mode.
#' @param header Optional named list of header attributes.
#' @param plugins Optional named list of plugin metadata.

#' @return Invisibly returns the `H5File` handle. When `checksum = "sha256"`
#'   the file is first written with a placeholder checksum attribute, the
#'   SHA256 digest is computed on that file, and then the attribute is updated
#'   with the final value.  The handle is closed during digest calculation and
#'   is therefore invalid when the function returns.
#' @import hdf5r
#' @keywords internal
materialise_plan <- function(h5, plan, checksum = c("none", "sha256"),
                             header = NULL, plugins = NULL) {
  checksum <- match.arg(checksum)
  stopifnot(inherits(h5, "H5File"))
  if (!h5$is_valid) {

    abort_lna(
      "Provided HDF5 handle is not open or valid",
      .subclass = "lna_error_validation",
      location = "materialise_plan:h5"
    )
  }
  stopifnot(inherits(plan, "Plan"))

  header <- validate_named_list(header, "header")
  plugins <- validate_named_list(plugins, "plugins")

  # Create core groups
  tf_group <- get_or_create_group(
    h5, "transforms", "materialise_plan:transforms", return_handle = TRUE
  )
  get_or_create_group(h5, "basis", "materialise_plan:basis", return_handle = FALSE)
  get_or_create_group(h5, "scans", "materialise_plan:scans", return_handle = FALSE)

  root <- h5[["/"]]
  h5_attr_write(root, "lna_spec", "LNA R v2.0")
  h5_attr_write(root, "creator", "lna R package v0.0.1")
  h5_attr_write(root, "required_transforms", character(0))

  # Write descriptors
  if (length(plan$descriptors) > 0) {
    for (nm in names(plan$descriptors)) {
      write_json_descriptor(tf_group, nm, plan$descriptors[[nm]])
    }
  }

  # Helper to write a single payload dataset with retries
  write_payload <- function(path, data, step_index, dtype_str = NA_character_) {
    comp_level <- lna_options("write.compression_level")[[1]]
    if (is.null(comp_level)) comp_level <- 0
    chunk_dims <- NULL

    attempt <- function(level, chunks) {
      effective_dtype <- if (is.na(dtype_str)) NULL else dtype_str
      h5_write_dataset(root, path, data, chunk_dims = chunks,
                       compression_level = level, dtype = effective_dtype)
      NULL
    }

    res <- tryCatch(attempt(comp_level, chunk_dims), error = function(e) e)
    if (inherits(res, "error")) {
      msg1 <- conditionMessage(res)
      if (!is.null(comp_level) && comp_level > 0 &&
          grepl("filter", msg1, ignore.case = TRUE)) {
        warning(sprintf("Compression failed for %s; retrying without compression", path))
        res <- tryCatch(attempt(0, chunk_dims), error = function(e) e)
      }
    }

    # Get data dimensions, handling vectors that don't have dim() set
    data_dims <- dim(data)
    if (is.null(data_dims)) {
      if (is.vector(data)) {
        data_dims <- length(data)
      } else {
        stop("Unable to determine dimensions for data")
      }
    }

    # Determine datatype size for chunk heuristics
    dtype <- if (!is.na(dtype_str)) map_dtype(dtype_str) else guess_h5_type(data)
    dtype_size <- dtype$get_size(variable_as_inf = FALSE)
    if (!is.finite(dtype_size) || dtype_size <= 0) {
      dtype_size <- 1L
    }
    if (inherits(dtype, "H5T") && is.na(dtype_str)) dtype$close()
    cdims <- if (is.null(chunk_dims)) guess_chunk_dims(data_dims, dtype_size) else as.integer(chunk_dims)

    if (inherits(res, "error")) {
      cdims1 <- reduce_chunk_dims(cdims, dtype_size, 1024^3)
      warning(sprintf(
        "Write failed for %s; retrying with smaller chunks (<1 GiB, ~%.1f MiB)",
        path, prod(cdims1) * dtype_size / 1024^2
      ))
      res <- tryCatch(attempt(0, cdims1), error = function(e) e)
    }

    if (inherits(res, "error")) {
      cdims2 <- reduce_chunk_dims(cdims1, dtype_size, 256 * 1024^2)
      warning(sprintf(
        "Write failed for %s; retrying with smaller chunks (<=256 MiB, ~%.1f MiB)",
        path, prod(cdims2) * dtype_size / 1024^2
      ))
      res <- tryCatch(attempt(0, cdims2), error = function(e) e)
    }

    if (inherits(res, "error")) {
      abort_lna(
        sprintf(
          "Failed to write dataset '%s' (step %d): %s",
          path, step_index, conditionMessage(res)
        ),
        .subclass = "lna_error_hdf5_write",
        location = sprintf("materialise_plan[%d]:%s", step_index, path),
        parent = res
      )
    }
  }

  # Write payload datasets
  if (nrow(plan$datasets) > 0) {
    idx <- seq_len(nrow(plan$datasets))
    has_payload <- plan$datasets$payload_key != "" & !is.na(plan$datasets$payload_key)
    steps <- sum(has_payload)
    progress_enabled <- steps > 1 && is_progress_globally_enabled()
    loop <- function() {
      p <- if (progress_enabled) progressr::progressor(steps = steps) else NULL
      for (i in idx) {
        row <- plan$datasets[i, ]
        key <- row$payload_key
        if (!nzchar(key)) next
        payload <- plan$payloads[[key]]
        if (is.null(payload)) {
          warning(sprintf("Payload '%s' missing; skipping dataset %s", key, row$path))
          next
        }
        if (!is.null(p)) p(message = row$path)

        write_payload(row$path, payload, row$step_index, row$dtype)
        if (row$producer == "quant" && row$role == "quantized") {
          bits_val <- tryCatch(jsonlite::fromJSON(row$params_json)$bits,
                               error = function(e) NULL)
          bits_val <- bits_val %||% 8L
          dset_obj <- root[[row$path]]
          h5_attr_write(dset_obj, "quant_bits", as.integer(bits_val))
          dset_obj$close()
        }
        #write_payload(row$path, payload, row$step_index, row$dtype)
        plan$datasets$write_mode_effective[i] <- "eager"
        plan$mark_payload_written(key)
      }
    }
    if (progress_enabled) {
      progressr::with_progress(loop())
    } else {
      loop()
    }
  }

  write_header_section(h5, header)
  write_plugins_section(h5, plugins)

  if (checksum == "sha256") {
    write_sha256_checksum(h5, root)
  }

  invisible(h5)
}

# Helper utilities -------------------------------------------------------

get_or_create_group <- function(h5, name, location, return_handle = TRUE) {
  if (h5$exists(name)) {
    obj <- h5[[name]]
    if (!inherits(obj, "H5Group")) {
      obj$close()
      abort_lna(
        sprintf("'/%s' already exists and is not a group", name),
        .subclass = "lna_error_validation",
        location = location
      )
    }
  } else {
    obj <- h5$create_group(name)
  }

  if (return_handle) {
    obj
  } else {
    obj$close()
    invisible(NULL)
  }
}

write_header_section <- function(h5, header) {
  if (is.null(header) || length(header) == 0) {
    return(invisible(NULL))
  }
  hdr_grp <- if (!h5$exists("header")) h5$create_group("header") else h5[["header"]]
  g <- if (hdr_grp$exists("global")) hdr_grp[["global"]] else hdr_grp$create_group("global")
  for (nm in names(header)) {
    h5_attr_write(g, nm, header[[nm]])
  }
  g$close()
  hdr_grp$close()
  invisible(NULL)
}

write_plugins_section <- function(h5, plugins) {
  if (is.null(plugins) || length(plugins) == 0) {
    return(invisible(NULL))
  }
  pl_grp <- if (!h5$exists("plugins")) h5$create_group("plugins") else h5[["plugins"]]
  for (nm in names(plugins)) {
    if (grepl("/", nm)) {
      abort_lna(
        sprintf("Plugin name '%s' contains '/' which is not allowed", nm),
        .subclass = "lna_error_validation",
        location = sprintf("materialise_plan:plugin[%s]", nm)
      )
    }
    write_json_descriptor(pl_grp, paste0(nm, ".json"), plugins[[nm]])
  }
  pl_grp$close()
  invisible(NULL)
}

write_sha256_checksum <- function(h5, root) {
  placeholder <- paste(rep("0", 64), collapse = "")
  h5_attr_write(root, "lna_checksum", placeholder)

  file_path <- h5$filename
  neuroarchive:::close_h5_safely(h5)

  if (is.character(file_path) && nzchar(file_path) && file.exists(file_path)) {
    hash_val <- digest::digest(file = file_path, algo = "sha256")
    h5_tmp <- NULL
    tryCatch({
      h5_tmp <- open_h5(file_path, mode = "r+")
      root_tmp <- h5_tmp[["/"]]
      h5_attr_write(root_tmp, "lna_checksum", hash_val)
    }, finally = {
      if (!is.null(h5_tmp) && inherits(h5_tmp, "H5File") && h5_tmp$is_valid) {
        neuroarchive:::close_h5_safely(h5_tmp)
      }
    })
  } else {
    warning(
      "Checksum requested but file path unavailable or invalid; skipping write of checksum attribute.")
  }
  invisible(NULL)
}
