
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
#' @return Invisibly returns the `H5File` handle. When `checksum = "sha256"` the
#'   handle is closed as part of digest computation and the returned object will
#'   no longer be valid.
#' @import hdf5r
#' @keywords internal
materialise_plan <- function(h5, plan, checksum = c("none", "sha256"),
                             header = NULL, plugins = NULL) {
  checksum <- match.arg(checksum)
  stopifnot(inherits(h5, "H5File"))
  if (!h5$is_valid()) {
    stop("`h5` handle is not valid", call. = FALSE)
  }
  stopifnot(inherits(plan, "Plan"))
  if (!is.null(header)) {
    stopifnot(is.list(header))
  }
  if (!is.null(plugins)) {
    stopifnot(is.list(plugins))
    if (is.null(names(plugins)) || any(names(plugins) == "")) {
      stop("plugins must be a named list", call. = FALSE)
    }
  }

  # Create core groups
  if (h5$exists("transforms")) {
    obj <- h5[["transforms"]]
    if (!inherits(obj, "H5Group")) {
      obj$close()
      stop("'/transforms' exists but is not a group", call. = FALSE)
    }
    tf_group <- obj
  } else {
    tf_group <- h5$create_group("transforms")
  }

  if (h5$exists("basis")) {
    obj <- h5[["basis"]]
    if (!inherits(obj, "H5Group")) {
      obj$close()
      stop("'/basis' exists but is not a group", call. = FALSE)
    }
    obj$close()
  } else {
    h5$create_group("basis")
  }

  if (h5$exists("scans")) {
    obj <- h5[["scans"]]
    if (!inherits(obj, "H5Group")) {
      obj$close()
      stop("'/scans' exists but is not a group", call. = FALSE)
    }
    obj$close()
  } else {
    h5$create_group("scans")
  }

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
  write_payload <- function(path, data, step_index) {
    comp_level <- lna_options("write.compression_level")[[1]]
    if (is.null(comp_level)) comp_level <- 0
    chunk_dims <- NULL

    attempt <- function(level, chunks) {
      ds <- h5_write_dataset(root, path, data, chunk_dims = chunks,
                             compression_level = level)
      if (inherits(ds, "H5D")) ds$close()
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

    if (is.integer(data)) {
      dtype_size <- 4L
    } else if (is.double(data)) {
      dtype_size <- 8L
    } else {
      t <- guess_h5_type(data)
      dtype_size <- t$get_size()
      if (inherits(t, "H5T")) t$close()
    }
    cdims <- if (is.null(chunk_dims)) guess_chunk_dims(dim(data), dtype_size) else as.integer(chunk_dims)

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
    progress_enabled <- steps > 1 && !progressr::handlers_is_empty()
    loop <- function() {
      p <- if (progress_enabled) progressr::progressor(steps = steps) else NULL
      for (i in idx) {
        row <- plan$datasets[i, ]
        key <- row$payload_key
        if (!nzchar(key)) next
        payload <- plan$payloads[[key]]
        if (is.null(payload)) {
          warning(sprintf("Payload '%s' missing for dataset '%s'", key, row$path))
          next
        }
        if (!is.null(p)) p(message = row$path)
        write_payload(row$path, payload, row$step_index)
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

  if (!is.null(header) && length(header) > 0) {
    hdr_grp <- if (!h5$exists("header")) h5$create_group("header") else h5[["header"]]
    g <- if (hdr_grp$exists("global")) hdr_grp[["global"]] else hdr_grp$create_group("global")
    for (nm in names(header)) {
      h5_attr_write(g, nm, header[[nm]])
    }
  }

  if (!is.null(plugins) && length(plugins) > 0) {
    pl_grp <- if (!h5$exists("plugins")) h5$create_group("plugins") else h5[["plugins"]]
    for (nm in names(plugins)) {
      if (grepl("/", nm, fixed = TRUE)) {
        stop(sprintf("Invalid plugin name '%s'", nm), call. = FALSE)
      }
      write_json_descriptor(pl_grp, paste0(nm, ".json"), plugins[[nm]])
    }
  }

  if (checksum == "sha256") {
    file_path <- h5$filename
    close_h5_safely(h5)
    if (is.character(file_path) && nzchar(file_path) && file.exists(file_path)) {
      hash_val <- digest::digest(file = file_path, algo = "sha256")
      h5_tmp <- open_h5(file_path, mode = "r+")
      root_tmp <- h5_tmp[["/"]]
      h5_attr_write(root_tmp, "lna_checksum", hash_val)
      close_h5_safely(h5_tmp)
    } else {
      warning("Checksum requested but file path unavailable; skipping")
    }
  }

  invisible(h5)
}
