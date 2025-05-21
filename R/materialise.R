
#' Materialise Plan to HDF5
#'
#' @description Writes transform descriptors and payload datasets to an open
#'   HDF5 file according to the provided `plan`. Implements basic retry logic
#'   for common HDF5 errors.
#' @param h5 An open `H5File` object.
#' @param plan A `Plan` R6 object produced by `core_write`.
#' @return Invisibly returns the modified `plan`.
#' @import hdf5r
#' @keywords internal
materialise_plan <- function(h5, plan, checksum = c("none", "sha256")) {
  checksum <- match.arg(checksum)
  stopifnot(inherits(h5, "H5File"))
  stopifnot(inherits(plan, "Plan"))

  # Create core groups
  tf_group <- if (h5$exists("transforms")) h5[["transforms"]] else h5$create_group("transforms")
  if (!h5$exists("basis")) h5$create_group("basis")
  if (!h5$exists("scans")) h5$create_group("scans")

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
    comp_level <- lna_options("write.compression")[[1]]
    if (is.null(comp_level)) comp_level <- 0
    chunk_dims <- NULL

    attempt <- function(level, chunks) {
      h5_write_dataset(root, path, data, chunk_dims = chunks,
                       compression_level = level)
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

    if (inherits(res, "error")) {
      warning(sprintf("Write failed for %s; retrying with smaller chunks", path))
      cdims <- if (is.null(chunk_dims)) pmin(dim(data), rep(128L, length(dim(data)))) else chunk_dims
      cdims[length(cdims)] <- pmax(1L, floor(cdims[length(cdims)] / 2))
      res <- tryCatch(attempt(0, cdims), error = function(e) e)
    }

    if (inherits(res, "error")) {
      abort_lna(sprintf("Failed to write dataset '%s' (step %d): %s",
                        path, step_index, conditionMessage(res)),
                .subclass = "lna_error_hdf5_write")
    }
  }

  # Write payload datasets
  if (nrow(plan$datasets) > 0) {
    for (i in seq_len(nrow(plan$datasets))) {
      row <- plan$datasets[i, ]
      key <- row$payload_key
      if (!nzchar(key)) next
      payload <- plan$payloads[[key]]
      if (is.null(payload)) next

      write_payload(row$path, payload, row$step_index)

      plan$datasets$write_mode_effective[i] <- "eager"
      plan$mark_payload_written(key)
    }
  }

  if (checksum == "sha256") {
    file_path <- h5$filename
    h5$close_all()
    if (is.character(file_path) && nzchar(file_path) && file.exists(file_path)) {
      hash_val <- digest::digest(file = file_path, algo = "sha256")
      h5_tmp <- hdf5r::H5File$new(file_path, mode = "r+")
      root_tmp <- h5_tmp[["/"]]
      h5_attr_write(root_tmp, "lna_checksum", hash_val)
      h5_tmp$close_all()
    } else {
      warning("Checksum requested but file path unavailable; skipping")
    }
  }

  invisible(h5)
}
