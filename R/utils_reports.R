#' Retrieve a transform report from an LNA file
#'
#' Opens the HDF5 file, reads the descriptor for the chosen transform and
#' returns the parsed JSON report associated with that transform.
#'
#' @param lna_file Path to an LNA file.
#' @param transform_index_or_name Integer index (0-based) or descriptor name
#'   identifying the transform.
#'
#' @return A list parsed from the JSON report.
#' @examples
#' tmp <- tempfile(fileext = ".h5")
#' arr <- array(runif(6), dim = c(2,3))
#' write_lna(arr, tmp, transforms = "quant")
#' rep <- lna_get_transform_report(tmp, 0)
#' @export
lna_get_transform_report <- function(lna_file, transform_index_or_name) {
  stopifnot(is.character(lna_file), length(lna_file) == 1)

  h5 <- open_h5(lna_file, mode = "r")
  on.exit(close_h5_safely(h5))

  tf_group <- h5[["transforms"]]
  meta <- discover_transforms(tf_group)

  desc_name <- NULL
  if (is.numeric(transform_index_or_name)) {
    idx <- as.integer(transform_index_or_name)
    row <- meta[meta$index == idx, , drop = FALSE]
    if (nrow(row) == 0) {
      stop(sprintf("Transform index %s not found", idx), call. = FALSE)
    }
    desc_name <- row$name[[1]]
  } else {
    nm <- as.character(transform_index_or_name)
    if (nm %in% meta$name) {
      desc_name <- nm
    } else {
      base <- tools::file_path_sans_ext(nm)
      row <- meta[tools::file_path_sans_ext(meta$name) == base, , drop = FALSE]
      if (nrow(row) == 0) {
        stop(sprintf("Transform '%s' not found", nm), call. = FALSE)
      }
      desc_name <- row$name[[1]]
    }
  }

  desc <- read_json_descriptor(tf_group, desc_name)
  params <- desc$params %||% list()
  base <- tools::file_path_sans_ext(desc_name)
  report_path <- params$report_path %||% paste0("/transforms/", base, "_report.json")

  root <- h5[["/"]]
  assert_h5_path(root, report_path)

  dset <- root[[report_path]]
  on.exit(if (!is.null(dset) && inherits(dset, "H5D")) dset$close(), add = TRUE)
  raw_data <- dset$read()
  
  # Handle case where raw data was stored as integers for hdf5r compatibility
  if (is.integer(raw_data) && h5_attr_exists(dset, "compression")) {
    raw_data <- as.raw(raw_data)
  }
  
  if (h5_attr_exists(dset, "compression")) {
    comp <- h5_attr_read(dset, "compression")
    if (identical(comp, "gzip")) {
      raw_data <- memDecompress(raw_data, type = "gzip")
    }
  }

  # Ensure json_str is a character string for jsonlite::fromJSON
  # If raw_data is raw (either originally or after decompression), convert to char.
  # If it was already character, it remains character.
  json_str <- if (is.raw(raw_data)) {
    rawToChar(raw_data)
  } else {
    as.character(raw_data)
  }
  
  # Remove potential leading/trailing whitespace
  json_str <- trimws(json_str)
  # Explicitly remove null characters which can cause parsing issues
  json_str <- gsub("\\0", "", json_str, fixed = TRUE)

  jsonlite::fromJSON(json_str, simplifyVector = TRUE,
                     simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
}

#' Retrieve a quantization report
#'
#' Thin wrapper around \code{lna_get_transform_report()} for convenience.
#'
#' @inheritParams lna_get_transform_report
#' @export
lna_get_quant_report <- function(lna_file, transform_index_or_name) {
  lna_get_transform_report(lna_file, transform_index_or_name)
}
