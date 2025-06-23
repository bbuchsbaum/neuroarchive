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
  transforms_root <- lna_options("paths.transforms_root")[[1]]
  report_path <- params$report_path %||% paste0(transforms_root, base, "_report.json")

  root <- h5[["/"]]
  assert_h5_path(root, report_path)

  dset <- root[[report_path]]
  on.exit(safe_h5_close(dset), add = TRUE)
  
  # When the HDF5 datatype is H5T_STD_U8LE (unsigned 8-bit int),
  # dset$read() returns a numeric vector of byte values.
  report_bytes_numeric <- dset$read()

  # Convert this numeric vector of bytes into an R raw vector.
  current_raw_vector <- as.raw(report_bytes_numeric)
  
  # Handle decompression if the attribute exists.
  if (h5_attr_exists(dset, "compression")) {
    comp <- h5_attr_read(dset, "compression")
    if (identical(comp, "gzip")) {
      current_raw_vector <- memDecompress(current_raw_vector, type = "gzip")
    }
  } else {
    # Check for gzip magic number (0x1f8b or 0x789c for different gzip formats)
    # 0x78 0x9c is the most common gzip header for default compression level
    if (length(current_raw_vector) >= 2) {
      if ((current_raw_vector[1] == as.raw(0x78) && current_raw_vector[2] == as.raw(0x9c)) ||
          (current_raw_vector[1] == as.raw(0x1f) && current_raw_vector[2] == as.raw(0x8b))) {
        current_raw_vector <- memDecompress(current_raw_vector, type = "gzip")
      }
    }
  }

  # Find the first null byte, if any, and truncate the raw vector there.
  # This prevents rawToChar from erroring on embedded nulls if the actual content ends before them.
  first_null_idx <- which(current_raw_vector == as.raw(0))
  if (length(first_null_idx) > 0) {
    current_raw_vector <- current_raw_vector[1:(min(first_null_idx) - 1)]
  }

  # Convert the (potentially decompressed and truncated) raw vector to a character string.
  json_str <- rawToChar(current_raw_vector)
  
  # Explicitly mark encoding as UTF-8 before iconv, in case rawToChar doesn't or marks it differently.
  Encoding(json_str) <- "UTF-8"
  
  # Clean the string:
  # 1. Sanitize UTF-8 encoding first, as subsequent operations assume valid UTF-8.
  #    Remove invalid sequences.
  json_str <- iconv(json_str, from = "UTF-8", to = "UTF-8", sub = "") 
  
  # 2. Trim whitespace.
  json_str <- trimws(json_str)

  # Parse the cleaned JSON string.
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
