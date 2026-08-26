# Fresh-process numerical court for an installed BIDS exchange restoration.

suppressPackageStartupMessages({
  library(neuroarchive)
  library(neuroim2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("expected values RDS, restored root, and receipt paths", call. = FALSE)
}
expected <- readRDS(args[[1L]])
restored_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
records <- lapply(names(expected), function(record_id) {
  item <- expected[[record_id]]
  path <- file.path(restored_root, item$restored_relative)
  if (!file.exists(path)) stop("missing restored NIfTI: ", path, call. = FALSE)
  restored <- neuroim2::read_vec(path)
  actual <- as.vector(as.array(restored))
  wanted <- as.vector(item$values)
  exact <- identical(
    writeBin(as.double(actual), raw(), size = 8L, endian = "little"),
    writeBin(as.double(wanted), raw(), size = 8L, endian = "little")
  )
  if (!isTRUE(exact)) {
    stop("lossless numerical mismatch for ", record_id, call. = FALSE)
  }
  list(
    record_id = record_id,
    status = "pass",
    dims = as.list(dim(restored)),
    n_values = length(actual),
    max_abs_error = max(abs(actual - wanted)),
    rmse = sqrt(mean((actual - wanted)^2))
  )
})
receipt <- list(
  schema_version = "1.0.0",
  status = "pass",
  package_version = as.character(utils::packageVersion("neuroarchive")),
  records = records
)
jsonlite::write_json(
  receipt, args[[3L]], auto_unbox = TRUE, null = "null", digits = NA,
  pretty = TRUE
)
message("BIDS exchange fresh-process numerical court: PASS")
