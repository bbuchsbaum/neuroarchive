#' Quant dataset paths helper
#'
#' Returns standard dataset paths for quantized data given a run ID.
#'
#' @param run_id Character run identifier.
#' @return Named list with `data`, `scale`, and `offset` paths.
#' @keywords internal
quant_paths <- function(run_id) {
  run_id <- sanitize_run_id(run_id)
  list(
    data = paste0("/scans/", run_id, "/quantized"),
    scale = paste0("/scans/", run_id, "/quant_scale"),
    offset = paste0("/scans/", run_id, "/quant_offset")
  )
}

