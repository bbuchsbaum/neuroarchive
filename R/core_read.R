#' Core LNA Read Routine
#'
#' @description Opens an LNA HDF5 file, discovers available transform
#'   descriptors and runs the inverse transform chain.
#'
#' @param file Path to an LNA file on disk.
#' @param run_id Character vector of run identifiers or glob patterns. If
#'   `NULL`, the first available run is used. Glob patterns are matched
#'   against available run groups under `/scans`.
#' @param allow_plugins Character. How to handle transforms requiring
#'   external packages. One of "installed" (default), "none", or "prompt".
#' @param validate Logical flag indicating if validation should be
#'   performed via `validate_lna()` before reading.
#' @param output_dtype Desired output data type. One of
#'   `"float32"`, `"float64"`, or `"float16"`.
#' @param roi_mask Optional logical mask to subset voxels.
#' @param time_idx Optional integer vector selecting time points.
#' @param lazy Logical. If `TRUE`, the HDF5 file handle remains open
#'   after return (for lazy reading).
#'
#' @return If a single run is selected, a `DataHandle` object. When
#'   multiple runs match and `lazy = FALSE`, a named list of `DataHandle`
#'   objects is returned.
#' @import hdf5r
#' @keywords internal
core_read <- function(file, run_id = NULL,
                      allow_plugins = c("installed", "none", "prompt"),
                      validate = FALSE,
                      output_dtype = c("float32", "float64", "float16"),
                      roi_mask = NULL, time_idx = NULL,
                      lazy = FALSE) {
  allow_plugins <- normalize_allow_plugins(allow_plugins)
  output_dtype <- match.arg(output_dtype)

  h5 <- open_h5(file, mode = "r")
  if (!lazy) {
    on.exit(neuroarchive:::close_h5_safely(h5))
  }

  runs <- resolve_runs_for_read(h5, run_id, lazy)

  if (validate) {
    validate_lna(file)
  }

  subset_params <- collect_subset_params(roi_mask, time_idx)
  tf_group <- h5[["transforms"]]
  transforms <- prepare_transforms_for_read(tf_group, allow_plugins, file)

  results <- lapply(runs, function(rid) {
    process_run_core_read(
      rid = rid,
      h5 = h5,
      runs = runs,
      subset_params = subset_params,
      transforms = transforms,
      tf_group = tf_group,
      validate = validate,
      output_dtype = output_dtype,
      allow_plugins = allow_plugins,
      file = file
    )
  })
  names(results) <- runs

  if (length(results) == 1) results[[1]] else results
}
