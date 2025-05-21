#' Materialise Plan to HDF5 (Stub)
#'
#' @description Minimal implementation that writes transform descriptors
#'   and core groups to an open HDF5 file. Numeric datasets are not written.
#' @param h5 An open `H5File` object.
#' @param plan A `Plan` R6 object produced by `core_write`.
#' @keywords internal
materialise_plan <- function(h5, plan) {
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

  # Mark datasets as eagerly written
  if (nrow(plan$datasets) > 0) {
    plan$datasets$write_mode_effective[] <- "eager"
  }

  invisible(h5)
}
