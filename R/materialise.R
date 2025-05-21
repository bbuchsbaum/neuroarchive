#' Materialise a Plan to an HDF5 file
#'
#' @description Creates the core group structure of an LNA file and writes
#'   transform descriptors according to the provided `Plan` object. Numeric
#'   payloads are not written at this stage.
#'
#' @param h5 An open `H5File` object from **hdf5r**.
#' @param plan A `Plan` R6 object describing the write operations.
#' @return Invisibly returns the modified `plan`.
#' @import hdf5r
#' @keywords internal
materialise_plan <- function(h5, plan) {
  stopifnot(inherits(h5, "H5File"))
  stopifnot(inherits(plan, "Plan"))

  # Ensure core groups exist
  tf_group    <- if (!h5$exists("transforms")) h5$create_group("transforms") else h5[["transforms"]]
  basis_group <- if (!h5$exists("basis"))       h5$create_group("basis")       else h5[["basis"]]
  scans_group <- if (!h5$exists("scans"))       h5$create_group("scans")       else h5[["scans"]]

  # Root attributes
  root <- h5[["/"]]
  h5_attr_write(root, "lna_spec", "LNA R v2.0")
  h5_attr_write(root, "creator", "lna R package v0.0.1")
  h5_attr_write(root, "required_transforms", character(0))

  # Write transform descriptors
  if (length(plan$descriptors) > 0) {
    for (nm in names(plan$descriptors)) {
      desc <- plan$descriptors[[nm]]
      write_json_descriptor(tf_group, nm, desc)
    }
  }

  # Update dataset definitions
  if (nrow(plan$datasets) > 0) {
    plan$datasets$write_mode_effective <- rep("eager", nrow(plan$datasets))
  }

  invisible(plan)
}
