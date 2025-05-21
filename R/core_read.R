#' Core LNA Read Routine
#'
#' @description Opens an LNA HDF5 file, discovers available transform
#'   descriptors and runs the inverse transform chain. This is a
#'   minimal skeleton used during early development.
#'
#' @param file Path to an LNA file on disk.
#' @param allow_plugins Placeholder argument for plugin policy.
#' @param validate Logical flag indicating if validation should be
#'   performed. Currently ignored.
#'
#' @return A `DataHandle` object representing the loaded data.
#' @import hdf5r
#' @keywords internal
core_read <- function(file, allow_plugins = c("warn", "off", "on"), validate = FALSE) {
  allow_plugins <- match.arg(allow_plugins)
  h5 <- H5File$new(file, mode = "r")
  on.exit({
    if (inherits(h5, "H5File") && h5$is_valid()) h5$close_all()
  })

  handle <- DataHandle$new(h5 = h5)
  tf_group <- h5[["transforms"]]

  transforms <- discover_transforms(tf_group)

  # TODO: make use of allow_plugins and validate

  if (nrow(transforms) > 0) {
    for (i in rev(seq_len(nrow(transforms)))) {
      name <- transforms$name[[i]]
      type <- transforms$type[[i]]
      desc <- read_json_descriptor(tf_group, name)
      handle <- invert_step(type, desc, handle)
    }
  }

  return(handle)
}
