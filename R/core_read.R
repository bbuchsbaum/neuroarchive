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
#' @param output_dtype Desired output data type. One of
#'   `"float32"`, `"float64"`, or `"float16"`.
#' @param lazy Logical. If `TRUE`, the HDF5 file handle remains open
#'   after return (for lazy reading).
#'
#' @return A `DataHandle` object representing the loaded data.
#' @import hdf5r
#' @keywords internal
core_read <- function(file, allow_plugins = c("warn", "off", "on"), validate = FALSE,
                      output_dtype = c("float32", "float64", "float16"),
                      lazy = FALSE) {
  allow_plugins <- match.arg(allow_plugins)
  output_dtype <- match.arg(output_dtype)
  h5 <- open_h5(file, mode = "r")
  if (!lazy) {
    on.exit(close_h5_safely(h5))
  }

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

  if (identical(output_dtype, "float16") && !has_float16_support()) {
    abort_lna(
      "float16 output not supported",
      .subclass = "lna_error_float16_unsupported",
      location = sprintf("core_read:%s", file)
    )
  }
  handle$meta$output_dtype <- output_dtype

  return(handle)
}
