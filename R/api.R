#' Write data to an LNA file
#'
#' A minimal public wrapper around `core_write` used during early
#' development. Creates the basic HDF5 structure via `materialise_plan`.
#'
#' @param x Input object passed to `core_write`.
#' @param file Path to output `.h5` file. If `NULL`, a temporary file is used.
#' @param transforms Character vector of transform types.
#' @param transform_params Named list of transform parameters.
#' @param mask Optional mask passed to `core_write`.
#' @return Invisibly returns a list with elements `file`, `plan`, and
#'   `header` with class `"lna_write_result"`.
#' @export
write_lna <- function(x, file = NULL, transforms = character(),
                      transform_params = list(), mask = NULL) {
  if (is.null(file)) {
    file <- tempfile(fileext = ".h5")
  }

  result <- core_write(x = x, transforms = transforms,
                       transform_params = transform_params,
                       mask = mask)

  h5 <- hdf5r::H5File$new(file, mode = "w")
  materialise_plan(h5, result$plan)
  h5$close_all()

  out <- list(file = file, plan = result$plan,
              header = result$handle$meta$header)
  class(out) <- c("lna_write_result", class(out))
  invisible(out)
}

#' Read data from an LNA file
#'
#' Thin wrapper around `core_read`.
#'
#' @param file Path to an LNA file on disk.
#' @param allow_plugins Forwarded to `core_read`.
#' @param validate Logical flag for validation; forwarded to `core_read`.
#' @return A `DataHandle` object from `core_read`.
#' @export
read_lna <- function(file, allow_plugins = c("warn", "off", "on"),
                     validate = FALSE) {
  core_read(file = file, allow_plugins = allow_plugins,
            validate = validate)
}
