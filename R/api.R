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
#' @param header Optional named list of header attributes.
#' @return Invisibly returns a list with elements `file`, `plan`, and
#'   `header` with class `"lna_write_result"`.
#' @export
write_lna <- function(x, file = NULL, transforms = character(),
                      transform_params = list(), mask = NULL,
                      header = NULL) {
  in_memory <- FALSE
  if (is.null(file)) {
    tmp <- tempfile(fileext = ".h5")
    file <- tmp
    in_memory <- TRUE
  }

  result <- core_write(x = x, transforms = transforms,
                       transform_params = transform_params,
                       mask = mask, header = header)

  if (in_memory) {
    h5 <- open_h5(file, mode = "w", driver = "core", backing_store = FALSE)
  } else {
    h5 <- open_h5(file, mode = "w")
  }
  materialise_plan(h5, result$plan, header = result$handle$meta$header)
  close_h5_safely(h5)

  out_file <- if (in_memory) NULL else file
  out <- list(file = out_file, plan = result$plan,
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
#' @param output_dtype Desired output data type. One of
#'   `"float32"`, `"float64"`, or `"float16"`.
#' @param lazy Logical. If `TRUE`, the HDF5 file remains open and the
#'   returned handle can be used for lazy reading (Phase 1 stub).
#' @return A `DataHandle` object from `core_read`.
#' @export
read_lna <- function(file, allow_plugins = c("warn", "off", "on"),
                     validate = FALSE,
                     output_dtype = c("float32", "float64", "float16"),
                     lazy = FALSE) {
  core_read(file = file, allow_plugins = allow_plugins,
            validate = validate, output_dtype = output_dtype,
            lazy = lazy)
}
