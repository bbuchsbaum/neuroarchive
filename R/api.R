#' Write data to an LNA file
#'
#' Compresses and stores a neuroimaging object using the specified
#' transform pipeline.  Parameter values are resolved by merging
#' transform schema defaults, package options set via
#' `lna_options()`, and the user supplied `transform_params`
#' (later values override earlier ones).
#'
#' @param x Input object passed to `core_write`.
#' @param file Path to output `.h5` file. If `NULL`, writing is performed
#'   in memory and the file is returned as a raw vector.
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
#' Loads data from an `.lna.h5` file using `core_read`.  When
#' `lazy = TRUE` the function returns an `lna_reader` object that keeps
#' the HDF5 handle open for on-demand reconstruction of the data.
#'
#' @param file Path to an LNA file on disk.
#' @param allow_plugins Forwarded to `core_read`.
#' @param validate Logical flag for validation; forwarded to `core_read`.
#' @param output_dtype Desired output data type. One of
#'   `"float32"`, `"float64"`, or `"float16"`.
#' @param lazy Logical. If `TRUE`, the HDF5 file remains open and the
#'   returned `lna_reader` can load data lazily.
#' @return A `DataHandle` object from `core_read`.
#' @export
read_lna <- function(file, allow_plugins = c("warn", "off", "on"),
                     validate = FALSE,
                     output_dtype = c("float32", "float64", "float16"),
                     lazy = FALSE) {
  output_dtype <- match.arg(output_dtype)
  allow_plugins <- match.arg(allow_plugins)

  if (lazy) {
    lna_reader$new(
      file = file,
      core_read_args = list(
        allow_plugins = allow_plugins,
        validate = validate,
        output_dtype = output_dtype
      )
    )
  } else {
    core_read(
      file = file,
      allow_plugins = allow_plugins,
      validate = validate,
      output_dtype = output_dtype,
      lazy = FALSE
    )
  }
}
