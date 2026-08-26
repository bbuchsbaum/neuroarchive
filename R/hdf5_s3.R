# S3 shims for hdf5r ----------------------------------------------------------

#' @keywords internal
`[[<-.H5File` <- function(x, i, value) {
  if (!is.character(i) || length(i) != 1L) {
    abort_lna("H5File dataset path must be a single string",
              .subclass = "lna_error_validation",
              location = "[[<-.H5File")
  }

  path <- i
  if (startsWith(path, "/")) {
    path <- substring(path, 2L)
  }
  if (!nzchar(path)) {
    abort_lna("H5File dataset path cannot be empty",
              .subclass = "lna_error_validation",
              location = "[[<-.H5File")
  }

  # hdf5r's `[[<-` can error when intermediate groups don't exist and the path
  # contains separators. Route through the package writer to create groups.
  h5_write_dataset(x[["/"]], path, value, compression_level = 0)
  x
}

