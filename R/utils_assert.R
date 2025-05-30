#' Assertion Helpers
#'
#' Provides internal assertion functions used across the package.
#'
#' @param x Object to test.
#' @param arg_name Name of the argument. Used in error messages.
#' @return Invisibly returns `TRUE` when the assertion passes.
#' @keywords internal
assert_scalar_character <- function(x, arg_name = deparse(substitute(x))) {
  if (!is.character(x) || length(x) != 1 || is.na(x)) {
    stop(sprintf("%s must be a non-missing character scalar", arg_name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Assert object is an H5Group
#'
#' @param obj Object to test.
#' @return Invisibly returns `TRUE` when `obj` is an `H5Group`.
#' @keywords internal
assert_h5group <- function(obj) {
  if (!inherits(obj, "H5Group")) {
    stop("Object must be an H5Group", call. = FALSE)
  }
  invisible(TRUE)
}
