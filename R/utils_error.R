
#' LNA Error Handling Helpers
#'
#' Provides a thin wrapper around `rlang::abort` for package specific
#' error classes used throughout the code base.
#'
#' @param message A character string describing the error.
#' @param ... Additional named data stored in the condition object.
#' @param .subclass Character string giving the LNA error subclass.
#' @return No return value. This function always throws an error.
#' @keywords internal
abort_lna <- function(message, ..., .subclass, location = NULL, parent = NULL) {
  stopifnot(is.character(message), length(message) == 1)
  stopifnot(is.character(.subclass))
  rlang::abort(
    message,
    ...,
    .subclass = .subclass,
    location = location,
    parent = parent
  )
}

#' LNA Warning Helper
#'
#' Provides a thin wrapper around `warning` for package specific
#' warnings. Mainly used for integration checks where execution should
#' continue but the user ought to be informed.
#'
#' @param message A character string describing the warning.
#' @param ... Additional named data stored in the condition object.
#' @param .subclass Character string giving the LNA warning subclass.
#' @return No return value. This function is called for its side effect of
#'   signalling a warning.
#' @keywords internal
warn_lna <- function(message, ..., .subclass = "lna_warning", location = NULL) {
  stopifnot(is.character(message), length(message) == 1)
  warning(message, call. = FALSE)
}

#' Error thrown when `lna_reader` methods are called after the reader is closed.
#'
#' @keywords internal
lna_error_closed_reader <- NULL
