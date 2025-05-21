
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

#' Error thrown when `lna_reader` methods are called after the reader is closed.
#'
#' @keywords internal
lna_error_closed_reader <- NULL
