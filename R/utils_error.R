#' Error Handling Utilities
#'
#' @description
#' Provides helper wrappers around `rlang::abort` for LNA specific
#' error classes.
#' @keywords internal
#' @import rlang
abort_lna <- function(message, ..., .subclass = NULL) {
  rlang::abort(message = message, ..., .subclass = .subclass)
}
