#' Default Parameter Retrieval (Stub)
#'
#' Provides a cached stub returning an empty list for any transform type.
#'
#' @param type Character transform type.
#' @return A list of default parameters (currently empty).
#' @keywords internal
.default_param_cache <- new.env(parent = emptyenv())

default_params <- function(type) {
  stopifnot(is.character(type), length(type) == 1)
  if (exists(type, envir = .default_param_cache, inherits = FALSE)) {
    return(.default_param_cache[[type]])
  }
  .default_param_cache[[type]] <- list()
  .default_param_cache[[type]]
}
