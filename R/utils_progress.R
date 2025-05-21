#' Check if progressr handlers are active and not void
#'
#' This function checks if progressr has any active handlers and
#' ensures that not all of them are "void" handlers (which silence output).
#'
#' @return Logical, TRUE if progress reporting is effectively enabled, FALSE otherwise.
#' @keywords internal
is_progress_globally_enabled <- function() {
  active_handlers_list <- progressr::handlers()
  if (length(active_handlers_list) == 0) {
    return(FALSE)
  }
  !all(sapply(active_handlers_list, function(h) inherits(h, "handler_void")))
} 