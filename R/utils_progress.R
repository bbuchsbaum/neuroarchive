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

#' Run a loop with optional progress reporting
#'
#' This helper creates a progressr progressor when progress reporting is
#' enabled and executes the provided loop function with that progressor.
#'
#' @param steps Integer number of steps to report progress for.
#' @param loop_fn Function accepting a single argument (the progressor or
#'   `NULL`) which performs the work and returns a value.
#' @return The value returned by `loop_fn`.
#' @keywords internal
with_progress_loop <- function(steps, loop_fn) {
  progress_enabled <- steps > 1 && is_progress_globally_enabled()

  run_loop <- function() {
    p <- if (progress_enabled) progressr::progressor(steps = steps) else NULL
    loop_fn(p)
  }

  if (progress_enabled) {
    progressr::with_progress(run_loop())
  } else {
    run_loop()
  }
}
