#' Initiate an LNA pipeline
#'
#' Creates a new `lna_pipeline` object and sets its input using
#' `lna_pipeline$set_input()`.
#'
#' @param x Data object or list of run data.
#' @param run_ids Optional character vector of run identifiers.
#' @param chunk_mb_suggestion Optional numeric hint for chunk size.
#'
#' @return A configured `lna_pipeline` object.
#' @export
as_pipeline <- function(x, run_ids = NULL, chunk_mb_suggestion = NULL) {
  pipe <- lna_pipeline$new()
  pipe$set_input(x, run_ids = run_ids, chunk_mb_suggestion = chunk_mb_suggestion)
  pipe
}
