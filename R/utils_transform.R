#' Check transform implementation for namespace collisions
#'
#' Warns if the provided transform type name collides with
#' core LNA transforms or with names of base R packages.
#'
#' @param type Character scalar transform name.
#' @return Logical `TRUE` invisibly. Called for side effects (warnings).
#' @export
check_transform_implementation <- function(type) {
  stopifnot(is.character(type), length(type) == 1)

  core <- c("quant", "basis", "embed", "temporal", "delta")
  base_pkgs <- rownames(installed.packages(priority = "base"))

  msgs <- character()
  if (type %in% core) {
    msgs <- c(msgs, "core LNA transform")
  }
  if (type %in% base_pkgs) {
    msgs <- c(msgs, "base R package")
  }
  if (length(msgs) > 0) {
    warning(sprintf(
      "Transform type '%s' collides with %s namespace",
      type,
      paste(msgs, collapse = " and ")
    ), call. = FALSE)
  }

  invisible(TRUE)
}
