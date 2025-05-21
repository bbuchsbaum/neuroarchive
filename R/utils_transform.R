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

#' Handle missing transform implementations
#'
#' Internal helper used by `core_read` and `lna_reader` to process cases
#' where a transform's S3 methods are unavailable. Behaviour depends on
#' the `allow_plugins` mode.
#'
#' @param missing_types Character vector of transform types lacking
#'   implementations.
#' @param allow_plugins One of "installed", "none", or "prompt".
#' @param location Optional string used in error conditions.
#' @keywords internal
handle_missing_methods <- function(missing_types, allow_plugins, location = NULL) {
  if (length(missing_types) == 0) return(character())

  msg <- paste0(
    "Missing invert_step implementation for transform(s): ",
    paste(unique(missing_types), collapse = ", ")
  )

  if (identical(allow_plugins, "none")) {
    abort_lna(msg, .subclass = "lna_error_no_method", location = location)
  } else if (identical(allow_plugins, "prompt") && rlang::is_interactive()) {
    response <- tolower(trimws(readline(paste0(msg, " Continue anyway? [y/N]: "))))
    if (!response %in% c("y", "yes")) {
      abort_lna(msg, .subclass = "lna_error_no_method", location = location)
    }
    warning(msg, call. = FALSE)
  } else {
    warning(msg, call. = FALSE)
  }

  invisible(missing_types)
}
