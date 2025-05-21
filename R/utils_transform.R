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

#' Execute a transform step with provenance
#'
#' Helper used by `core_write` and `core_read` to ensure that errors
#' include the step index and transform name in their location field.
#'
#' @param direction Either "forward" or "invert".
#' @param type Transform type.
#' @param desc Descriptor list.
#' @param handle DataHandle object.
#' @param step_index Integer step index.
#' @keywords internal
run_transform_step <- function(direction, type, desc, handle, step_index) {
  stopifnot(direction %in% c("forward", "invert"))
  fun <- if (direction == "forward") forward_step else invert_step
  tryCatch(
    fun(structure(type, class = c(type, "character")), desc, handle),
    error = function(e) {
      prefix <- sprintf("%s_step.%s[%d]", direction, type, step_index)
      loc <- if (!is.null(e$location)) paste0(prefix, ":", e$location) else prefix
      abort_lna(
        conditionMessage(e),
        .subclass = class(e),
        location = loc,
        step_index = step_index,
        transform = type,
        parent = e
      )
    }
  )
}
