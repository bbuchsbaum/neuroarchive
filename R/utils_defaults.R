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

#' Resolve Transform Parameters
#'
#' Merges transform parameters from schema defaults, package options, and
#' user supplied values (in that order). Performs a deep merge using
#' `utils::modifyList` with left-to-right precedence.
#'
#' @param transforms Character vector of transform types.
#' @param transform_params Named list of user-supplied parameters.
#' @return Named list of merged parameter lists.
#' @keywords internal
resolve_transform_params <- function(transforms, transform_params = list()) {
  stopifnot(is.character(transforms))
  stopifnot(is.list(transform_params))

  if (length(transform_params) > 0) {
    if (is.null(names(transform_params)) || any(names(transform_params) == "")) {
      abort_lna(
        "transform_params must be a named list",
        .subclass = "lna_error_validation"
      )
    }

    unknown <- setdiff(names(transform_params), transforms)
    if (length(unknown) > 0) {
      abort_lna(
        paste0(
          "Unknown transform(s) in transform_params: ",
          paste(unknown, collapse = ", ")
        ),
        .subclass = "lna_error_validation"
      )
    }
  }

  pkg_opts <- lna_options()
  merged <- setNames(vector("list", length(transforms)), transforms)

  for (type in transforms) {
    defaults <- default_params(type)
    pkg_default <- pkg_opts[[type]]
    user <- transform_params[[type]]

    params <- defaults
    if (is.list(pkg_default)) {
      params <- utils::modifyList(params, pkg_default, keep.null = TRUE)
    }
    if (is.list(user)) {
      params <- utils::modifyList(params, user, keep.null = TRUE)
    }
    merged[[type]] <- params
  }

  merged
}
