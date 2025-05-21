#' Default Parameter Retrieval
#'
#' Loads parameter defaults from a JSON schema for the given transform type and
#' caches both the schema and the resulting defaults. If a schema cannot be
#' found, an empty list is cached and returned with a warning.
#'
#' @param type Character transform type.
#' @return A list of default parameters extracted from the schema, or an empty
#'   list if none are defined or the schema is missing.
#' @keywords internal
.default_param_cache <- new.env(parent = emptyenv())


# Recursively extract `default` values from a parsed JSON schema list
.extract_schema_defaults <- function(node) {
  if (!is.list(node)) {
    return(NULL)
  }

  if (!is.null(node$default)) {
    return(node$default)
  }

  defaults <- list()
  if (is.list(node$properties)) {
    for (nm in names(node$properties)) {
      val <- .extract_schema_defaults(node$properties[[nm]])
      if (!is.null(val)) {
        defaults[[nm]] <- val
      }
    }
  }

  if (length(defaults) > 0) defaults else NULL
}

#' Null-coalescing helper
#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b


default_params <- function(type) {
  stopifnot(is.character(type), length(type) == 1)

  cache <- .default_param_cache
  if (exists(type, envir = cache, inherits = FALSE)) {
    return(cache[[type]])
  }

  schema_file <- system.file("schemas", paste0(type, ".schema.json"),
                             package = "neuroarchive")
  defaults <- list()
  if (nzchar(schema_file) && file.exists(schema_file)) {
    schema <- jsonlite::read_json(schema_file, simplifyVector = TRUE)
    if (is.list(schema$properties)) {
      for (nm in names(schema$properties)) {
        prop <- schema$properties[[nm]]
        if (!is.null(prop$default)) {
          defaults[[nm]] <- prop$default
        }
      }
    }
  }

  schema_path <- system.file("schemas", paste0(type, ".schema.json"),
                              package = "neuroarchive")

  if (schema_path == "") {
    warning(sprintf("Schema for transform '%s' not found", type), call. = FALSE)
    defaults <- list()
  } else {
    schema <- jsonlite::read_json(schema_path, simplifyVector = FALSE)
    assign(type, schema, envir = .schema_cache)
    defaults <- .extract_schema_defaults(schema)
  }

  assign(type, defaults, envir = .default_param_cache)

  cache[[type]] <- defaults
  defaults
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

#' Default parameters for the 'quant' transform
#'
#' Convenience wrapper around `default_params("quant")`.
#' @export
lna_default.quant <- function() {
  default_params("quant")
}
