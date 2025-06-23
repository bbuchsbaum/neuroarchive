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
.required_param_cache <- new.env(parent = emptyenv())

#' Clear the default parameter cache
#'
#' Removes all cached default parameter lists.
#'
#' @return invisible(NULL)
#' @keywords internal
default_param_cache_clear <- function() {
  rm(list = ls(envir = .default_param_cache, all.names = TRUE), envir = .default_param_cache)
  memoise::forget(default_params)
  invisible(NULL)
}

#' Clear the required parameter cache
#'
#' Removes all cached required parameter vectors.
#'
#' @return invisible(NULL)
#' @keywords internal
required_param_cache_clear <- function() {
  rm(list = ls(envir = .required_param_cache, all.names = TRUE), envir = .required_param_cache)
  invisible(NULL)
}


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

  # Handle oneOf patterns where actual properties are in oneOf[0]
  if (is.list(node$oneOf) && length(node$oneOf) > 0) {
    first_option <- node$oneOf[[1]]
    if (is.list(first_option$properties)) {
      for (nm in names(first_option$properties)) {
        val <- .extract_schema_defaults(first_option$properties[[nm]])
        if (!is.null(val)) {
          defaults[[nm]] <- val
        }
      }
    }
  }

  if (is.list(node$items)) {
    if (is.null(names(node$items))) {
      item_vals <- lapply(node$items, .extract_schema_defaults)
      if (any(vapply(item_vals, Negate(is.null), logical(1)))) {
        defaults$items <- item_vals
      }
    } else {
      val <- .extract_schema_defaults(node$items)
      if (!is.null(val)) {
        defaults$items <- val
      }
    }
  }

  if (length(defaults) > 0) defaults else NULL
}

# Null-coalescing helper (internal, undocumented)
`%||%` <- function(a, b) if (!is.null(a)) a else b


default_params_impl <- function(type) {
  stopifnot(is.character(type), length(type) == 1)

  cache <- .default_param_cache
  if (exists(type, envir = cache, inherits = FALSE)) {
    return(cache[[type]])
  }

  pkgs <- unique(c("neuroarchive", loadedNamespaces()))
  schema_path <- ""
  for (pkg in pkgs) {
    path <- system.file("schemas", paste0(type, ".schema.json"), package = pkg)
    if (nzchar(path) && file.exists(path)) {
      schema_path <- path
      break
    }
  }

  if (!nzchar(schema_path)) {
    warning(sprintf("Schema for transform '%s' not found", type), call. = FALSE)
    defaults <- list()
  } else {
    schema <- jsonlite::read_json(schema_path, simplifyVector = FALSE)
    assign(type, schema, envir = .schema_cache)
    
    defaults <- .extract_schema_defaults(schema) %||% list()
    
    # Special handling for params.oneOf pattern (e.g., delta.schema.json)
    # If the extraction found defaults nested under 'params', flatten them to top level
    if (length(defaults) == 1 && !is.null(defaults$params) && is.list(defaults$params)) {
      defaults <- defaults$params
    }
  }

  assign(type, defaults, envir = cache)
  defaults
}

# Memoised wrapper -----------------------------------------------------------

default_params <- memoise::memoise(default_params_impl)

#' Required parameters for a transform
#'
#' Retrieves the `required` fields from a transform's JSON schema. Results
#' are cached for efficiency.
#'
#' @param type Character scalar transform type.
#' @return Character vector of required parameter names (may be empty).
#' @keywords internal
required_params <- function(type) {
  stopifnot(is.character(type), length(type) == 1)

  cache <- .required_param_cache
  if (exists(type, envir = cache, inherits = FALSE)) {
    return(cache[[type]])
  }

  pkgs <- unique(c("neuroarchive", loadedNamespaces()))
  schema_path <- ""
  for (pkg in pkgs) {
    path <- system.file("schemas", paste0(type, ".schema.json"), package = pkg)
    if (nzchar(path) && file.exists(path)) {
      schema_path <- path
      break
    }
  }

  if (!nzchar(schema_path)) {
    warning(sprintf("Schema for transform '%s' not found", type), call. = FALSE)
    req <- character()
  } else {
    schema <- jsonlite::read_json(schema_path, simplifyVector = FALSE)
    assign(type, schema, envir = .schema_cache)
    req <- schema$required %||% character()
  }

  assign(type, req, envir = cache)
  req
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
        .subclass = "lna_error_validation",
        location = "resolve_transform_params"
      )
    }

    unknown <- setdiff(names(transform_params), transforms)
    if (length(unknown) > 0) {
      abort_lna(
        paste0(
          "Unknown transform(s) in transform_params: ",
          paste(unknown, collapse = ", ")
        ),
        .subclass = "lna_error_validation",
        location = "resolve_transform_params"
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

#' @title Default parameters for the 'quant' transform
#' @description Convenience wrapper around `default_params("quant")`.
#' @seealso default_params
#' @export
lna_default.quant <- function() {
  default_params("quant")
}

#' @title Default parameters for the 'basis' transform
#' @description Convenience wrapper around `default_params("basis")`.
#' @seealso default_params
#' @export
lna_default.basis <- function() {
  default_params("basis")
}

#' @title Default parameters for the 'embed' transform
#' @description Convenience wrapper around `default_params("embed")`.
#' @seealso default_params
#' @export
lna_default.embed <- function() {
  default_params("embed")
}


#' @title Default parameters for the 'delta' transform
#' @description Convenience wrapper around `default_params("delta")`.
#' @seealso default_params
#' @export
lna_default.delta <- function() {
  default_params("delta")
}

#' @title Default parameters for the 'temporal' transform
#' @description Convenience wrapper around `default_params("temporal")`.
#' @seealso default_params
#' @export
lna_default.temporal <- function() {
  default_params("temporal")
}
