#' DSL Template Registry
#'
#' Provides functions to register and apply pipeline templates.
#'
#' @param template_name A single string naming the template.
#' @param template_function A function taking `(pipeline_obj, ...)` and
#'   returning a modified `lna_pipeline`.
#' @param force Overwrite an existing registration.
#' @return Invisibly returns a list with the registered name.
#' @export
lna_template_registry_env <- new.env(parent = emptyenv())
assign(".template_registry", new.env(parent = emptyenv()), envir = lna_template_registry_env)




# internal environment holding the template registry
lna_template_registry_env <- new.env(parent = emptyenv())
assign(".template_registry", new.env(parent = emptyenv()), envir = lna_template_registry_env)

#' Register a pipeline template
#'
#' @param template_name Single string naming the template.
#' @param template_function Function implementing the template. It must
#'   accept a `lna_pipeline` as its first argument.
#' @param force Overwrite existing template with the same name.
#'
#' @return Invisibly returns a list with the registered `name`.
#' @export
register_lna_template <- function(template_name, template_function, force = FALSE) {
  if (missing(template_name) || !is.character(template_name) || length(template_name) != 1) {
    abort_lna(
      "template_name must be a single character string",
      .subclass = "lna_error_validation",
      location = "register_lna_template:template_name"
    )
  }

  if (!is.function(template_function)) {
    abort_lna(
      "template_function must be a function",
      .subclass = "lna_error_validation",
      location = "register_lna_template:template_function"
    )
  }

  reg <- get(".template_registry", envir = lna_template_registry_env)
  if (!force && exists(template_name, envir = reg, inherits = FALSE)) {
    warning(
      sprintf("Template '%s' already registered; use force=TRUE to replace", template_name),
      call. = FALSE
    )
    return(invisible(list(name = template_name)))
  }

  assign(template_name, template_function, envir = reg)
  invisible(list(name = template_name))
}

#' Apply a registered pipeline template
#'
#' @param pipeline_obj An `lna_pipeline` object.
#' @param template_name Name of a registered template.
#' @param ... Passed to the template function and also interpreted as
#'   parameter overrides (e.g. `pca.k = 120`).
#' @return The modified `lna_pipeline` object.
#' @export
apply_template <- function(pipeline_obj, template_name, ...) {
  if (!inherits(pipeline_obj, "lna_pipeline")) {
    abort_lna(
      "pipeline_obj must be an lna_pipeline",
      .subclass = "lna_error_validation",
      location = "apply_template:pipeline_obj"
    )
  }
  if (!is.character(template_name) || length(template_name) != 1) {
    abort_lna(
      "template_name must be a single string",

#' @param data_or_pipe Data object or `lna_pipeline`.
#' @param template_name Name of the template to apply.
#' @param ... Additional arguments passed to the template function and
#'   interpreted as parameter overrides.
#'
#' @return The modified `lna_pipeline` object.
#' @export
apply_template <- function(data_or_pipe, template_name, ...) {
  pipe <- if (inherits(data_or_pipe, "lna_pipeline")) {
    data_or_pipe
  } else {
    as_pipeline(data_or_pipe)
  }

  if (missing(template_name) || !is.character(template_name) || length(template_name) != 1) {
    abort_lna(
      "template_name must be a single character string",

      .subclass = "lna_error_validation",
      location = "apply_template:template_name"
    )
  }
  reg <- get(".template_registry", envir = lna_template_registry_env)
  if (!exists(template_name, envir = reg, inherits = FALSE)) {
    abort_lna(
      sprintf("Template '%s' not registered", template_name),


  reg <- get(".template_registry", envir = lna_template_registry_env)
  if (!exists(template_name, envir = reg, inherits = FALSE)) {
    abort_lna(
      sprintf("Template '%s' not found", template_name),

      .subclass = "lna_error_validation",
      location = "apply_template:template_name"
    )
  }

  fun <- get(template_name, envir = reg)
  args <- list(...)
  pipe <- fun(pipeline_obj, ...)
  if (length(args)) {
    for (nm in names(args)) {
      val <- args[[nm]]
      if (is.null(nm) || nm == "") next
      if (grepl("\\.", nm, fixed = TRUE)) {
        parts <- strsplit(nm, ".", fixed = TRUE)[[1]]
        if (length(parts) >= 2) {
          type <- parts[1]
          param <- parts[2]
          pipe$modify_step(type, setNames(list(val), param))
        }
      } else if (is.list(val)) {
        pipe$modify_step(nm, val)
      }
    }
  }
  pipe
}


  fun <- get(template_name, envir = reg, inherits = FALSE)
  dots <- list(...)
  pipe <- do.call(fun, c(list(pipe), dots))

  if (length(dots) > 0) {
    overrides <- list()

    # dotted name style transform.param=value
    named <- names(dots)
    if (!is.null(named)) {
      for (nm in named) {
        if (grepl("\\.", nm, fixed = TRUE)) {
          parts <- strsplit(nm, "\\.", fixed = TRUE)[[1]]
          if (length(parts) >= 2) {
            step_type <- paste(parts[-length(parts)], collapse = ".")
            param <- parts[length(parts)]
            overrides[[step_type]][[param]] <- dots[[nm]]
          }
        } else if (is.list(dots[[nm]])) {
          overrides[[nm]] <- utils::modifyList(overrides[[nm]] %||% list(), dots[[nm]])
        }
      }
    }

    if (length(overrides) > 0) {
      for (st in names(overrides)) {
        pipe$modify_step(st, overrides[[st]])
      }
    }
  }

  pipe
}

