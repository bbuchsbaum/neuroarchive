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
      .subclass = "lna_error_validation",
      location = "apply_template:template_name"
    )
  }
  reg <- get(".template_registry", envir = lna_template_registry_env)
  if (!exists(template_name, envir = reg, inherits = FALSE)) {
    abort_lna(
      sprintf("Template '%s' not registered", template_name),
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
