#' DSL Verb Registry
#'
#' Provides a registry for DSL verbs so that external packages can
#' expose additional pipeline verbs.  The registry simply maps a verb
#' function name to the corresponding LNA transform type.
#'
#' @format An environment containing the verb registry.
#' @export
lna_verb_registry_env <- new.env(parent = emptyenv())
assign(".verb_registry", new.env(parent = emptyenv()), envir = lna_verb_registry_env)

#' Register LNA Verb
#'
#' Register a DSL verb that maps to an LNA transform type.
#'
#' @param verb_name A symbol or single string naming the verb.  If
#'   `NULL` and `default_slug` is `TRUE`, the name is derived from
#'   `lna_transform_type` by replacing non-alphanumeric characters with
#'   underscores.
#' @param lna_transform_type Character string identifying the LNA
#'   transform.
#' @param default_slug Logical flag controlling whether a missing
#'   `verb_name` is slugged from `lna_transform_type`.
#' @param force Overwrite an existing registration with the same name.
#'
#' @return Invisibly returns a list with the registered `name` and
#'   `type`.
#' @keywords internal
register_lna_verb <- function(verb_name = NULL, lna_transform_type,
                              default_slug = TRUE, force = FALSE) {
  if (missing(lna_transform_type) ||
      !is.character(lna_transform_type) || length(lna_transform_type) != 1) {
    abort_lna(
      "lna_transform_type must be a single character string",
      .subclass = "lna_error_validation",
      location = "register_lna_verb:lna_transform_type"
    )
  }

  if (rlang::is_symbol(verb_name)) {
    verb_name <- rlang::as_string(verb_name)
  }

  if (is.null(verb_name) || !nzchar(verb_name)) {
    if (!default_slug) {
      abort_lna(
        "verb_name must be provided when default_slug = FALSE",
        .subclass = "lna_error_validation",
        location = "register_lna_verb:verb_name"
      )
    }
    verb_name <- gsub("[^A-Za-z0-9]+", "_", lna_transform_type)
  }

  if (!is.character(verb_name) || length(verb_name) != 1) {
    abort_lna(
      "verb_name must be a single string or symbol",
      .subclass = "lna_error_validation",
      location = "register_lna_verb:verb_name"
    )
  }

  reg <- get(".verb_registry", envir = lna_verb_registry_env)
  if (!force && exists(verb_name, envir = reg, inherits = FALSE)) {
    warning(
      sprintf("Verb '%s' already registered; use force=TRUE to replace", verb_name),
      call. = FALSE
    )
    return(invisible(list(name = verb_name, type = lna_transform_type)))
  }

  assign(verb_name, lna_transform_type, envir = reg)
  invisible(list(name = verb_name, type = lna_transform_type))
}
