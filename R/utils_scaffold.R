#' Scaffold files for a custom transform
#'
#' Creates skeleton R code, JSON schema, and unit test for a new transform type.
#'
#' @param type Character scalar name of the transform.
#' @return Invisibly returns a list with created file paths.
#' @export
scaffold_transform <- function(type) {
  assert_scalar_character(type, "type")
  if (!nzchar(type)) {
    stop("type must be a non-empty string", call. = FALSE)
  }

  check_transform_implementation(type)
  r_path <- file.path("R", sprintf("transform_%s.R", type))
  schema_path <- file.path("inst", "schemas", sprintf("%s.schema.json", type))
  test_path <- file.path("tests", "testthat", sprintf("test-transform_%s.R", type))

  if (file.exists(r_path) || file.exists(schema_path) || file.exists(test_path)) {
    stop("Transform files already exist for type: ", type, call. = FALSE)
  }

  dir.create(dirname(r_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(schema_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(test_path), recursive = TRUE, showWarnings = FALSE)

  pkg <- utils::packageName(environment(scaffold_transform))
  r_template <- sprintf(
"forward_step.%1$s <- function(type, desc, handle) {
  params <- %2$s:::default_params('%1$s')
  ## TODO: implement forward transform
  handle
}

invert_step.%1$s <- function(type, desc, handle) {
  params <- %2$s:::default_params('%1$s')
  ## TODO: implement inverse transform
  handle
}
",
    type, pkg)

  writeLines(r_template, r_path)

  schema_template <- "{\n  \"$schema\": \"http://json-schema.org/draft-07/schema#\",\n  \"type\": \"object\",\n  \"properties\": {},\n  \"required\": [],\n  \"additionalProperties\": true\n}\n"
  writeLines(schema_template, schema_path)

  test_template <- sprintf("test_that('scaffold %1$s transform', { expect_true(TRUE) })\n", type)
  writeLines(test_template, test_path)

  invisible(list(r_file = r_path, schema = schema_path, test = test_path))
}
