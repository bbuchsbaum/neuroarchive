#' Validate an LNA file
#'
#' @description Basic validator that checks the LNA specification version,
#' optional SHA256 checksum and, if available, validates transform
#' descriptors against their JSON schemas.
#'
#' @param file Path to the `.h5` file to validate.
#' @param strict Logical. If `TRUE` (default) validation failures abort with
#'   `lna_error_validation`. If `FALSE`, the function returns `FALSE` on the
#'   first failure and issues a warning.
#' @param checksum Logical. If `TRUE` (default) verify the `lna_checksum`
#'   attribute when present.
#'
#' @return `TRUE` if validation succeeds. If `strict = FALSE` and a problem is
#'   detected the function returns `FALSE`.
#' @export
validate_lna <- function(file, strict = TRUE, checksum = TRUE) {
  stopifnot(is.character(file), length(file) == 1)

  h5 <- open_h5(file, mode = "r")
  on.exit(close_h5_safely(h5))
  root <- h5[["/"]]

  fail <- function(msg) {
    if (strict) {
      abort_lna(msg, .subclass = "lna_error_validation")
    } else {
      warning(msg)
      return(FALSE)
    }
  }

  if (!h5_attr_exists(root, "lna_spec")) {
    res <- fail("Missing lna_spec attribute")
    if (!is.null(res)) return(res)
  } else {
    spec <- h5_attr_read(root, "lna_spec")
    if (!identical(spec, "LNA R v2.0")) {
      res <- fail(sprintf("Unsupported lna_spec '%s'", spec))
      if (!is.null(res)) return(res)
    }
  }

  if (checksum && h5_attr_exists(root, "lna_checksum")) {
    stored <- h5_attr_read(root, "lna_checksum")
    if (is.character(file) && file.exists(file)) {
      calc <- digest::digest(file = file, algo = "sha256")
      if (!identical(calc, stored)) {
        res <- fail("Checksum does not match")
        if (!is.null(res)) return(res)
      }
    } else {
      warning("Checksum requested but file path unavailable; skipping")
    }
  }

  if (h5$exists("transforms")) {
    tf_group <- h5[["transforms"]]
    tf_names <- tf_group$ls()$name
    for (nm in tf_names) {
      desc <- read_json_descriptor(tf_group, nm)
      if (is.list(desc) && !is.null(desc$type)) {
        schema_path <- system.file("schemas", paste0(desc$type, ".schema.json"),
                                  package = utils::packageName())
        if (nzchar(schema_path) && file.exists(schema_path)) {
          json <- jsonlite::toJSON(desc, auto_unbox = TRUE)
          valid <- jsonvalidate::json_validate(json, schema_path, verbose = TRUE)
          if (!isTRUE(valid)) {
            res <- fail(sprintf("Descriptor %s failed schema validation", nm))
            if (!is.null(res)) return(res)
          }
        }
      }
    }
  }

  TRUE
}

