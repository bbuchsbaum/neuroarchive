#' Validate an LNA file
#'
#' @description Basic validator that checks the LNA specification version,
#' optional SHA256 checksum and, if available, validates transform
#' descriptors against their JSON schemas.
#'
#' @param file Path to the `.h5` file to validate.
#' @param strict Logical. If `TRUE` (default) validation failures abort with
#'   `lna_error_validation`. If `FALSE`, all validation issues are collected and
#'   returned. A warning is issued for each problem found.
#' @param checksum Logical. If `TRUE` (default) verify the `lna_checksum`
#'   attribute when present.
#'
#' @return `TRUE` if validation succeeds. If `strict = FALSE` and problems are
#'   found, a character vector of issue messages is returned instead.
#' @examples
#' validate_lna("example.lna.h5")
#' @export
validate_lna <- function(file, strict = TRUE, checksum = TRUE) {
  stopifnot(is.character(file), length(file) == 1)

  h5 <- open_h5(file, mode = "r")
  on.exit(close_h5_safely(h5))
  root <- h5[["/"]]

  issues <- character()

  fail <- function(msg) {
    if (strict) {
      abort_lna(
        msg,
        .subclass = "lna_error_validation",
        location = sprintf("validate_lna:%s", file)
      )
    } else {
      warning(msg)
      issues <<- c(issues, msg)
      invisible(NULL)
    }
  }

  if (!h5_attr_exists(root, "lna_spec")) {
    fail("Missing lna_spec attribute")
  } else {
    spec <- h5_attr_read(root, "lna_spec")
    if (!identical(spec, "LNA R v2.0")) {
      fail(sprintf("Unsupported lna_spec '%s'", spec))
    }
  }

  if (checksum && h5_attr_exists(root, "lna_checksum")) {
    stored <- h5_attr_read(root, "lna_checksum")
    if (is.character(file) && file.exists(file)) {
      calc <- digest::digest(file = file, algo = "sha256")
      if (!identical(calc, stored)) {
        fail("Checksum does not match")
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
        pkgs <- unique(c("neuroarchive", loadedNamespaces()))
        schema_path <- ""
        for (pkg in pkgs) {
          path <- system.file(
            "schemas",
            paste0(desc$type, ".schema.json"),
            package = pkg
          )
          if (nzchar(path) && file.exists(path)) {
            schema_path <- path
            break
          }
        }

        if (!nzchar(schema_path)) {
          fail(sprintf("Schema for transform '%s' not found", desc$type))
          next
        }

        json <- jsonlite::toJSON(desc, auto_unbox = TRUE)
        valid <- jsonvalidate::json_validate(json, schema_path, verbose = TRUE)
        if (!isTRUE(valid)) {
          fail(sprintf("Descriptor %s failed schema validation", nm))
        }
      }
    }
  }

  if (length(issues) > 0) {
    return(issues)
  }

  TRUE
}

