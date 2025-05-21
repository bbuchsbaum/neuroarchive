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
#' @seealso write_lna, read_lna
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

  for (attr_name in c("creator", "required_transforms")) {
    if (!h5_attr_exists(root, attr_name)) {
      fail(sprintf("Missing %s attribute", attr_name))
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

  for (grp in c("transforms", "basis", "scans")) {
    if (!h5$exists(grp)) {
      fail(sprintf("Required group '%s' missing", grp))
    }
  }

  optional_groups <- c("spatial", "plugins")
  for (grp in optional_groups) {
    if (h5$exists(grp)) {
      NULL  # presence noted but no action; placeholder for future checks
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

        if (!is.null(desc$datasets)) {
          for (ds in desc$datasets) {
            path <- ds$path
            if (is.null(path) || !nzchar(path)) next
            if (!h5$exists(path)) {
              fail(sprintf("Dataset '%s' referenced in %s missing", path, nm))
              next
            }

            dset <- h5[[path]]
            if (!is.null(ds$dims)) {
              if (!identical(as.integer(ds$dims), as.integer(dset$dims))) {
                fail(sprintf("Dimensions mismatch for dataset '%s'", path))
              }
            }

            if (!is.null(ds$dtype)) {
              dt <- dset$get_type()
              class_id <- dt$get_class()
              size <- dt$get_size()
              actual <- switch(as.character(class_id),
                `1` = paste0(ifelse(dt$get_sign() == "H5T_SGN_NONE", "u", ""),
                              "int", size * 8),
                `0` = paste0("float", size * 8),
                "unknown" )
              if (!identical(tolower(ds$dtype), actual)) {
                fail(sprintf("Dtype mismatch for dataset '%s'", path))
              }
            }

            data <- tryCatch(
              h5_read(root, path),
              error = function(e) {
                fail(sprintf("Error reading dataset '%s': %s", path, e$message))
                NULL
              }
            )
            if (is.numeric(data)) {
              if (all(is.na(data)) || all(data == 0)) {
                fail(sprintf("Dataset '%s' contains only zeros/NaN", path))
              }
            }
          }
        }
      }
    }
  }

  if (length(issues) > 0) {
    return(issues)
  }

  TRUE
}

#' Runtime validation for a transform step
#'
#' Checks dataset paths referenced in a descriptor and verifies that all
#' required parameters are present before a transform is executed.
#'
#' @param type Transform type name.
#' @param desc Descriptor list parsed from JSON.
#' @param h5 An open `H5File` object.
#' @return Invisibly `TRUE` or throws an error on validation failure.
#' @keywords internal
runtime_validate_step <- function(type, desc, h5) {
  stopifnot(is.character(type), length(type) == 1)
  stopifnot(is.list(desc))
  stopifnot(inherits(h5, "H5File"))

  root <- h5[["/"]]
  if (!is.null(desc$datasets)) {
    for (ds in desc$datasets) {
      if (!is.null(ds$path)) {
        assert_h5_path(root, ds$path)
      }
    }
  }

  req <- required_params(type)
  params <- desc$params %||% list()
  missing <- setdiff(req, names(params))
  if (length(missing) > 0) {
    abort_lna(
      paste0(
        "Descriptor for transform '", type,
        "' missing required parameter(s): ",
        paste(missing, collapse = ", ")
      ),
      .subclass = "lna_error_descriptor",
      location = sprintf("runtime_validate_step:%s", type)
    )
  }

  invisible(TRUE)
}

