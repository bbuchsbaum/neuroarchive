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
#' When a checksum is present it was computed on the file with the attribute
#' temporarily set to a 64 character placeholder of zeros.  Validation
#' reproduces that state and compares the digest to the stored value.
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
  on.exit(neuroarchive:::close_h5_safely(h5))
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
      warning(msg, call. = FALSE)
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
    stored_checksum_value <- h5_attr_read(root, "lna_checksum")
    # To validate correctly, we need the hash of the file *with* the placeholder in place of the actual checksum.
    # This matches how materialise_plan calculates the hash initially.
    
    # Create the same placeholder used by materialise_plan
    placeholder_checksum <- paste(rep("0", 64), collapse = "")
    
    # Close the original file, make a copy, overwrite the checksum with the placeholder in the copy, then hash
    current_file_path <- h5$filename
    neuroarchive:::close_h5_safely(h5) # Close the original file handle

    temp_copy_path <- tempfile(fileext = ".h5")
    file.copy(current_file_path, temp_copy_path, overwrite = TRUE)

    h5_temp_copy <- NULL
    calculated_checksum_on_copy <- NULL
    
    tryCatch({
      h5_temp_copy <- open_h5(temp_copy_path, mode = "r+")
      root_temp_copy <- h5_temp_copy[["/"]]
      if (h5_attr_exists(root_temp_copy, "lna_checksum")) {
        h5_attr_delete(root_temp_copy, "lna_checksum")
      }
      # Now, lna_checksum attribute is guaranteed not to exist or has been deleted.
      # Write the placeholder anew.
      h5_attr_write(root_temp_copy, "lna_checksum", placeholder_checksum)
      
      # Important: close the temp file *before* hashing it
      neuroarchive:::close_h5_safely(h5_temp_copy)
      h5_temp_copy <- NULL # Mark as closed for finally block
      
      calculated_checksum_on_copy <- digest::digest(file = temp_copy_path, algo = "sha256")
    }, finally = {
      if (!is.null(h5_temp_copy) && inherits(h5_temp_copy, "H5File") && h5_temp_copy$is_valid) {
        neuroarchive:::close_h5_safely(h5_temp_copy)
      }
      if (file.exists(temp_copy_path)) {
        unlink(temp_copy_path)
      }
    })

    # Reopen the original file for subsequent validation steps if any
    h5 <- open_h5(file, mode = "r") # Re-open the original file
    root <- h5[["/"]] # Re-assign root based on the new h5 handle

    if (!is.null(calculated_checksum_on_copy) && !identical(calculated_checksum_on_copy, stored_checksum_value)) {
      fail("Checksum does not match")
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
        valid <- jsonvalidate::json_validate(json, schema_path, verbose = TRUE, engine = "ajv")
        if (!isTRUE(valid)) {
          errors <- attr(valid, "errors")
          error_message <- sprintf("Descriptor %s failed schema validation", nm)
          if (!is.null(errors) && length(errors) > 0) {
            error_details <- paste(utils::capture.output(print(errors)), collapse = "\n")
            error_message <- paste(error_message, "Details:", error_details, sep = "\n")
          }
          fail(error_message)
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
            on.exit(if (inherits(dset, "H5D")) dset$close(), add = TRUE)
            if (!is.null(ds$dims)) {
              if (!identical(as.integer(ds$dims), as.integer(dset$dims))) {
                fail(sprintf("Dimensions mismatch for dataset '%s'", path))
              }
            }

            if (!is.null(ds$dtype)) {
              dt <- dset$get_type()
              on.exit(if (inherits(dt, "H5T")) dt$close(), add = TRUE)
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
              # Skip zero/NaN check for datasets that can legitimately be zero
              is_special_dataset <- !is.null(ds$role) && ds$role %in% c("singular_values", "basis_matrix", "temporal_basis")
              # Also skip if the data is genuinely empty (e.g. 0xN or Nx0 matrix)
              is_empty_data <- length(data) == 0
              if (!is_special_dataset && !is_empty_data && (all(is.na(data)) || all(data == 0))) {
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
  
  # Check both top-level descriptor keys and params sublist for required parameters
  all_keys <- unique(c(names(desc), names(params)))
  missing <- setdiff(req, all_keys)
  
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

