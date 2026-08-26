#' Quick validation of an LNA file
#'
#' @description Performs fast structural validation checks without computing
#' checksums. Verifies that the file exists, is a valid HDF5 file, contains
#' the required /lna group, and has at least one transform descriptor.
#'
#' @param file Path to the `.h5` file to validate.
#'
#' @return Invisibly returns `TRUE` if validation succeeds. Throws an error
#' with class `lna_error_validation` if any check fails.
#' @seealso validate_lna, read_lna
#' @keywords internal
quick_validate_lna <- function(file) {
  # Check 1: File exists
  if (!file.exists(file)) {
    abort_lna(
      sprintf("File does not exist: '%s'", file),
      .subclass = "lna_error_validation",
      location = sprintf("quick_validate_lna:%s", file)
    )
  }

  # Check 2: Valid HDF5 file (can be opened)
  h5 <- tryCatch(
    open_h5(file, mode = "r"),
    error = function(e) {
      abort_lna(
        sprintf("Not a valid HDF5 file: '%s'", file),
        .subclass = "lna_error_validation",
        location = sprintf("quick_validate_lna:%s", file)
      )
    }
  )

  on.exit(close_h5_safely(h5), add = TRUE)

  # Check 3: Has /lna group (legacy) or root-level lna_spec attribute
  root <- h5[["/"]]
  has_lna_group <- h5$exists("lna")
  has_lna_spec <- h5_attr_exists(root, "lna_spec")

  if (!has_lna_group && !has_lna_spec) {
    abort_lna(
      sprintf("File is missing LNA structure (no /lna group or lna_spec attribute): '%s'", file),
      .subclass = "lna_error_validation",
      location = sprintf("quick_validate_lna:%s", file)
    )
  }

  # Check 4: Has /transforms group (descriptors may be empty for edge cases)
  if (!h5$exists("transforms")) {
    abort_lna(
      sprintf("File is missing /transforms group: '%s'", file),
      .subclass = "lna_error_validation",
      location = sprintf("quick_validate_lna:%s", file)
    )
  }

  quick_fail <- function(msg, code = "validation", path = NA_character_,
                         layer = NA_character_,
                         subclass = "lna_error_validation",
                         location = sprintf("quick_validate_lna:%s", file)) {
    abort_lna(msg, .subclass = subclass, location = location,
              code = code, path = path, layer = layer)
  }
  .validate_required_transforms(h5, root, quick_fail)
  .validate_qpp_archives(
    h5, root, quick_fail, mode = "structural",
    sample_tiles = 3L, quick = TRUE
  )

  invisible(TRUE)
}

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
#' @param qpp QPP validation depth: `"all"` decodes every tile, `"sample"`
#'   decodes a deterministic sample, and `"structural"` validates metadata,
#'   index coverage, steps, maps, and receipts without payload decoding.
#' @param qpp_sample_tiles Maximum number of deterministic tiles decoded when
#'   `qpp = "sample"`.
#' @param qpp_limits Optional named resource-cap overrides.
#' @param receipt If `TRUE`, return a reproducible validation receipt instead
#'   of scalar `TRUE` on success.
#' @details
#' When a checksum is present it was computed on the file with the attribute
#' temporarily set to a 64 character placeholder of zeros.  Validation
#' reproduces that state and compares the digest to the stored value.
#'
#' @return `TRUE` if validation succeeds. If `strict = FALSE` and problems are
#'   found, an `lna_validation_report` object is returned.
#' @section Lifecycle:
#' Stable. This is the LNA file-format trust boundary for archives written by
#' `neuroarchive`.
#' @seealso write_lna, read_lna
#' @examples
#' tmp <- tempfile(fileext = ".h5")
#' arr <- array(rnorm(16), dim = c(4, 4, 1, 1))
#' write_lna(arr, tmp, transforms = "quant")
#' validate_lna(tmp)
#' unlink(tmp)
#' @export
validate_lna <- function(file, strict = TRUE, checksum = TRUE,
                         qpp = c("all", "sample", "structural"),
                         qpp_sample_tiles = 16L, qpp_limits = NULL,
                         receipt = FALSE) {
  stopifnot(is.character(file), length(file) == 1)
  qpp <- match.arg(qpp)
  qpp_sample_tiles <- as.numeric(qpp_sample_tiles)
  if (length(qpp_sample_tiles) != 1L || !is.finite(qpp_sample_tiles) ||
      qpp_sample_tiles != trunc(qpp_sample_tiles) || qpp_sample_tiles < 1) {
    abort_lna("qpp_sample_tiles must be a positive exact integer",
              .subclass = "lna_error_validation",
              location = "validate_lna:qpp_sample_tiles")
  }
  if (!is.logical(receipt) || length(receipt) != 1L || is.na(receipt)) {
    abort_lna("receipt must be one TRUE or FALSE value",
              .subclass = "lna_error_validation",
              location = "validate_lna:receipt")
  }

  h5 <- open_h5(file, mode = "r")
  on.exit(close_h5_safely(h5), add = TRUE)
  root <- h5[["/"]]

  issues <- list()

  fail <- function(msg, code = "validation", path = NA_character_,
                   layer = NA_character_, subclass = "lna_error_validation",
                   location = sprintf("validate_lna:%s", file)) {
    if (strict) {
      abort_lna(
        msg,
        .subclass = subclass,
        location = location,
        code = code, path = path, layer = layer
      )
    } else {
      warning(msg, call. = FALSE)
      issues[[length(issues) + 1L]] <<- list(
        level = "error",
        code = code,
        path = path,
        layer = layer,
        message = msg
      )
      invisible(NULL)
    }
  }

  .validate_lna_structure(h5, root, fail)
  .validate_lna_run_geometry(h5, fail)
  .validate_lna_descriptors(h5, root, fail)
  .validate_required_transforms(h5, root, fail)
  qpp_receipts <- .validate_qpp_archives(
    h5, root, fail, mode = qpp,
    sample_tiles = as.integer(qpp_sample_tiles), limits = qpp_limits
  )
  validation_receipt <- .new_lna_validation_receipt(
    file, root, qpp, qpp_receipts
  )

  # Checksum validation is last so byte-level integrity issues cannot mask
  # structural or descriptor-reference failures.
  close_h5_safely(h5)
  h5 <- NULL
  if (checksum) {
    .validate_lna_checksum(file, fail)
  }

  if (length(issues) > 0) {
    report <- .new_lna_validation_report(file, issues)
    report$receipt <- validation_receipt
    return(report)
  }

  if (isTRUE(receipt)) validation_receipt else TRUE
}

.new_lna_validation_receipt <- function(file, root, qpp_mode, qpp_receipts) {
  required <- if (h5_attr_exists(root, "required_transforms")) {
    as.character(h5_attr_read(root, "required_transforms"))
  } else {
    character()
  }
  checksum <- if (h5_attr_exists(root, "lna_checksum")) {
    h5_attr_read(root, "lna_checksum")
  } else {
    NULL
  }
  structure(list(
    receipt_version = "1.0.0",
    file_size_bytes = unname(file.info(file)$size),
    lna_checksum = checksum,
    required_transforms = required,
    qpp_validation_mode = qpp_mode,
    qpp = qpp_receipts,
    validator = list(
      package = "neuroarchive",
      version = tryCatch(
        as.character(utils::packageVersion("neuroarchive")),
        error = function(e) "0.1.1"
      )
    )
  ), class = c("lna_validation_receipt", "list"))
}

.new_lna_validation_report <- function(file, issues) {
  if (length(issues) == 0L) {
    issue_df <- data.frame(
      level = character(),
      code = character(),
      path = character(),
      layer = character(),
      message = character(),
      stringsAsFactors = FALSE
    )
  } else {
    issue_df <- do.call(
      rbind,
      lapply(issues, function(x) {
        data.frame(
          level = x$level %||% "error",
          code = x$code %||% "validation",
          path = x$path %||% NA_character_,
          layer = x$layer %||% NA_character_,
          message = x$message %||% "",
          stringsAsFactors = FALSE
        )
      })
    )
    rownames(issue_df) <- NULL
  }

  structure(
    list(
      ok = nrow(issue_df) == 0L,
      file = file,
      issues = issue_df
    ),
    class = c("lna_validation_report", "list")
  )
}

#' @export
as.character.lna_validation_report <- function(x, ...) {
  x$issues$message
}

#' @export
print.lna_validation_report <- function(x, ...) {
  status <- if (isTRUE(x$ok)) "OK" else "FAILED"
  cat(sprintf("LNA validation report: %s\n", status))
  cat(sprintf("File: %s\n", x$file))
  cat(sprintf("Issues: %d\n", nrow(x$issues)))
  if (nrow(x$issues) > 0L) {
    utils::str(head(x$issues, 10L))
  }
  invisible(x)
}

.validate_lna_structure <- function(h5, root, fail) {
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

  for (grp in c("transforms", "basis", "scans")) {
    if (!h5$exists(grp)) {
      fail(sprintf("Required group '%s' missing", grp))
      next
    }
    obj <- h5[[grp]]
    if (!inherits(obj, "H5Group")) {
      fail(sprintf("Required group '%s' is not an HDF5 group", grp))
    }
    safe_h5_close(obj)
  }

  invisible(TRUE)
}

.validate_lna_run_geometry <- function(h5, fail) {
  if (!h5$exists("header")) return(invisible(TRUE))
  header_group <- h5[["header"]]
  on.exit(safe_h5_close(header_group), add = TRUE)
  policy <- if (h5_attr_exists(header_group, "geometry_policy")) {
    as.character(h5_attr_read(header_group, "geometry_policy"))
  } else {
    "legacy"
  }

  scan_ids <- discover_run_ids(h5)
  if (!header_group$exists("runs")) {
    if (identical(policy, "fmri")) {
      fail(
        "fMRI archive is missing /header/runs geometry",
        code = "run_geometry_missing",
        path = "/header/runs",
        layer = "geometry"
      )
    }
    if (header_group$exists("global")) {
      global <- header_group[["global"]]
      on.exit(safe_h5_close(global), add = TRUE)
      legacy <- .read_lna_header_group(global)
      if (.lna_header_has_geometry(legacy)) {
        error <- tryCatch({
          .lna_normalize_run_geometry(
            legacy,
            if (length(scan_ids) == 1L) scan_ids[[1L]] else "legacy-global"
          )
          NULL
        }, error = identity)
        if (inherits(error, "error")) {
          fail(
            paste0("legacy global geometry is malformed: ", conditionMessage(error)),
            code = "run_geometry_malformed",
            path = "/header/global",
            layer = "geometry"
          )
        }
      }
    }
    return(invisible(TRUE))
  }

  runs_group <- header_group[["runs"]]
  on.exit(safe_h5_close(runs_group), add = TRUE)
  geometry_ids <- names(runs_group)
  missing <- setdiff(scan_ids, geometry_ids)
  extra <- setdiff(geometry_ids, scan_ids)
  if (length(missing)) {
    fail(
      paste0("run geometry is missing for: ", paste(missing, collapse = ", ")),
      code = "run_geometry_missing",
      path = "/header/runs",
      layer = "geometry"
    )
  }
  if (length(extra)) {
    fail(
      paste0("run geometry has no matching scan for: ", paste(extra, collapse = ", ")),
      code = "run_geometry_orphan",
      path = "/header/runs",
      layer = "geometry"
    )
  }

  validated <- list()
  for (run_id in intersect(scan_ids, geometry_ids)) {
    group <- runs_group[[run_id]]
    geometry <- .read_lna_header_group(group)
    safe_h5_close(group)
    required <- c(
      "geometry_schema_version", "run_id", "dims", "spacing", "origin",
      "transform", "spatial_units", "voxel_ordering", "n_volumes", "grid_id"
    )
    absent <- setdiff(required, names(geometry))
    if (length(absent)) {
      fail(
        paste0(
          "run geometry for ", run_id, " is missing: ",
          paste(absent, collapse = ", ")
        ),
        code = "run_geometry_missing_field",
        path = paste0("/header/runs/", run_id),
        layer = "geometry"
      )
      next
    }
    if (!identical(as.character(geometry$geometry_schema_version),
                   .lna_run_geometry_schema_version) ||
        !identical(as.character(geometry$run_id), run_id)) {
      fail(
        paste0("run geometry identity/schema mismatch for ", run_id),
        code = "run_geometry_identity",
        path = paste0("/header/runs/", run_id),
        layer = "geometry"
      )
      next
    }
    normalized <- tryCatch(
      .lna_normalize_run_geometry(geometry, run_id),
      error = identity
    )
    if (inherits(normalized, "error")) {
      fail(
        paste0(
          "run geometry for ", run_id, " is malformed: ",
          conditionMessage(normalized)
        ),
        code = "run_geometry_malformed",
        path = paste0("/header/runs/", run_id),
        layer = "geometry"
      )
      next
    }
    validated[[run_id]] <- normalized
  }

  if (length(validated) > 1L && header_group$exists("global")) {
    global_group <- header_group[["global"]]
    global <- .read_lna_header_group(global_group)
    safe_h5_close(global_group)
    if (.lna_header_has_geometry(global)) {
      grids <- vapply(validated, `[[`, character(1), "grid_id")
      if (length(unique(grids)) != 1L) {
        fail(
          "global spatial geometry is ambiguous because run grids differ",
          code = "run_geometry_ambiguous_global",
          path = "/header/global",
          layer = "geometry"
        )
      } else {
        normalized_global <- tryCatch(
          .lna_normalize_run_geometry(global, names(validated)[[1L]]),
          error = identity
        )
        if (inherits(normalized_global, "error") ||
            !identical(normalized_global$grid_id, grids[[1L]])) {
          fail(
            "global spatial geometry does not match the run-scoped grid",
            code = "run_geometry_global_mismatch",
            path = "/header/global",
            layer = "geometry"
          )
        }
      }
    }
  }
  invisible(TRUE)
}

.validate_lna_descriptors <- function(h5, root, fail) {
  if (!h5$exists("transforms")) {
    return(invisible(TRUE))
  }

  tf_group <- h5[["transforms"]]
  on.exit(safe_h5_close(tf_group), add = TRUE)

  tf_names <- tf_group$ls()$name
  tf_names <- tf_names[grepl("\\.json$", tf_names)]
  tf_names <- tf_names[!grepl("_report\\.json$", tf_names)]

  for (nm in tf_names) {
    desc <- tryCatch(
      read_json_descriptor(tf_group, nm),
      error = function(e) e
    )
    if (inherits(desc, "error")) {
      fail(conditionMessage(desc))
      next
    }
    if (!is.list(desc) || is.null(desc$type) ||
        !is.character(desc$type) || length(desc$type) != 1 || !nzchar(desc$type)) {
      fail(sprintf("Descriptor %s is missing a scalar transform type", nm))
      next
    }
    if (!is.null(desc$run_id)) {
      valid_run_id <- is.character(desc$run_id) && length(desc$run_id) == 1L &&
        grepl("^run-[0-9]{2}$", desc$run_id)
      if (!valid_run_id) {
        fail(sprintf("Descriptor %s has an invalid run_id", nm))
      } else if (!h5$exists(paste0("scans/", desc$run_id))) {
        fail(sprintf(
          "Descriptor %s references missing run group /scans/%s",
          nm,
          desc$run_id
        ))
      }
    }

    if (identical(desc$type, "qpp") && !identical(desc$version, "1.1.0")) {
      fail(
        paste0("unsupported required qpp transform version: ", desc$version %||% "<missing>"),
        code = if (identical(desc$version, "1.0.0")) {
          "qpp_legacy_profile"
        } else {
          "required_transform"
        },
        path = nm,
        layer = "decoder",
        subclass = "lna_error_required_transform",
        location = "validate_lna:qpp:decoder_version"
      )
    }
    if (desc$type %in% c("qpp", "parcel.mean")) {
      contract_error <- tryCatch({
        expected <- if (identical(desc$type, "qpp")) {
          .qpp_analysis_contract(desc$params %||% list())
        } else {
          .parcel_mean_analysis_contract(desc$params %||% list())
        }
        .analysis_assert_contract_match(
          desc$analysis_contract, expected,
          paste0(desc$type, "@", desc$version %||% "<missing>")
        )
        NULL
      }, error = function(error) error)
      if (inherits(contract_error, "error")) {
        subclass <- class(contract_error)[grepl("^lna_error_", class(contract_error))][1L]
        if (length(subclass) != 1L || is.na(subclass)) {
          subclass <- "lna_error_analysis_contract"
        }
        fail(
          conditionMessage(contract_error),
          code = contract_error$code %||% "analysis_contract",
          path = nm,
          layer = "descriptor",
          subclass = subclass,
          location = contract_error$location %||% paste0(
            "validate_lna:", desc$type, ":analysis_contract"
          )
        )
      }
    }

    .validate_lna_descriptor_schema(desc, nm, fail)
    .validate_lna_dataset_references(h5, root, desc, nm, fail)
  }

  invisible(TRUE)
}

.validate_lna_descriptor_schema <- function(desc, nm, fail) {
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
    return(invisible(FALSE))
  }

  json <- jsonlite::toJSON(desc, auto_unbox = TRUE, null = "null")
  valid <- jsonvalidate::json_validate(json, schema_path,
                                       verbose = TRUE, engine = "ajv")
  if (!isTRUE(valid)) {
    errors <- attr(valid, "errors")
    error_message <- sprintf("Descriptor %s failed schema validation", nm)
    if (!is.null(errors) && length(errors) > 0) {
      error_details <- paste(utils::capture.output(print(errors)), collapse = "\n")
      error_message <- paste(error_message, "Details:", error_details, sep = "\n")
    }
    fail(error_message)
  }

  invisible(TRUE)
}

.validate_lna_dataset_references <- function(h5, root, desc, nm, fail) {
  if (is.null(desc$datasets)) {
    return(invisible(TRUE))
  }

  for (ds in desc$datasets) {
    path <- ds$path
    if (is.null(path) || !nzchar(path)) next
    if (!path_exists_safely(h5, path)) {
      fail(sprintf("Dataset '%s' referenced in %s missing", path, nm))
      next
    }

    dset <- h5[[path]]
    on.exit(safe_h5_close(dset), add = TRUE)
    if (!inherits(dset, "H5D")) {
      fail(sprintf("Path '%s' referenced in %s is not an HDF5 dataset", path, nm))
      next
    }

    if (!is.null(ds$dims)) {
      if (!identical(as.integer(ds$dims), as.integer(dset$dims))) {
        fail(sprintf("Dimensions mismatch for dataset '%s'", path))
      }
    }

    if (!is.null(ds$dtype)) {
      dt <- dset$get_type()
      on.exit(if (inherits(dt, "H5T")) dt$close(), add = TRUE)
      actual <- .h5_dtype_name(dt)
      compatible <- identical(tolower(ds$dtype), actual) ||
        (identical(tolower(ds$dtype), "qpp_index_v1") && .h5_is_qpp_index_v1(dt))
      if (!compatible) {
        fail(sprintf("Dtype mismatch for dataset '%s'", path))
      }
    }

    qpp_bounded_role <- identical(desc$type, "qpp") &&
      !is.null(ds$role) && ds$role %in% c(
        "qpp_payload", "qpp_tile_index", "qpp_parcel_steps",
        "qpp_quality_receipt"
      )
    data <- if (qpp_bounded_role) {
      NULL
    } else {
      tryCatch(
        h5_read(root, path),
        error = function(e) {
          fail(sprintf("Error reading dataset '%s': %s", path, e$message))
          NULL
        }
      )
    }
    if (is.numeric(data)) {
      is_special_dataset <- !is.null(ds$role) &&
        ds$role %in% c("singular_values", "basis_matrix", "temporal_basis")
      is_empty_data <- length(data) == 0
      if (!is_special_dataset && !is_empty_data &&
          (all(is.na(data)) || all(data == 0))) {
        fail(sprintf("Dataset '%s' contains only zeros/NaN", path))
      }
    }
  }

  invisible(TRUE)
}

.h5_dtype_name <- function(dt) {
  class_id <- dt$get_class()
  size <- dt$get_size()
  class_name <- as.character(class_id)
  switch(class_name,
    H5T_INTEGER = paste0(
      ifelse(identical(as.character(dt$get_sign()), "H5T_SGN_NONE"), "u", ""),
      "int", size * 8
    ),
    H5T_FLOAT = paste0("float", size * 8),
    H5T_COMPOUND = "compound",
    `0` = paste0("float", size * 8),
    `1` = paste0(ifelse(identical(as.character(dt$get_sign()), "H5T_SGN_NONE"),
                        "u", ""), "int", size * 8),
    "unknown"
  )
}

.h5_is_qpp_index_v1 <- function(dt) {
  if (!identical(as.character(dt$get_class()), "H5T_COMPOUND")) return(FALSE)
  labels <- c(
    "tile_id", "run_ordinal", "unit_ordinal", "time_start", "time_count",
    "parcel_start", "parcel_count", "byte_offset", "byte_length",
    "uncompressed_length", "inner_header_offset", "predictor_table_offset",
    "outer_codec", "flags", "digest_row"
  )
  widths <- c(4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 4, 4, 1, 2, 4)
  offsets <- c(0, 4, 8, 12, 16, 20, 24, 32, 40, 48, 56, 60, 64, 66, 68)
  if (!identical(dt$get_cpd_labels(), labels) ||
      !identical(as.numeric(dt$get_cpd_offsets()), as.numeric(offsets)) ||
      dt$get_size() != 72) return(FALSE)
  types <- dt$get_cpd_types()
  on.exit(lapply(types, safe_h5_close), add = TRUE)
  sizes <- vapply(types, function(type) type$get_size(), numeric(1))
  classes <- vapply(types, function(type) as.character(type$get_class()), character(1))
  signs <- vapply(types, function(type) as.character(type$get_sign()), character(1))
  identical(as.numeric(sizes), as.numeric(widths)) &&
    all(classes == "H5T_INTEGER") && all(signs == "H5T_SGN_NONE")
}

.validate_lna_checksum <- function(file, fail) {
  h5 <- open_h5(file, mode = "r")
  on.exit(close_h5_safely(h5), add = TRUE)
  root <- h5[["/"]]

  if (!h5_attr_exists(root, "lna_checksum")) {
    return(invisible(TRUE))
  }

  stored_checksum_value <- h5_attr_read(root, "lna_checksum")
  placeholder_checksum <- paste(rep("0", 64), collapse = "")
  is_valid_hash <- nchar(stored_checksum_value) == 64L &&
    grepl("^[0-9a-f]{64}$", stored_checksum_value)

  if (!is_valid_hash) {
    fail("Checksum does not match")
    return(invisible(FALSE))
  }

  close_h5_safely(h5)
  h5 <- NULL

  temp_copy_path <- tempfile(fileext = ".h5")
  on.exit(if (file.exists(temp_copy_path)) unlink(temp_copy_path), add = TRUE)
  if (!file.copy(file, temp_copy_path, overwrite = TRUE)) {
    fail("Checksum does not match")
    return(invisible(FALSE))
  }

  calculated_checksum_on_copy <- NULL
  tryCatch({
    raw_bytes <- readBin(temp_copy_path, "raw", file.info(temp_copy_path)$size)
    stored_raw <- charToRaw(stored_checksum_value)
    placeholder_raw <- charToRaw(placeholder_checksum)
    pos <- grepRaw(stored_raw, raw_bytes, fixed = TRUE)
    if (length(pos) == 1L) {
      raw_bytes[pos:(pos + 63L)] <- placeholder_raw
      calculated_checksum_on_copy <- digest::digest(
        raw_bytes,
        algo = "sha256",
        serialize = FALSE
      )
    }
  }, error = function(e) {
    calculated_checksum_on_copy <<- NULL
  })

  if (is.null(calculated_checksum_on_copy) ||
      !identical(calculated_checksum_on_copy, stored_checksum_value)) {
    fail("Checksum does not match")
  }

  invisible(TRUE)
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
