# Fixed parcel-mean spatial transform.

#' Forward step for parcel.mean
#'
#' Reduces a time-by-voxel input to deterministic time-by-parcel means while
#' persisting only the spatial reconstruction assets. The parcel series remains
#' in the transform stash for a following codec such as `qpp`.
#'
#' @param type Transform type string.
#' @param desc Transform descriptor.
#' @param handle DataHandle for the current run.
#' @return Updated DataHandle.
#' @method forward_step parcel.mean
#' @export
#' @keywords internal
forward_step.parcel.mean <- function(type, desc, handle) {
  p <- desc$params %||% list()
  .parcel_mean_validate_provenance(p, "forward_step.parcel.mean")

  if (is.null(handle$mask_info) || is.null(handle$mask_info$mask)) {
    abort_lna("parcel.mean requires an explicit spatial mask",
              .subclass = "lna_error_validation",
              location = "forward_step.parcel.mean:mask")
  }
  inp <- handle$pull_first(c(desc$inputs[[1]] %||% "input", "input_dense_mat",
                             "dense_mat", "input"))
  if (is.list(inp$value) && !is.array(inp$value)) {
    abort_lna(
      "parcel.mean accepts one run at a time; encode multiple runs independently",
      .subclass = "lna_error_validation",
      location = "forward_step.parcel.mean:runs"
    )
  }

  mask <- as.logical(as.array(handle$mask_info$mask))
  mask_dims <- dim(as.array(handle$mask_info$mask))
  dim(mask) <- mask_dims
  if (!any(mask)) {
    abort_lna("parcel.mean mask has no active voxels",
              .subclass = "lna_error_validation",
              location = "forward_step.parcel.mean:mask")
  }
  mask_checksum <- .compute_mask_checksum(mask)
  map <- .parcel_mean_active_vector(
    p$labels,
    mask,
    p$map_mask_checksum,
    mask_checksum,
    "labels"
  )
  .qpp_assert_source_labels(map, "forward_step.parcel.mean:labels")

  source_ids <- sort(unique(as.numeric(map)))
  if (!is.null(p$source_labels)) {
    expected <- sort(unique(as.numeric(p$source_labels)))
    .qpp_assert_source_labels(expected, "forward_step.parcel.mean:source_labels")
    if (!identical(source_ids, expected)) {
      abort_lna(
        "source_labels must exactly match represented labels; empty parcels are forbidden",
        .subclass = "lna_error_validation",
        location = "forward_step.parcel.mean:empty_parcel"
      )
    }
  }
  dense_labels <- match(map, source_ids)
  k <- length(source_ids)
  label_dtype <- if (k <= 65535L) "uint16" else "uint32"

  weights_supplied <- !is.null(p$weights)
  weights <- if (weights_supplied) {
    .parcel_mean_active_vector(
      p$weights,
      mask,
      p$map_mask_checksum,
      mask_checksum,
      "weights"
    )
  } else {
    rep(1, length(dense_labels))
  }
  if (!is.numeric(weights) || any(!is.finite(weights)) || any(weights < 0)) {
    abort_lna("parcel weights must be finite and non-negative",
              .subclass = "lna_error_validation",
              location = "forward_step.parcel.mean:weights")
  }
  counts <- tabulate(dense_labels, nbins = k)
  weight_sums <- as.numeric(rowsum(weights, dense_labels, reorder = TRUE))
  if (any(counts == 0L) || any(weight_sums <= 0)) {
    abort_lna("every parcel must contain voxels with positive total weight",
              .subclass = "lna_error_validation",
              location = "forward_step.parcel.mean:empty_or_zero_weight")
  }

  vox_time <- convert_to_masked_vox_time_matrix(inp$value, mask)
  if (nrow(vox_time) != length(dense_labels)) {
    abort_lna("parcel map and masked input voxel counts differ",
              .subclass = "lna_error_validation",
              location = "forward_step.parcel.mean:dims")
  }
  if (!is.numeric(vox_time) || any(!is.finite(vox_time))) {
    abort_lna("parcel.mean input must contain only finite numeric values",
              .subclass = "lna_error_validation",
              location = "forward_step.parcel.mean:finite")
  }
  weighted <- sweep(vox_time, 1L, weights, "*")
  sums <- rowsum(weighted, dense_labels, reorder = TRUE)
  parcel_series <- t(sweep(sums, 1L, weight_sums, "/"))
  dimnames(parcel_series) <- list(NULL, as.character(source_ids))

  geometry_identity <- .parcel_mean_validate_affine(
    p,
    handle$meta$header %||% list()
  )
  map_checksum <- .parcel_mean_map_checksum(dense_labels, source_ids, mask_checksum)
  run_id <- sanitize_run_id(handle$current_run_id %||% "run-01")
  plan <- handle$plan
  fname <- plan$get_next_filename(type)
  step_index <- plan$next_index
  output_key <- desc$outputs[[1]] %||% "parcel_series"
  weight_checksum <- if (weights_supplied) {
    .parcel_mean_numeric_checksum(weights, "lna_parcel_weights_v1")
  } else {
    NULL
  }
  asset_id <- .qpp_canonical_checksum(list(
    map_checksum = map_checksum,
    mask_checksum = mask_checksum,
    grid_id = geometry_identity$grid_id,
    affine_checksum = geometry_identity$affine_id,
    voxel_ordering = p$voxel_ordering,
    spatial_units = p$spatial_units,
    weight_policy = if (weights_supplied) "fixed-voxel-weights" else "uniform",
    weight_checksum = weight_checksum,
    source_id = p$source_id,
    fit_scope = p$fit_scope,
    calibration_checksum = p$calibration_checksum %||% NULL
  ), "lna_parcel_asset_v1")
  base <- paste0(
    "/spatial/qpp_assets/", sub("^sha256:", "", asset_id)
  )
  paths <- list(
    labels = paste0(base, "/parcel_labels"),
    source = paste0(base, "/source_parcel_ids"),
    counts = paste0(base, "/parcel_counts"),
    weights = paste0(base, "/parcel_weights")
  )

  p_record <- p
  p_record$labels <- NULL
  p_record$weights <- NULL
  p_record$source_labels <- NULL
  p_record$matrix_order <- "time-by-parcel"
  p_record$parcel_count <- k
  p_record$label_dtype <- label_dtype
  p_record$weight_policy <- if (weights_supplied) "fixed-voxel-weights" else "uniform"
  p_record$weight_checksum <- weight_checksum
  p_record$mask_checksum <- mask_checksum
  p_record$affine_checksum <- geometry_identity$affine_id
  p_record$grid_id <- geometry_identity$grid_id
  p_record$map_checksum <- map_checksum
  p_record$asset_id <- asset_id
  p_record$asset_path <- base
  p_record$step_index <- as.integer(step_index)
  p_record$lossy_support <- "parcel-at-same-TR"
  p_record$coding_dependency <- "none"
  p_record$normalization_provenance <- .analysis_normalization_record(
    p$normalization_provenance
  )
  params_json <- as.character(jsonlite::toJSON(p_record, auto_unbox = TRUE,
                                                null = "null"))

  datasets <- list(
    list(path = paths$labels, role = "parcel_labels", dtype = label_dtype,
         dims = length(dense_labels)),
    list(path = paths$source, role = "source_parcel_ids", dtype = "uint64",
         dims = k),
    list(path = paths$counts, role = "parcel_counts", dtype = "uint32", dims = k)
  )
  if (weights_supplied) {
    datasets[[length(datasets) + 1L]] <- list(
      path = paths$weights, role = "parcel_weights", dtype = "float64",
      dims = length(weights)
    )
  }
  desc$version <- "1.1.0"
  desc$run_id <- run_id
  desc$params <- p_record
  desc$analysis_contract <- .parcel_mean_analysis_contract(p_record)
  desc$inputs <- c(inp$key)
  desc$outputs <- c(output_key)
  desc$datasets <- datasets
  plan$add_descriptor(fname, desc)

  .parcel_mean_add_asset(plan, paths$labels, dense_labels, label_dtype,
                         "parcel_labels", run_id, step_index, params_json)
  .parcel_mean_add_asset(plan, paths$source, source_ids, "uint64",
                         "source_parcel_ids", run_id, step_index, params_json)
  .parcel_mean_add_asset(plan, paths$counts, counts, "uint32",
                         "parcel_counts", run_id, step_index, params_json)
  if (weights_supplied) {
    .parcel_mean_add_asset(plan, paths$weights, weights, "float64",
                           "parcel_weights", run_id, step_index, params_json)
  }

  meta <- handle$meta
  meta$parcel_mean_assets <- list(
    dense_labels = as.integer(dense_labels), source_ids = source_ids,
    counts = counts, weights = if (weights_supplied) weights else NULL,
    descriptor = fname, transform = "parcel.mean@1.1.0",
    run_id = run_id, step_index = as.integer(step_index),
    asset_id = asset_id, asset_path = base, labels_path = paths$labels,
    mask_checksum = mask_checksum,
    affine_checksum = geometry_identity$affine_id,
    grid_id = geometry_identity$grid_id,
    map_checksum = map_checksum, signal_units = p$signal_units,
    spatial_units = p$spatial_units,
    normalization_provenance = p_record$normalization_provenance
  )
  handle <- handle$with(meta = meta)
  handle$update_stash(inp$key, setNames(list(parcel_series), output_key))
}

#' Inverse step for parcel.mean
#'
#' @param type Transform type string.
#' @param desc Transform descriptor.
#' @param handle DataHandle for the current run.
#' @return Updated DataHandle containing a time-by-active-voxel matrix.
#' @method invert_step parcel.mean
#' @export
#' @keywords internal
invert_step.parcel.mean <- function(type, desc, handle) {
  p <- desc$params %||% list()
  .parcel_mean_validate_provenance(p, "invert_step.parcel.mean")
  parcel_key <- desc$outputs[[1]] %||% "parcel_series"
  input_key <- desc$inputs[[1]] %||% "input"
  if (!handle$has_key(parcel_key)) {
    abort_lna(
      "parcel series is absent; parcel.mean stores only reconstruction assets and requires a downstream codec",
      .subclass = "lna_error_contract",
      location = "invert_step.parcel.mean:parcel_series"
    )
  }
  parcel_series <- as.matrix(handle$get_inputs(parcel_key)[[1]])
  assets <- handle$meta$parcel_mean_assets
  applied <- handle$subset_applied %||% list(roi_mask = FALSE, time_idx = FALSE)
  time_idx <- handle$subset$time_idx %||% handle$subset$time
  if (!isTRUE(applied$time_idx) && !is.null(time_idx)) {
    idx <- as.integer(time_idx)
    if (anyNA(idx) || any(idx < 1L) || any(idx > nrow(parcel_series))) {
      abort_lna("time_idx is outside the parcel series",
                .subclass = "lna_error_validation",
                location = "invert_step.parcel.mean:time_idx")
    }
    parcel_series <- parcel_series[idx, , drop = FALSE]
    handle$subset_applied$time_idx <- TRUE
  }

  roi <- handle$subset$roi_mask %||% handle$subset$roi
  partial_labels <- handle$meta$qpp_roi_dense_labels
  selected_parcels <- handle$meta$qpp_selected_parcels
  if (is.null(selected_parcels)) selected_parcels <- seq_len(ncol(parcel_series))
  selected_parcels <- as.integer(selected_parcels)
  if (length(selected_parcels) != ncol(parcel_series) || anyNA(selected_parcels) ||
      any(selected_parcels < 1L)) {
    abort_lna("decoded parcel columns have invalid parcel identities",
              .subclass = "lna_error_validation",
              location = "invert_step.parcel.mean:parcel_selection")
  }

  if (!is.null(partial_labels)) {
    if (is.null(roi)) {
      abort_lna("partial parcel labels require an roi_mask",
                .subclass = "lna_error_validation",
                location = "invert_step.parcel.mean:partial_labels")
    }
    dense_labels <- as.integer(partial_labels)
    if (anyNA(dense_labels) || any(dense_labels < 1L)) {
      abort_lna("partial parcel labels are invalid",
                .subclass = "lna_error_validation",
                location = "invert_step.parcel.mean:partial_labels")
    }
    handle$subset_applied$roi_mask <- TRUE
  } else {
    dense_labels <- assets$dense_labels %||% NULL
    if (is.null(dense_labels)) {
      entry <- .parcel_mean_dataset(desc, "parcel_labels")
      root <- handle$h5[["/"]]
      on.exit(safe_h5_close(root), add = TRUE)
      dense_labels <- as.integer(h5_read(root, entry$path))
    }
    if (length(dense_labels) == 0L || anyNA(dense_labels) ||
        any(dense_labels < 1L)) {
      abort_lna("parcel labels are incompatible with decoded parcel series",
                .subclass = "lna_error_validation",
                location = "invert_step.parcel.mean:dims")
    }
  }

  if (!is.null(roi) && is.null(partial_labels)) {
    roi <- as.logical(roi)
    file_mask <- handle$mask_info$mask %||% NULL
    if (!is.null(file_mask) && length(roi) == length(file_mask)) {
      roi <- roi[as.logical(file_mask)]
    }
    if (length(roi) != length(dense_labels)) {
      abort_lna("roi_mask does not match active voxel order",
                .subclass = "lna_error_validation",
                location = "invert_step.parcel.mean:roi_mask")
    }
    dense_labels <- dense_labels[roi]
    handle$subset_applied$roi_mask <- TRUE
  }
  gather <- match(dense_labels, selected_parcels)
  if (anyNA(gather)) {
    abort_lna("decoded parcel subset does not cover requested voxel labels",
              .subclass = "lna_error_validation",
              location = "invert_step.parcel.mean:gather")
  }
  dense <- parcel_series[, gather, drop = FALSE]
  dimnames(dense) <- NULL
  meta <- handle$meta
  meta$parcel_reconstruction <- list(
    matrix_order = "time-by-active-voxel",
    mask_checksum = p$mask_checksum,
    affine_checksum = p$affine_checksum,
    signal_units = p$signal_units,
    spatial_units = p$spatial_units
  )
  handle <- handle$with(meta = meta)
  handle$update_stash(parcel_key, setNames(list(dense), input_key))
}

.parcel_mean_validate_provenance <- function(p, location) {
  required <- c("labels", "fit_scope", "source_id", "signal_units", "spatial_units")
  # Inverse descriptors intentionally omit the large label array.
  if (!is.null(p$map_checksum)) required <- setdiff(required, "labels")
  missing <- required[vapply(required, function(nm) {
    is.null(p[[nm]]) || (is.character(p[[nm]]) && (!length(p[[nm]]) || !nzchar(p[[nm]][1L])))
  }, logical(1))]
  if (length(missing)) {
    abort_lna(paste("parcel.mean is missing required provenance:",
                    paste(missing, collapse = ", ")),
              .subclass = "lna_error_validation", location = location)
  }
  scopes <- c("fixed", "external-calibration", "fold", "run", "session", "subject")
  if (!is.character(p$fit_scope) || length(p$fit_scope) != 1L ||
      !p$fit_scope %in% scopes) {
    abort_lna("parcel.mean fit_scope is invalid",
              .subclass = "lna_error_validation", location = location)
  }
  if (identical(p$fit_scope, "external-calibration") &&
      (!is.character(p$calibration_id) || length(p$calibration_id) != 1L ||
       !nzchar(p$calibration_id) ||
       !is.character(p$calibration_checksum) ||
       length(p$calibration_checksum) != 1L ||
       !nzchar(p$calibration_checksum))) {
    abort_lna("external-calibration parcel maps require calibration_id and calibration_checksum",
              .subclass = "lna_error_validation", location = location)
  }
  if (identical(p$fit_scope, "external-calibration")) {
    .lna_canonical_sha256(
      p$calibration_checksum, "calibration_checksum", location = location
    )
  } else if (!is.null(p$calibration_id) || !is.null(p$calibration_checksum)) {
    abort_lna(
      "parcel calibration identity requires external-calibration fit_scope",
      .subclass = "lna_error_analysis_contract", location = location
    )
  }
  for (field in c(
    "map_mask_checksum", "map_grid_id", "map_affine_checksum",
    "mask_checksum", "affine_checksum", "grid_id", "map_checksum",
    "weight_checksum", "asset_id"
  )) {
    if (!is.null(p[[field]])) {
      .lna_canonical_sha256(p[[field]], field, location = location)
    }
  }
  .analysis_normalization_record(p$normalization_provenance)
  if (p$fit_scope %in% c("fold", "run", "session", "subject") &&
      (!is.list(p$fit_provenance) || !length(p$fit_provenance) ||
       is.null(names(p$fit_provenance)) || any(!nzchar(names(p$fit_provenance))))) {
    abort_lna("fitted parcel maps require named fit_provenance",
              .subclass = "lna_error_validation", location = location)
  }
  invisible(TRUE)
}

.qpp_assert_source_labels <- function(x, location) {
  if (!is.numeric(x) || any(!is.finite(x)) || any(x != trunc(x)) ||
      any(x < 1) || any(x > 2^53 - 1)) {
    abort_lna("active source parcel labels must be positive exact integers",
              .subclass = "lna_error_validation", location = location)
  }
  invisible(TRUE)
}

.parcel_mean_active_vector <- function(x, mask, supplied_checksum,
                                       current_checksum, what) {
  if (is.null(x)) {
    abort_lna(paste0("parcel.mean ", what, " are required"),
              .subclass = "lna_error_validation",
              location = paste0("forward_step.parcel.mean:", what))
  }
  if (is.array(x) && length(x) == length(mask)) {
    if (!identical(dim(x), dim(mask))) {
      abort_lna(paste0(what, " dimensions do not match the archive mask"),
                .subclass = "lna_error_validation",
                location = paste0("forward_step.parcel.mean:", what))
    }
    inactive <- as.vector(x)[!as.vector(mask)]
    if (identical(what, "labels") && any(!is.na(inactive) & inactive != 0)) {
      abort_lna("full-volume parcel labels must be zero outside the mask",
                .subclass = "lna_error_validation",
                location = "forward_step.parcel.mean:foreign_map")
    }
    value <- as.vector(x)[as.vector(mask)]
  } else if (is.atomic(x) && is.null(dim(x)) && length(x) == sum(mask)) {
    if (is.null(supplied_checksum) || !identical(supplied_checksum, current_checksum)) {
      abort_lna(
        paste0("active-order ", what, " require the matching map_mask_checksum"),
        .subclass = "lna_error_validation",
        location = "forward_step.parcel.mean:mask_order"
      )
    }
    value <- x
  } else {
    abort_lna(paste0(what, " must match full mask dimensions or active-mask length"),
              .subclass = "lna_error_validation",
              location = paste0("forward_step.parcel.mean:", what))
  }
  if (anyNA(value)) {
    abort_lna(paste0(what, " contain missing values inside the active mask"),
              .subclass = "lna_error_validation",
              location = paste0("forward_step.parcel.mean:", what))
  }
  as.numeric(value)
}

.parcel_mean_validate_affine <- function(p, header) {
  transform <- header$transform %||% NULL
  if (is.null(transform)) {
    abort_lna(
      "parcel map cannot be validated because input affine is missing",
      asset = p$source_id %||% "<unknown-parcel-map>",
      field = "input_transform",
      expected_identity = p$map_grid_id %||% p$map_affine_checksum,
      actual_identity = NA_character_,
      .subclass = "lna_error_spatial_mismatch",
      location = "forward_step.parcel.mean:affine"
    )
  }
  if (is.null(header$dims)) {
    abort_lna(
      "parcel map cannot be validated because input grid dimensions are missing",
      asset = p$source_id %||% "<unknown-parcel-map>",
      field = "input_dims",
      expected_identity = p$map_grid_id %||% NA_character_,
      actual_identity = NA_character_,
      .subclass = "lna_error_spatial_mismatch",
      location = "forward_step.parcel.mean:grid"
    )
  }
  identity <- .lna_grid_identity_from_header(
    header,
    units = p$spatial_units %||% header$spatial_units %||% "mm",
    voxel_ordering = p$voxel_ordering %||%
      header$voxel_ordering %||% "first-axis-fastest",
    location = "forward_step.parcel.mean:grid"
  )
  expected_grid <- p$map_grid_id %||% ""
  if (!nzchar(expected_grid)) {
    abort_lna(
      "parcel map is missing required grid identity",
      asset = p$source_id %||% "<unknown-parcel-map>",
      field = "map_grid_id",
      expected_identity = NA_character_,
      actual_identity = identity$grid_id,
      .subclass = "lna_error_spatial_mismatch",
      location = "forward_step.parcel.mean:grid"
    )
  }
  if (!identical(expected_grid, identity$grid_id)) {
    abort_lna(
      "parcel map grid identity does not match input geometry",
      asset = p$source_id %||% "<unknown-parcel-map>",
      field = "grid_id",
      expected_identity = expected_grid,
      actual_identity = identity$grid_id,
      .subclass = "lna_error_spatial_mismatch",
      location = "forward_step.parcel.mean:grid"
    )
  }
  if (!is.null(p$map_affine_checksum) &&
      !identical(p$map_affine_checksum, identity$affine_id)) {
    abort_lna(
      "parcel map affine checksum does not match input geometry",
      asset = p$source_id %||% "<unknown-parcel-map>",
      field = "affine_id",
      expected_identity = p$map_affine_checksum,
      actual_identity = identity$affine_id,
      .subclass = "lna_error_spatial_mismatch",
      location = "forward_step.parcel.mean:affine"
    )
  }
  identity
}

.parcel_mean_numeric_checksum <- function(x, prefix) {
  if (!is.numeric(x) || any(!is.finite(x))) {
    abort_lna("checksum input must be finite numeric data",
              .subclass = "lna_error_validation",
              location = ".parcel_mean_numeric_checksum")
  }
  .lna_numeric_id(x, prefix)
}

.parcel_mean_map_checksum <- function(dense_labels, source_ids, mask_checksum) {
  text <- paste("lna_parcel_map_v1", mask_checksum,
                paste(dense_labels, collapse = ","),
                paste(format(source_ids, scientific = FALSE, trim = TRUE), collapse = ","),
                sep = "|")
  paste0("sha256:", digest::digest(text, algo = "sha256", serialize = FALSE))
}

.parcel_mean_add_asset <- function(plan, path, value, dtype, role, run_id,
                                   step_index, params_json) {
  if (path %in% names(plan$payloads)) {
    if (!identical(plan$payloads[[path]], value)) {
      abort_lna(
        paste0("parcel asset identity collision at ", path),
        code = "parcel_asset_collision",
        .subclass = "lna_error_spatial_identity",
        location = ".parcel_mean_add_asset"
      )
    }
    return(invisible(plan))
  }
  plan$add_payload(path, value)
  plan$add_dataset_def(path, role, "parcel.mean", run_id, as.integer(step_index),
                       params_json, path, "eager", dtype = dtype)
  invisible(plan)
}

.parcel_mean_dataset <- function(desc, role) {
  entries <- Filter(function(x) identical(x$role, role), desc$datasets %||% list())
  if (length(entries) != 1L) {
    abort_lna(paste0("parcel.mean descriptor requires one ", role, " dataset"),
              .subclass = "lna_error_descriptor",
              location = "invert_step.parcel.mean:descriptor")
  }
  entries[[1L]]
}

#' Default parcel.mean parameters
#' @return A list of defaults. Required scientific provenance has no default.
#' @keywords internal
lna_default.parcel.mean <- function() {
  list(weights = NULL, source_labels = NULL, map_mask_checksum = NULL,
       map_grid_id = NULL, map_affine_checksum = NULL,
       voxel_ordering = "first-axis-fastest", calibration_id = NULL,
       calibration_checksum = NULL, fit_provenance = NULL,
       normalization_provenance = NULL)
}
