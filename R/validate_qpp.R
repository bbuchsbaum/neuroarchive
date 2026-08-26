# QPP archive validation and deterministic validation receipts.

.qpp_validation_default_limits <- list(
  max_time = 10000000,
  max_parcels = 65535,
  max_tiles = 1000000,
  max_payload_bytes = 1024^4,
  max_index_bytes = 1024^3,
  max_map_values = 10000000,
  max_tile_values = .qpp_default_max_values,
  max_tile_bytes = .qpp_default_max_bytes,
  max_uncompressed_tile_bytes = .qpp_default_max_bytes
)

.qpp_validation_limits <- function(overrides = NULL) {
  if (is.null(overrides)) return(.qpp_validation_default_limits)
  if (!is.list(overrides) || is.null(names(overrides)) ||
      any(!nzchar(names(overrides)))) {
    abort_lna("qpp_limits must be a named list",
              .subclass = "lna_error_validation",
              location = "validate_lna:qpp_limits")
  }
  unknown <- setdiff(names(overrides), names(.qpp_validation_default_limits))
  if (length(unknown)) {
    abort_lna(paste0("unknown qpp_limits: ", paste(unknown, collapse = ", ")),
              .subclass = "lna_error_validation",
              location = "validate_lna:qpp_limits")
  }
  out <- utils::modifyList(.qpp_validation_default_limits, overrides)
  valid <- vapply(out, function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 &&
      x == trunc(x)
  }, logical(1))
  if (!all(valid)) {
    abort_lna("all qpp validation limits must be positive exact integers",
              .subclass = "lna_error_validation",
              location = "validate_lna:qpp_limits")
  }
  out
}

.qpp_validation_abort <- function(message, location, code = "qpp_corrupt",
                                  subclass = "lna_error_qpp_corrupt") {
  if (subclass %in% c("lna_error_qpp_integrity", "lna_error_qpp_binding")) {
    subclass <- c(subclass, "lna_error_qpp_corrupt")
  }
  abort_lna(message, .subclass = subclass, location = location, code = code)
}

.qpp_validation_exact_integer <- function(x, name, min = 1, max = Inf) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != trunc(x) || x < min || x > max) {
    .qpp_validation_abort(
      paste0("qpp ", name, " exceeds its resource cap or is not an exact integer"),
      paste0("validate_lna:qpp:", name), "qpp_resource_limit"
    )
  }
  as.numeric(x)
}

.qpp_validation_entry <- function(desc, role) {
  entries <- Filter(function(x) identical(x$role, role), desc$datasets %||% list())
  if (length(entries) != 1L || !is.character(entries[[1L]]$path) ||
      length(entries[[1L]]$path) != 1L || !nzchar(entries[[1L]]$path)) {
    .qpp_validation_abort(
      paste0("qpp descriptor requires exactly one ", role, " dataset"),
      "validate_lna:qpp:descriptor", "qpp_descriptor"
    )
  }
  entries[[1L]]
}

.qpp_validation_extent <- function(root, entry, role) {
  value <- tryCatch(
    .qpp_dataset_length(root, entry$path),
    error = function(e) {
      .qpp_validation_abort(
        conditionMessage(e), paste0("validate_lna:qpp:", role),
        "qpp_dataset_extent"
      )
    }
  )
  as.numeric(value)
}

.qpp_validation_dims <- function(root, entry, role) {
  tryCatch(
    .qpp_dataset_dims(root, entry$path),
    error = function(error) {
      .qpp_validation_abort(
        conditionMessage(error), paste0("validate_lna:qpp:", role),
        "qpp_dataset_extent"
      )
    }
  )
}

.qpp_read_quality_receipt <- function(root, entry) {
  text <- tryCatch(
    h5_read(root, entry$path),
    error = function(e) {
      .qpp_validation_abort(
        paste0("cannot read qpp quality receipt: ", conditionMessage(e)),
        "validate_lna:qpp:quality_receipt", "qpp_receipt"
      )
    }
  )
  if (!is.character(text) || length(text) != 1L || !jsonlite::validate(text)) {
    .qpp_validation_abort("qpp quality receipt is not one valid JSON document",
                          "validate_lna:qpp:quality_receipt", "qpp_receipt")
  }
  jsonlite::fromJSON(
    text, simplifyVector = TRUE, simplifyDataFrame = FALSE,
    simplifyMatrix = FALSE
  )
}

.qpp_validate_receipt_metadata <- function(receipt, desc, index, steps,
                                           payload_length) {
  required <- c(
    "receipt_version", "codec_version", "decoder_version", "tile_count",
    "logical_values", "stored_payload_bytes", "maximum_absolute_error",
    "parcelwise_RMSE", "parcel_steps_checksum", "payload_checksum",
    "index_checksum", "transform_settings_checksum", "q_roundtrip_evidence",
    "tile_payload_checksums", "software_versions", "receipt_checksum"
  )
  missing <- required[vapply(required, function(nm) is.null(receipt[[nm]]), logical(1))]
  if (length(missing)) {
    .qpp_validation_abort(
      paste0("qpp quality receipt is missing: ", paste(missing, collapse = ", ")),
      "validate_lna:qpp:quality_receipt", "qpp_receipt"
    )
  }
  if (!identical(receipt$receipt_version, "1.0.0") ||
      !identical(receipt$codec_version, "1.0.0") ||
      !identical(receipt$decoder_version, "neuroarchive-qpp-cpp-1")) {
    .qpp_validation_abort("unsupported qpp quality receipt or decoder version",
                          "validate_lna:qpp:quality_receipt", "qpp_version")
  }
  stored_receipt_checksum <- receipt$receipt_checksum
  receipt$receipt_checksum <- NULL
  actual_receipt_checksum <- .qpp_canonical_checksum(
    receipt, "lna_qpp_quality_receipt_v1"
  )
  if (!identical(stored_receipt_checksum, actual_receipt_checksum)) {
    .qpp_validation_abort("qpp quality receipt checksum mismatch",
                          "validate_lna:qpp:quality_receipt_checksum",
                          "qpp_receipt_checksum")
  }
  p <- desc$params
  expected_values <- as.numeric(p$time_count) * as.numeric(p$parcel_count)
  if (as.numeric(receipt$tile_count) != nrow(index) ||
      as.numeric(receipt$logical_values) != expected_values ||
      as.numeric(receipt$stored_payload_bytes) != payload_length ||
      length(receipt$tile_payload_checksums) != nrow(index) ||
      length(receipt$q_roundtrip_evidence$tile_q_checksums) != nrow(index) ||
      length(receipt$parcelwise_RMSE) != as.numeric(p$parcel_count)) {
    .qpp_validation_abort("qpp quality receipt dimensions contradict the archive",
                          "validate_lna:qpp:quality_receipt_shape",
                          "qpp_receipt")
  }
  if (!is.numeric(receipt$maximum_absolute_error) ||
      length(receipt$maximum_absolute_error) != 1L ||
      !is.finite(receipt$maximum_absolute_error) ||
      receipt$maximum_absolute_error < 0 ||
      any(!is.finite(receipt$parcelwise_RMSE)) ||
      any(receipt$parcelwise_RMSE < 0)) {
    .qpp_validation_abort("qpp quality metrics are invalid",
                          "validate_lna:qpp:quality_metrics", "qpp_receipt")
  }
  if (!is.null(steps)) {
    tolerance <- 64 * .Machine$double.eps * max(1, max(steps))
    if (receipt$maximum_absolute_error > max(steps) / 2 + tolerance ||
        any(receipt$parcelwise_RMSE > max(steps) / 2 + tolerance)) {
      .qpp_validation_abort("qpp quality metrics exceed the quantizer error bound",
                            "validate_lna:qpp:quality_metrics", "qpp_receipt")
    }
  }
  if (!identical(
    receipt$transform_settings_checksum,
    .qpp_canonical_checksum(p, "lna_qpp_transform_settings_v1")
  )) {
    .qpp_validation_abort("qpp transform-settings checksum mismatch",
                          "validate_lna:qpp:settings_checksum",
                          "qpp_metadata_checksum")
  }
  if (!identical(
    receipt$index_checksum,
    .qpp_canonical_checksum(index, "lna_qpp_index_v1")
  )) {
    .qpp_validation_abort("qpp index checksum mismatch",
                          "validate_lna:qpp:index_checksum",
                          "qpp_metadata_checksum")
  }
  if (!is.null(steps) && !identical(
    receipt$parcel_steps_checksum,
    .parcel_mean_numeric_checksum(steps, "lna_qpp_steps_v1")
  )) {
    .qpp_validation_abort("qpp parcel-step checksum mismatch",
                          "validate_lna:qpp:steps_checksum",
                          "qpp_metadata_checksum")
  }
  invisible(TRUE)
}

.qpp_validate_index <- function(index, rows, p, payload_length, limits,
                                full_index = FALSE) {
  if (!is.data.frame(index) || nrow(index) != length(rows)) {
    .qpp_validation_abort("qpp index read did not return the requested rows",
                          "validate_lna:qpp:index", "qpp_index")
  }
  uncompressed <- as.numeric(index$uncompressed_length)
  if (any(!is.finite(uncompressed)) || any(uncompressed < 1) ||
      any(uncompressed > limits$max_uncompressed_tile_bytes)) {
    .qpp_validation_abort("qpp declared uncompressed tile length exceeds its cap",
                          "validate_lna:qpp:uncompressed_length",
                          "qpp_decompression_bomb")
  }
  if (any(as.numeric(index$outer_codec) != 0)) {
    .qpp_validation_abort("qpp index uses an unsupported outer codec",
                          "validate_lna:qpp:outer_codec", "qpp_outer_codec")
  }
  time_chunk <- as.integer(p$time_chunk)
  parcel_chunk <- as.integer(p$parcel_chunk)
  n_parcel_tiles <- ceiling(as.numeric(p$parcel_count) / parcel_chunk)
  for (i in seq_along(rows)) {
    tryCatch(
      .qpp_validate_index_row(
        index[i, , drop = FALSE], rows[[i]], as.integer(p$time_count),
        as.integer(p$parcel_count), time_chunk, parcel_chunk,
        n_parcel_tiles, payload_length
      ),
      error = function(e) {
        .qpp_validation_abort(conditionMessage(e),
                              "validate_lna:qpp:index_row", "qpp_index")
      }
    )
    if (as.numeric(index$flags[[i]]) != 1 ||
        as.numeric(index$digest_row[[i]]) != rows[[i]] - 1L) {
      .qpp_validation_abort(
        "qpp index digest reference is not canonical",
        "validate_lna:qpp:index_digest", "qpp_digest_index",
        "lna_error_qpp_integrity"
      )
    }
  }
  tryCatch(
    .qpp_validate_selected_offsets(
      index, rows,
      if (isTRUE(full_index)) length(rows) else max(rows),
      payload_length
    ),
    error = function(e) {
      .qpp_validation_abort(conditionMessage(e),
                            "validate_lna:qpp:index_ranges", "qpp_index")
    }
  )
  invisible(TRUE)
}

.qpp_validation_sample_rows <- function(tile_count, sample_tiles) {
  count <- min(tile_count, as.integer(sample_tiles))
  sort(unique(as.integer(round(seq(1, tile_count, length.out = count)))))
}

.qpp_validate_tiles <- function(root, payload_entry, digest_entry, index, rows,
                                receipt, p, mode) {
  q_checksums <- character(length(rows))
  payload_checksums <- character(length(rows))
  checked_bytes <- 0
  stored_digests <- .qpp_read_digest_rows(
    root, digest_entry$path, rows, .qpp_new_io_stats(), "tile_digests"
  )
  for (i in seq_along(rows)) {
    row <- index[rows[[i]], , drop = FALSE]
    offset <- as.numeric(row$byte_offset)
    byte_length <- as.numeric(row$byte_length)
    bytes <- .qpp_payload_as_raw(h5_read_subset(
      root, payload_entry$path,
      list(seq.int(offset + 1, offset + byte_length))
    ))
    checked_bytes <- checked_bytes + byte_length
    if (!identical(
      .qpp_sha256_raw(bytes), as.raw(stored_digests[i, , drop = TRUE])
    )) {
      .qpp_validation_abort(
        "qpp tile stored-byte digest mismatch",
        "validate_lna:qpp:tile_digest", "qpp_tile_digest_mismatch",
        "lna_error_qpp_integrity"
      )
    }
    tile <- tryCatch(
      qpp_decode_tile(
        bytes,
        max_time = as.integer(row$time_count),
        max_parcels = as.integer(row$parcel_count),
        max_values = p$max_tile_values %||% .qpp_default_max_values,
        max_bytes = p$max_tile_bytes %||% .qpp_default_max_bytes
      ),
      error = function(e) {
        .qpp_validation_abort(
          paste0("qpp tile decode failed: ", conditionMessage(e)),
          "validate_lna:qpp:tile_decode", "qpp_tile_corrupt"
        )
      }
    )
    payload_checksum <- .qpp_raw_checksum(bytes, "lna_qpp_tile_payload_v1")
    if (!identical(payload_checksum, receipt$tile_payload_checksums[[rows[[i]]]])) {
      .qpp_validation_abort("qpp tile payload checksum mismatch",
                            "validate_lna:qpp:tile_payload_checksum",
                            "qpp_payload_checksum")
    }
    q_checksum <- .parcel_mean_numeric_checksum(tile, "lna_qpp_tile_q_v1")
    if (!identical(q_checksum,
                   receipt$q_roundtrip_evidence$tile_q_checksums[[rows[[i]]]])) {
      .qpp_validation_abort("qpp decoded-q checksum mismatch",
                            "validate_lna:qpp:q_checksum",
                            "qpp_q_roundtrip")
    }
    reencoded <- qpp_encode_tile(
      tile,
      max_values = p$max_tile_values %||% .qpp_default_max_values,
      max_bytes = p$max_tile_bytes %||% .qpp_default_max_bytes
    )
    if (!identical(reencoded, bytes)) {
      .qpp_validation_abort("qpp decode/re-encode roundtrip changed tile bytes",
                            "validate_lna:qpp:q_roundtrip",
                            "qpp_q_roundtrip")
    }
    q_checksums[[i]] <- q_checksum
    payload_checksums[[i]] <- payload_checksum
  }
  if (identical(mode, "all")) {
    if (!identical(
      .qpp_canonical_checksum(payload_checksums, "lna_qpp_payload_v1"),
      receipt$payload_checksum
    )) {
      .qpp_validation_abort("qpp complete payload checksum mismatch",
                            "validate_lna:qpp:payload_checksum",
                            "qpp_payload_checksum")
    }
    if (!identical(
      .qpp_canonical_checksum(q_checksums, "lna_qpp_all_tile_q_v1"),
      receipt$q_roundtrip_evidence$combined_checksum
    )) {
      .qpp_validation_abort("qpp complete decoded-q evidence mismatch",
                            "validate_lna:qpp:q_evidence",
                            "qpp_q_roundtrip")
    }
  }
  list(
    tiles_checked = length(rows),
    payload_bytes_checked = checked_bytes,
    q_roundtrip_verified = TRUE
  )
}

.qpp_validate_parcel_map <- function(root, descriptors, qpp_desc, limits,
                                     quick = FALSE) {
  if (!identical(qpp_desc$params$profile_id, "analysis-safe")) {
    return(invisible(TRUE))
  }
  binding <- qpp_desc$params$parcel_binding %||% NULL
  tryCatch(
    .qpp_validate_parcel_binding_fields(
      binding, run_id = qpp_desc$run_id,
      location = "validate_lna:qpp:parcel_binding"
    ),
    error = function(error) {
      .qpp_validation_abort(
        conditionMessage(error),
        error$location %||% "validate_lna:qpp:parcel_binding",
        error$code %||% "qpp_parcel_binding",
        "lna_error_qpp_binding"
      )
    }
  )
  semantic_matches <- vapply(descriptors, function(parcel) {
    p <- parcel$params %||% list()
    labels <- Filter(
      function(dataset) identical(dataset$role, "parcel_labels"),
      parcel$datasets %||% list()
    )
    identical(parcel$type, "parcel.mean") &&
      identical(parcel$version, "1.1.0") &&
      identical(parcel$run_id, binding$run_id) &&
      identical(as.integer(p$step_index), as.integer(binding$step_index)) &&
      identical(p$asset_id, binding$asset_id) &&
      identical(p$map_checksum, binding$map_checksum) &&
      identical(p$mask_checksum, binding$mask_checksum) &&
      length(labels) == 1L && identical(labels[[1L]]$path, binding$labels_path)
  }, logical(1))
  match_names <- names(descriptors)[semantic_matches]
  if (length(match_names) != 1L || !identical(match_names, binding$descriptor)) {
    code <- if (!length(match_names)) "qpp_parcel_missing" else "qpp_parcel_duplicate"
    .qpp_validation_abort(
      "analysis-safe qpp binding does not resolve exactly one parcel descriptor",
      "validate_lna:qpp:parcel_binding", code, "lna_error_qpp_binding"
    )
  }
  parcel_desc <- descriptors[[binding$descriptor]]
  label_entry <- .qpp_validation_entry(parcel_desc, "parcel_labels")
  source_entry <- .qpp_validation_entry(parcel_desc, "source_parcel_ids")
  label_count <- .qpp_validation_extent(root, label_entry, "parcel_labels")
  source_count <- .qpp_validation_extent(root, source_entry, "source_parcel_ids")
  if (label_count > limits$max_map_values ||
      source_count != as.numeric(qpp_desc$params$parcel_count)) {
    .qpp_validation_abort("parcel map extent contradicts qpp parcel_count or cap",
                          "validate_lna:qpp:parcel_map_extent",
                          "qpp_parcel_map")
  }
  if (isTRUE(quick)) return(invisible(TRUE))
  labels <- as.integer(h5_read(root, label_entry$path))
  source_ids <- as.numeric(h5_read(root, source_entry$path))
  k <- as.integer(qpp_desc$params$parcel_count)
  if (length(labels) != label_count || anyNA(labels) || any(labels < 1L) ||
      any(labels > k) || !identical(sort(unique(labels)), seq_len(k))) {
    .qpp_validation_abort("parcel map labels do not exactly cover qpp parcels",
                          "validate_lna:qpp:parcel_map_labels",
                          "qpp_parcel_map")
  }
  expected <- .parcel_mean_map_checksum(
    labels, source_ids, parcel_desc$params$mask_checksum
  )
  if (!identical(expected, parcel_desc$params$map_checksum)) {
    .qpp_validation_abort("parcel map checksum mismatch",
                          "validate_lna:qpp:parcel_map_checksum",
                          "qpp_parcel_map")
  }
  invisible(TRUE)
}

.validate_one_qpp <- function(root, desc, descriptors, name, mode,
                              sample_tiles, limits, quick = FALSE) {
  if (!identical(desc$version, "1.1.0") ||
      !identical(desc$codec_version, "1.0.0") ||
      !identical(desc$decoder_version, "neuroarchive-qpp-cpp-1")) {
    .qpp_validation_abort("unsupported required qpp decoder version",
                          "validate_lna:qpp:decoder_version",
                          "required_transform",
                          "lna_error_required_transform")
  }
  p <- desc$params %||% list()
  if (identical(p$profile_id, "analysis-safe")) {
    required <- if (h5_attr_exists(root, "required_transforms")) {
      as.character(h5_attr_read(root, "required_transforms"))
    } else {
      character()
    }
    expected_required <- c("parcel.mean@1.1.0", "qpp@1.1.0")
    if (!all(expected_required %in% required)) {
      .qpp_validation_abort(
        "analysis-safe qpp requires exact parcel.mean and qpp decoder declarations",
        "validate_lna:qpp:required_transforms", "required_transform",
        "lna_error_required_transform"
      )
    }
  }
  tryCatch({
    expected_contract <- .qpp_analysis_contract(p)
    .analysis_assert_contract_match(
      desc$analysis_contract, expected_contract, "qpp@1.1.0"
    )
  }, error = function(e) {
    .qpp_validation_abort(
      paste0("invalid qpp calibration or analysis contract: ", conditionMessage(e)),
      "validate_lna:qpp:calibration_contract", "qpp_calibration"
    )
  })
  time_count <- .qpp_validation_exact_integer(
    p$time_count, "time_count", max = limits$max_time
  )
  parcel_count <- .qpp_validation_exact_integer(
    p$parcel_count, "parcel_count", max = limits$max_parcels
  )
  time_chunk <- .qpp_validation_exact_integer(
    p$time_chunk, "time_chunk", max = 65535
  )
  parcel_chunk <- .qpp_validation_exact_integer(
    p$parcel_chunk, "parcel_chunk", max = 65535
  )
  tile_count <- ceiling(time_count / time_chunk) *
    ceiling(parcel_count / parcel_chunk)
  tile_count <- as.integer(.qpp_validation_exact_integer(
    tile_count, "tile_count", max = limits$max_tiles
  ))
  if (tile_count * 72 > limits$max_index_bytes) {
    .qpp_validation_abort("qpp index exceeds validation byte cap",
                          "validate_lna:qpp:index_cap", "qpp_resource_limit")
  }
  payload_entry <- .qpp_validation_entry(desc, "qpp_payload")
  index_entry <- .qpp_validation_entry(desc, "qpp_tile_index")
  step_entry <- .qpp_validation_entry(desc, "qpp_parcel_steps")
  receipt_entry <- .qpp_validation_entry(desc, "qpp_quality_receipt")
  tile_digest_entry <- .qpp_validation_entry(desc, "qpp_tile_sha256")
  step_digest_entry <- .qpp_validation_entry(desc, "qpp_parcel_step_sha256")
  payload_length <- .qpp_validation_extent(root, payload_entry, "payload")
  index_length <- .qpp_validation_extent(root, index_entry, "index")
  step_length <- .qpp_validation_extent(root, step_entry, "steps")
  tile_digest_dims <- .qpp_validation_dims(
    root, tile_digest_entry, "tile_digests"
  )
  step_digest_dims <- .qpp_validation_dims(
    root, step_digest_entry, "step_digests"
  )
  if (payload_length < 1 || payload_length > limits$max_payload_bytes ||
      index_length != tile_count || step_length != parcel_count ||
      !identical(tile_digest_dims, c(as.numeric(tile_count), 32)) ||
      !identical(step_digest_dims, c(as.numeric(parcel_count), 32))) {
    .qpp_validation_abort("qpp dataset extents contradict dimensions or resource caps",
                          "validate_lna:qpp:dataset_extents",
                          "qpp_resource_limit")
  }

  rows <- if (isTRUE(quick)) {
    sort(unique(c(1L, as.integer(ceiling(tile_count / 2)), as.integer(tile_count))))
  } else {
    seq_len(tile_count)
  }
  index <- h5_read_subset(root, index_entry$path, list(rows))
  .qpp_validate_index(index, rows, p, payload_length, limits,
                      full_index = !isTRUE(quick))
  .qpp_validate_parcel_map(root, descriptors, desc, limits, quick = quick)
  if (isTRUE(quick)) {
    return(list(
      transform = paste0(desc$type, "@", desc$version), descriptor = name,
      validation_mode = "quick-structural", tiles_checked = length(rows),
      tile_count = tile_count, payload_bytes_checked = 0,
      q_roundtrip_verified = FALSE
    ))
  }

  steps <- as.numeric(h5_read(root, step_entry$path))
  if (length(steps) != parcel_count || any(!is.finite(steps)) || any(steps <= 0)) {
    .qpp_validation_abort("qpp parcel steps are invalid",
                          "validate_lna:qpp:steps", "qpp_steps")
  }
  step_digests <- .qpp_read_digest_rows(
    root, step_digest_entry$path, seq_len(parcel_count),
    .qpp_new_io_stats(), "step_digests"
  )
  for (parcel in seq_len(parcel_count)) {
    if (!identical(
      .qpp_step_digest_raw(steps[[parcel]], parcel - 1L),
      as.raw(step_digests[parcel, , drop = TRUE])
    )) {
      .qpp_validation_abort(
        "qpp parcel-step digest mismatch",
        "validate_lna:qpp:step_digest", "qpp_step_digest_mismatch",
        "lna_error_qpp_integrity"
      )
    }
  }
  receipt <- .qpp_read_quality_receipt(root, receipt_entry)
  .qpp_validate_receipt_metadata(receipt, desc, index, steps, payload_length)
  if (identical(mode, "structural")) {
    decoded <- list(tiles_checked = 0L, payload_bytes_checked = 0,
                    q_roundtrip_verified = FALSE)
  } else {
    decode_rows <- if (identical(mode, "all")) {
      seq_len(tile_count)
    } else {
      .qpp_validation_sample_rows(tile_count, sample_tiles)
    }
    decoded <- .qpp_validate_tiles(
      root, payload_entry, tile_digest_entry, index, decode_rows,
      receipt, p, mode
    )
  }
  list(
    transform = paste0(desc$type, "@", desc$version), descriptor = name,
    validation_mode = mode, tile_count = tile_count,
    tiles_checked = decoded$tiles_checked,
    payload_bytes_checked = decoded$payload_bytes_checked,
    payload_bytes = payload_length,
    q_roundtrip_verified = decoded$q_roundtrip_verified,
    maximum_absolute_error = receipt$maximum_absolute_error,
    parcelwise_RMSE = receipt$parcelwise_RMSE,
    checksums = list(
      receipt = receipt$receipt_checksum,
      payload = receipt$payload_checksum,
      index = receipt$index_checksum,
      steps = receipt$parcel_steps_checksum,
      settings = receipt$transform_settings_checksum
    ),
    software_versions = receipt$software_versions
  )
}

.validate_qpp_archives <- function(h5, root, fail, mode = "all",
                                   sample_tiles = 16L, limits = NULL,
                                   quick = FALSE) {
  limits <- .qpp_validation_limits(limits)
  if (!h5$exists("transforms")) return(list())
  group <- h5[["transforms"]]
  on.exit(safe_h5_close(group), add = TRUE)
  discovered <- discover_transforms(group)
  descriptors <- lapply(discovered$name, function(name) {
    read_json_descriptor(group, name)
  })
  names(descriptors) <- discovered$name
  qpp_names <- discovered$name[discovered$type == "qpp"]
  receipts <- list()
  for (name in qpp_names) {
    result <- tryCatch(
      .validate_one_qpp(
        root, descriptors[[name]], descriptors, name, mode,
        sample_tiles, limits, quick
      ),
      error = function(e) e
    )
    if (inherits(result, "error")) {
      subclass <- class(result)[grepl("^lna_error_", class(result))][1L]
      if (length(subclass) != 1L || is.na(subclass) || !nzchar(subclass)) {
        subclass <- "lna_error_qpp_corrupt"
      }
      fail(
        conditionMessage(result),
        code = result$code %||% "qpp_corrupt",
        path = name,
        layer = "qpp",
        subclass = subclass,
        location = result$location %||% "validate_lna:qpp"
      )
    } else {
      receipts[[name]] <- result
    }
  }
  receipts
}

.validate_required_transforms <- function(h5, root, fail) {
  if (!h5_attr_exists(root, "required_transforms") || !h5$exists("transforms")) {
    return(invisible(TRUE))
  }
  required <- as.character(h5_attr_read(root, "required_transforms"))
  required <- required[nzchar(required)]
  if (!length(required)) return(invisible(TRUE))
  group <- h5[["transforms"]]
  on.exit(safe_h5_close(group), add = TRUE)
  discovered <- discover_transforms(group)
  descriptors <- lapply(discovered$name, function(name) {
    read_json_descriptor(group, name)
  })
  available <- vapply(descriptors, function(desc) {
    paste0(desc$type %||% "", "@", desc$version %||% "")
  }, character(1))
  for (id in required) {
    if (!grepl("^[A-Za-z][A-Za-z0-9_.]*@[0-9]+\\.[0-9]+\\.[0-9]+$", id)) {
      fail(paste0("invalid required transform id: ", id),
           code = "required_transform", layer = "decoder",
           subclass = "lna_error_required_transform",
           location = "validate_lna:required_transforms")
      next
    }
    match_idx <- which(available == id)
    if (!length(match_idx)) {
      fail(paste0("required decoder is unavailable: ", id),
           code = "required_transform", layer = "decoder",
           subclass = "lna_error_required_transform",
           location = "validate_lna:required_transforms")
      next
    }
    desc <- descriptors[[match_idx[[1L]]]]
    built_in <- id %in% c(
      "parcel.mean@1.0.0", "qpp@1.0.0",
      "parcel.mean@1.1.0", "qpp@1.1.0"
    )
    method <- .resolve_transform_method("invert_step", desc$type)
    if (!built_in && isTRUE(desc$required_decoder)) {
      package <- desc$decoder_package %||% NULL
      version <- desc$decoder_package_version %||% NULL
      installed <- is.character(package) && length(package) == 1L &&
        requireNamespace(package, quietly = TRUE)
      exact_version <- installed && is.character(version) && length(version) == 1L &&
        identical(as.character(utils::packageVersion(package)), version)
      if (!exact_version) method <- NULL
    }
    if (is.null(method)) {
      fail(paste0("no exact installed decoder for required transform: ", id),
           code = "required_transform", layer = "decoder",
           subclass = "lna_error_required_transform",
           location = "validate_lna:required_transforms")
    }
  }
  invisible(TRUE)
}
