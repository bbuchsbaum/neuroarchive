# Fused pointwise-quantize, predict, and pack transform.

#' Forward step for qpp
#'
#' @param type Transform type string.
#' @param desc Transform descriptor.
#' @param handle DataHandle for one run.
#' @return Updated DataHandle whose terminal representation is registered in the Plan.
#' @method forward_step qpp
#' @export
#' @keywords internal
forward_step.qpp <- function(type, desc, handle) {
  p <- desc$params %||% list()
  input_key <- desc$inputs[[1]] %||% "parcel_series"
  input <- handle$get_inputs(input_key)[[1]]
  if (is.list(input) && !is.matrix(input)) {
    abort_lna("qpp accepts one run/independence unit at a time",
              .subclass = "lna_error_validation",
              location = "forward_step.qpp:runs")
  }
  y <- as.matrix(input)
  if (!is.numeric(y) || !nrow(y) || !ncol(y) || any(!is.finite(y))) {
    abort_lna("qpp input must be a non-empty finite time-by-parcel matrix",
              .subclass = "lna_error_validation",
              location = "forward_step.qpp:input")
  }
  policy <- p$quantizer_policy %||% "fixed"
  allowed <- c("fixed", "external-calibration", "fitted-run-maxabs")
  if (!is.character(policy) || length(policy) != 1L || !policy %in% allowed) {
    abort_lna(paste0("unsupported qpp quantizer_policy; use one of ",
                     paste(allowed, collapse = ", ")),
              .subclass = "lna_error_validation",
              location = "forward_step.qpp:quantizer_policy")
  }
  units <- p$units
  if (!is.character(units) || length(units) != 1L || !nzchar(units)) {
    abort_lna("qpp units are required",
              .subclass = "lna_error_validation",
              location = "forward_step.qpp:units")
  }
  parcel_assets <- handle$meta$parcel_mean_assets %||% NULL
  parcel_units <- parcel_assets$signal_units %||% NULL
  if (!is.null(parcel_units) && !identical(units, parcel_units)) {
    abort_lna("qpp units do not match parcel.mean signal units",
              .subclass = "lna_error_validation",
              location = "forward_step.qpp:units")
  }
  if (isTRUE(p$allow_clip %||% FALSE)) {
    abort_lna("qpp never permits clipping; change delta or use a different codec",
              .subclass = "lna_error_validation",
              location = "forward_step.qpp:clipping")
  }

  safe_policy <- policy %in% c("fixed", "external-calibration")
  has_parcel_binding <- is.list(parcel_assets) && all(vapply(
    c(
      "descriptor", "transform", "run_id", "step_index", "asset_id",
      "map_checksum", "mask_checksum", "labels_path"
    ),
    function(field) !is.null(parcel_assets[[field]]), logical(1)
  ))
  analysis_safe <- safe_policy && has_parcel_binding
  if (safe_policy) {
    delta <- .qpp_validate_delta(p$delta, ncol(y))
    fit_scope <- if (identical(policy, "fixed")) "fixed" else "external-calibration"
    if (identical(policy, "external-calibration")) {
      for (field in c("calibration_id", "calibration_checksum")) {
        if (!is.character(p[[field]]) || length(p[[field]]) != 1L || !nzchar(p[[field]])) {
          abort_lna(paste0("external-calibration qpp requires ", field),
                    .subclass = "lna_error_validation",
                    location = "forward_step.qpp:calibration")
        }
      }
      .lna_canonical_sha256(
        p$calibration_checksum, "calibration_checksum",
        location = "forward_step.qpp:calibration_checksum"
      )
    } else if (!is.null(p$calibration_id) || !is.null(p$calibration_checksum)) {
      abort_lna(
        "qpp calibration identity requires external-calibration policy",
        .subclass = "lna_error_analysis_contract",
        location = "forward_step.qpp:calibration"
      )
    }
  } else {
    if (!is.null(p$delta)) {
      abort_lna("fitted-run-maxabs derives delta and forbids an explicit delta",
                .subclass = "lna_error_validation",
                location = "forward_step.qpp:fitted_delta")
    }
    if (!is.character(p$fit_id) || length(p$fit_id) != 1L || !nzchar(p$fit_id)) {
      abort_lna("fitted-run-maxabs requires fit_id provenance",
                .subclass = "lna_error_validation",
                location = "forward_step.qpp:fit_id")
    }
    target_q <- as.numeric(p$target_q %||% 32767)
    zero_step <- as.numeric(p$zero_step %||% 1)
    if (length(target_q) != 1L || !is.finite(target_q) ||
        target_q != trunc(target_q) || target_q < 1 || target_q > 2^31 - 1 ||
        length(zero_step) != 1L || !is.finite(zero_step) || zero_step <= 0) {
      abort_lna("fitted-run-maxabs target_q/zero_step are invalid",
                .subclass = "lna_error_validation",
                location = "forward_step.qpp:fitted_params")
    }
    maxima <- apply(abs(y), 2L, max)
    delta <- ifelse(maxima == 0, zero_step, maxima / target_q)
    fit_scope <- "run"
  }

  time_chunk <- .qpp_chunk_size(p$time_chunk %||% 128L, "time_chunk")
  parcel_chunk <- .qpp_chunk_size(p$parcel_chunk %||% 128L, "parcel_chunk")
  max_tile_values <- as.numeric(p$max_tile_values %||% .qpp_default_max_values)
  max_tile_bytes <- as.numeric(p$max_tile_bytes %||% .qpp_default_max_bytes)
  time_starts <- seq.int(1L, nrow(y), by = time_chunk)
  parcel_starts <- seq.int(1L, ncol(y), by = parcel_chunk)

  connection <- rawConnection(raw(), open = "wb")
  on.exit(close(connection), add = TRUE)
  index_rows <- vector("list", length(time_starts) * length(parcel_starts))
  tile_id <- 0L
  byte_offset <- 0
  maximum_error <- 0
  logical_values <- 0
  row_id <- 0L
  parcel_squared_error <- numeric(ncol(y))
  parcel_observations <- numeric(ncol(y))
  tile_q_checksums <- character(length(index_rows))
  tile_payload_checksums <- character(length(index_rows))
  tile_payload_digests <- vector("list", length(index_rows))
  for (time_start in time_starts) {
    time_end <- min(nrow(y), time_start + time_chunk - 1L)
    for (parcel_start in parcel_starts) {
      parcel_end <- min(ncol(y), parcel_start + parcel_chunk - 1L)
      tile_y <- y[time_start:time_end, parcel_start:parcel_end, drop = FALSE]
      tile_delta <- delta[parcel_start:parcel_end]
      q <- .qpp_quantize_pointwise(tile_y, tile_delta)
      decoded <- sweep(q, 2L, tile_delta, "*")
      error <- decoded - tile_y
      maximum_error <- max(maximum_error, max(abs(error)))
      parcel_squared_error[parcel_start:parcel_end] <-
        parcel_squared_error[parcel_start:parcel_end] + colSums(error^2)
      parcel_observations[parcel_start:parcel_end] <-
        parcel_observations[parcel_start:parcel_end] + nrow(error)
      encoded <- qpp_encode_tile(q, max_values = max_tile_values,
                                 max_bytes = max_tile_bytes)
      writeBin(encoded, connection, useBytes = TRUE)
      stored_length <- length(encoded)
      row_id <- row_id + 1L
      tile_q_checksums[[row_id]] <- .parcel_mean_numeric_checksum(
        q, "lna_qpp_tile_q_v1"
      )
      tile_payload_checksums[[row_id]] <- .qpp_raw_checksum(
        encoded, "lna_qpp_tile_payload_v1"
      )
      tile_payload_digests[[row_id]] <- .qpp_sha256_raw(encoded)
      index_rows[[row_id]] <- data.frame(
        tile_id = tile_id,
        run_ordinal = as.integer(handle$meta$run_ordinal %||% 0L),
        unit_ordinal = as.integer(handle$meta$run_ordinal %||% 0L),
        time_start = time_start - 1L,
        time_count = nrow(q),
        parcel_start = parcel_start - 1L,
        parcel_count = ncol(q),
        byte_offset = byte_offset,
        byte_length = stored_length,
        uncompressed_length = stored_length,
        inner_header_offset = 0,
        predictor_table_offset = 8,
        outer_codec = 0,
        flags = 1,
        digest_row = row_id - 1L,
        check.names = FALSE
      )
      byte_offset <- byte_offset + stored_length
      logical_values <- logical_values + length(q)
      tile_id <- tile_id + 1L
    }
  }
  payload <- rawConnectionValue(connection)
  index <- do.call(rbind, index_rows)
  rownames(index) <- NULL

  plan <- handle$plan
  fname <- plan$get_next_filename(type)
  base_name <- tools::file_path_sans_ext(fname)
  run_id <- sanitize_run_id(handle$current_run_id %||% "run-01")
  root <- paste0("/scans/", run_id, "/qpp/", base_name)
  paths <- list(
    payload = paste0(root, "/payload"),
    index = paste0(root, "/tile_index"),
    steps = paste0(root, "/parcel_steps"),
    tile_digests = paste0(root, "/tile_sha256"),
    step_digests = paste0(root, "/parcel_step_sha256"),
    receipt = paste0(root, "/quality_receipt.json")
  )
  tile_digest_matrix <- .qpp_digest_matrix(tile_payload_digests)
  step_digest_matrix <- .qpp_digest_matrix(lapply(
    seq_along(delta),
    function(parcel) .qpp_step_digest_raw(delta[[parcel]], parcel - 1L)
  ))
  step_checksum <- .parcel_mean_numeric_checksum(delta, "lna_qpp_steps_v1")
  p_record <- p
  p_record$delta <- NULL
  p_record$allow_clip <- FALSE
  p_record$profile_id <- if (analysis_safe) {
    "analysis-safe"
  } else if (safe_policy) {
    "parcel-series-qpp"
  } else {
    "fitted-run"
  }
  p_record$fit_scope <- fit_scope
  p_record$lossy_support <- "pointwise"
  p_record$coding_dependency <- "chunk"
  p_record$rounding_mode <- "nearest-ties-to-even"
  p_record$integer_domain <- "signed-int32-q"
  p_record$predictor_set <- c("P0", "P1", "P2")
  p_record$predictor_cost <- "anchors32-plus-residual-width-v1"
  p_record$predictor_tie_order <- c("P0", "P1", "P2")
  p_record$byte_order <- "little-endian"
  p_record$bit_order <- "lsb-first"
  p_record$outer_codecs <- "raw"
  p_record$time_count <- nrow(y)
  p_record$parcel_count <- ncol(y)
  p_record$time_chunk <- time_chunk
  p_record$parcel_chunk <- parcel_chunk
  p_record$parcel_steps_checksum <- step_checksum
  p_record$matrix_order <- "time-by-parcel"
  p_record$integrity_policy <- "selected-v1"
  p_record$tile_digest_algorithm <- "sha256-stored-bytes-v1"
  p_record$step_digest_algorithm <- "sha256-binary64-le-v1"
  p_record$parcel_binding <- if (analysis_safe) list(
    descriptor = parcel_assets$descriptor,
    transform = parcel_assets$transform,
    run_id = run_id,
    step_index = as.integer(parcel_assets$step_index),
    asset_id = parcel_assets$asset_id,
    map_checksum = parcel_assets$map_checksum,
    mask_checksum = parcel_assets$mask_checksum,
    labels_path = parcel_assets$labels_path
  ) else NULL
  p_record$normalization_provenance <- .analysis_normalization_record(
    p$normalization_provenance %||%
      parcel_assets$normalization_provenance
  )
  parcelwise_rmse <- sqrt(parcel_squared_error / parcel_observations)
  receipt <- list(
    receipt_version = "1.0.0",
    codec_version = "1.0.0",
    decoder_version = "neuroarchive-qpp-cpp-1",
    units = units,
    quantizer_policy = policy,
    fit_scope = fit_scope,
    maximum_absolute_error = maximum_error,
    parcelwise_RMSE = as.numeric(parcelwise_rmse),
    logical_values = logical_values,
    tile_count = nrow(index),
    stored_payload_bytes = length(payload),
    total_inner_bits = length(payload) * 8,
    effective_bits_per_parcel_tr = length(payload) * 8 / logical_values,
    parcel_steps_checksum = step_checksum,
    payload_checksum = .qpp_canonical_checksum(
      tile_payload_checksums, "lna_qpp_payload_v1"
    ),
    index_checksum = .qpp_canonical_checksum(index, "lna_qpp_index_v1"),
    transform_settings_checksum = .qpp_canonical_checksum(
      p_record, "lna_qpp_transform_settings_v1"
    ),
    q_roundtrip_evidence = list(
      tile_q_checksums = tile_q_checksums,
      combined_checksum = .qpp_canonical_checksum(
        tile_q_checksums, "lna_qpp_all_tile_q_v1"
      )
    ),
    tile_payload_checksums = tile_payload_checksums,
    outer_codec_counts = list(raw = nrow(index)),
    software_versions = list(
      neuroarchive = tryCatch(
        as.character(utils::packageVersion("neuroarchive")),
        error = function(e) "0.1.1"
      ),
      R = paste(R.version$major, R.version$minor, sep = "."),
      decoder = "neuroarchive-qpp-cpp-1"
    )
  )
  receipt$receipt_checksum <- .qpp_canonical_checksum(
    receipt, "lna_qpp_quality_receipt_v1"
  )
  receipt_json <- as.character(jsonlite::toJSON(receipt, auto_unbox = TRUE,
                                                null = "null", digits = NA))
  params_json <- as.character(jsonlite::toJSON(p_record, auto_unbox = TRUE,
                                                null = "null", digits = NA))
  output_key <- desc$outputs[[1]] %||% "qpp_packed"
  desc$version <- "1.1.0"
  desc$codec_version <- "1.0.0"
  desc$run_id <- run_id
  desc$required_decoder <- TRUE
  desc$profile_id <- p_record$profile_id
  desc$decoder_version <- "neuroarchive-qpp-cpp-1"
  desc$params <- p_record
  desc$analysis_contract <- .qpp_analysis_contract(p_record)
  desc$inputs <- c(input_key)
  desc$outputs <- c(output_key)
  desc$datasets <- list(
    list(path = paths$payload, role = "qpp_payload", dtype = "uint8",
         dims = length(payload)),
    list(path = paths$index, role = "qpp_tile_index", dtype = "qpp_index_v1",
         dims = nrow(index)),
    list(path = paths$steps, role = "qpp_parcel_steps", dtype = "float64",
         dims = length(delta)),
    list(path = paths$tile_digests, role = "qpp_tile_sha256", dtype = "uint8",
         dims = dim(tile_digest_matrix)),
    list(path = paths$step_digests, role = "qpp_parcel_step_sha256", dtype = "uint8",
         dims = dim(step_digest_matrix)),
    list(path = paths$receipt, role = "qpp_quality_receipt", dims = 1L)
  )
  plan$add_descriptor(fname, desc)
  .qpp_add_payload(plan, paths$payload, payload, "uint8", "qpp_payload",
                   run_id, plan$next_index - 1L, params_json)
  .qpp_add_payload(plan, paths$index, index, "qpp_index_v1", "qpp_tile_index",
                   run_id, plan$next_index - 1L, params_json)
  .qpp_add_payload(plan, paths$steps, delta, "float64", "qpp_parcel_steps",
                   run_id, plan$next_index - 1L, params_json)
  .qpp_add_payload(plan, paths$tile_digests, tile_digest_matrix, "uint8",
                   "qpp_tile_sha256", run_id, plan$next_index - 1L,
                   params_json)
  .qpp_add_payload(plan, paths$step_digests, step_digest_matrix, "uint8",
                   "qpp_parcel_step_sha256", run_id, plan$next_index - 1L,
                   params_json)
  .qpp_add_payload(plan, paths$receipt, receipt_json, NA_character_,
                   "qpp_quality_receipt", run_id, plan$next_index - 1L,
                   params_json)

  meta <- handle$meta
  meta$qpp_receipt <- receipt
  meta$qpp_analysis_contract <- desc$analysis_contract
  # Remove the uncompressed parcel series. Only terminal payloads remain in Plan.
  handle <- handle$with(meta = meta)
  handle$update_stash(input_key, setNames(list(payload), output_key))
}

#' Inverse step for qpp
#'
#' @param type Transform type string.
#' @param desc Transform descriptor.
#' @param handle DataHandle for one run.
#' @return Updated DataHandle with decoded parcel series.
#' @method invert_step qpp
#' @export
#' @keywords internal
invert_step.qpp <- function(type, desc, handle) {
  legacy <- identical(desc$version, "1.0.0")
  current <- identical(desc$version, "1.1.0") &&
    identical(desc$codec_version, "1.0.0")
  if ((!legacy && !current) ||
      !identical(desc$decoder_version, "neuroarchive-qpp-cpp-1")) {
    abort_lna("unsupported qpp codec or decoder version",
              .subclass = "lna_error_descriptor",
              location = "invert_step.qpp:version")
  }
  if (legacy) {
    warn_lna(
      paste(
        "qpp@1.0.0 lacks exact parcel binding and selected-read integrity;",
        "this decode cannot establish the analysis-safe profile"
      ),
      .subclass = "lna_warning_qpp_legacy_integrity",
      location = "invert_step.qpp:version"
    )
  }
  p <- desc$params %||% list()
  if (current && identical(p$profile_id, "analysis-safe")) {
    .qpp_resolve_bound_parcel(desc, handle)
  } else if (legacy && identical(p$profile_id, "analysis-safe")) {
    p$profile_id <- "parcel-series-qpp"
  }
  entries <- setNames(lapply(
    c("qpp_payload", "qpp_tile_index", "qpp_parcel_steps"),
    function(role) .qpp_dataset(desc, role)
  ), c("payload", "index", "steps"))
  if (current) {
    entries$tile_digests <- .qpp_dataset(desc, "qpp_tile_sha256")
    entries$step_digests <- .qpp_dataset(desc, "qpp_parcel_step_sha256")
  }
  root <- handle$h5[["/"]]
  on.exit(safe_h5_close(root), add = TRUE)
  time_count <- as.integer(p$time_count)
  parcel_count <- as.integer(p$parcel_count)
  if (!is.finite(time_count) || !is.finite(parcel_count) ||
      time_count < 1L || parcel_count < 1L ||
      as.numeric(time_count) * parcel_count > .qpp_default_max_values) {
    abort_lna("qpp decoded dimensions are invalid",
              .subclass = "lna_error_validation",
              location = "invert_step.qpp:dims")
  }
  stats <- .qpp_new_io_stats()
  stats$payload_total_bytes <- .qpp_dataset_length(root, entries$payload$path)
  stats$step_values_total <- .qpp_dataset_length(root, entries$steps$path)
  if (stats$step_values_total != parcel_count) {
    abort_lna("qpp parcel-step extent does not match parcel_count",
              .subclass = "lna_error_validation",
              location = "invert_step.qpp:steps")
  }

  requested_time <- .qpp_normalize_subset_index(
    handle$subset$time_idx %||% handle$subset$time,
    time_count, "time_idx"
  )
  parcel_request <- .qpp_requested_parcels(
    handle, root, parcel_count, stats, desc
  )
  requested_parcels <- parcel_request$requested
  unique_time <- sort(unique(requested_time))
  unique_parcels <- sort(unique(requested_parcels))

  time_chunk <- .qpp_chunk_size(p$time_chunk, "time_chunk")
  parcel_chunk <- .qpp_chunk_size(p$parcel_chunk, "parcel_chunk")
  n_time_tiles <- ceiling(time_count / time_chunk)
  n_parcel_tiles <- ceiling(parcel_count / parcel_chunk)
  index_total <- .qpp_dataset_length(root, entries$index$path)
  if (index_total != n_time_tiles * n_parcel_tiles) {
    abort_lna("qpp index extent contradicts deterministic tile layout",
              .subclass = "lna_error_validation",
              location = "invert_step.qpp:index_extent")
  }
  stats$index_rows_total <- index_total
  if (current) {
    tile_digest_dims <- .qpp_dataset_dims(root, entries$tile_digests$path)
    step_digest_dims <- .qpp_dataset_dims(root, entries$step_digests$path)
    if (!identical(tile_digest_dims, c(as.numeric(index_total), 32)) ||
        !identical(step_digest_dims, c(as.numeric(parcel_count), 32))) {
      abort_lna(
        "qpp digest dataset extents contradict tile or parcel counts",
        code = "qpp_digest_shape",
        .subclass = "lna_error_qpp_integrity",
        location = "invert_step.qpp:digest_extent"
      )
    }
    stats$tile_digest_rows_total <- index_total
    stats$step_digest_rows_total <- parcel_count
  }
  tile_rows <- .qpp_required_tile_rows(
    unique_time, unique_parcels, time_chunk, parcel_chunk,
    n_parcel_tiles
  )
  stats$index_rows_read <- length(tile_rows)
  stats$touched_tile_rows <- as.integer(tile_rows)
  q_unique <- matrix(NA_real_, length(unique_time), length(unique_parcels))

  cache <- handle$meta$qpp_tile_cache
  if (!inherits(cache, "environment")) cache <- .qpp_new_tile_cache()
  if (length(tile_rows)) {
    index <- .qpp_read_subset(
      root, entries$index$path, list(tile_rows), stats, "index"
    )
    if (!is.data.frame(index) || nrow(index) != length(tile_rows)) {
      abort_lna("qpp index subset has an invalid storage type",
                .subclass = "lna_error_validation",
                location = "invert_step.qpp:index_storage")
    }
    .qpp_validate_selected_offsets(
      index, tile_rows, index_total, stats$payload_total_bytes
    )
    tile_digests <- if (current) {
      .qpp_read_digest_rows(
        root, entries$tile_digests$path, tile_rows, stats, "tile_digests"
      )
    } else {
      NULL
    }
    if (current) stats$tile_digest_rows_read <- length(tile_rows)
    for (i in seq_along(tile_rows)) {
      row <- index[i, , drop = FALSE]
      row_number <- tile_rows[[i]]
      .qpp_validate_index_row(
        row, row_number, time_count, parcel_count, time_chunk,
        parcel_chunk, n_parcel_tiles, stats$payload_total_bytes
      )
      if (current &&
          (as.numeric(row$flags) != 1 ||
           as.numeric(row$digest_row) != row_number - 1L)) {
        abort_lna(
          "qpp index digest reference is invalid",
          code = "qpp_digest_index",
          .subclass = "lna_error_qpp_integrity",
          location = "invert_step.qpp:index_digest"
        )
      }
      expected_digest <- if (current) {
        as.raw(tile_digests[i, , drop = TRUE])
      } else {
        raw()
      }
      cache_key <- paste(
        entries$payload$path, as.numeric(row$tile_id),
        as.numeric(row$byte_offset), as.numeric(row$byte_length),
        if (current) paste(sprintf("%02x", as.integer(expected_digest)), collapse = "") else "legacy",
        sep = "|"
      )
      tile <- .qpp_cache_get(cache, cache_key)
      if (is.null(tile)) {
        offset <- as.numeric(row$byte_offset)
        byte_length <- as.numeric(row$byte_length)
        bytes <- .qpp_payload_as_raw(.qpp_read_subset(
          root, entries$payload$path,
          list(seq.int(offset + 1, offset + byte_length)),
          stats, "payload"
        ))
        stats$payload_bytes_read <- stats$payload_bytes_read + byte_length
        stats$payload_read_calls <- stats$payload_read_calls + 1L
        if (current && !identical(.qpp_sha256_raw(bytes), expected_digest)) {
          abort_lna(
            "qpp selected tile digest mismatch",
            code = "qpp_tile_digest_mismatch",
            .subclass = "lna_error_qpp_integrity",
            location = "invert_step.qpp:tile_digest"
          )
        }
        tile <- qpp_decode_tile(
          bytes,
          max_time = as.integer(row$time_count),
          max_parcels = as.integer(row$parcel_count),
          max_values = p$max_tile_values %||% .qpp_default_max_values,
          max_bytes = p$max_tile_bytes %||% .qpp_default_max_bytes
        )
        stats$tiles_decoded <- stats$tiles_decoded + 1L
        .qpp_cache_put(cache, cache_key, tile)
      } else {
        stats$tile_cache_hits <- stats$tile_cache_hits + 1L
      }
      t_start <- as.integer(row$time_start) + 1L
      p_start <- as.integer(row$parcel_start) + 1L
      t_values <- seq.int(t_start, length.out = nrow(tile))
      p_values <- seq.int(p_start, length.out = ncol(tile))
      wanted_t <- which(unique_time %in% t_values)
      wanted_p <- which(unique_parcels %in% p_values)
      if (length(wanted_t) && length(wanted_p)) {
        local_t <- match(unique_time[wanted_t], t_values)
        local_p <- match(unique_parcels[wanted_p], p_values)
        q_unique[wanted_t, wanted_p] <- tile[local_t, local_p, drop = FALSE]
      }
    }
  }
  if (length(q_unique) && anyNA(q_unique)) {
    abort_lna("qpp selected tiles do not cover the requested parcel matrix",
              .subclass = "lna_error_validation",
              location = "invert_step.qpp:coverage")
  }

  delta_unique <- if (length(unique_parcels)) {
    as.numeric(.qpp_read_subset(
      root, entries$steps$path, list(unique_parcels), stats, "steps"
    ))
  } else {
    numeric()
  }
  stats$step_values_read <- length(unique_parcels)
  if (length(delta_unique) != length(unique_parcels) ||
      any(!is.finite(delta_unique)) || any(delta_unique <= 0)) {
    abort_lna("qpp selected parcel steps are invalid",
              .subclass = "lna_error_validation",
              location = "invert_step.qpp:steps")
  }
  if (current && length(unique_parcels)) {
    step_digests <- .qpp_read_digest_rows(
      root, entries$step_digests$path, unique_parcels, stats, "step_digests"
    )
    stats$step_digest_rows_read <- length(unique_parcels)
    for (i in seq_along(unique_parcels)) {
      actual <- .qpp_step_digest_raw(
        delta_unique[[i]], unique_parcels[[i]] - 1L
      )
      expected <- as.raw(step_digests[i, , drop = TRUE])
      if (!identical(actual, expected)) {
        abort_lna(
          "qpp selected parcel step digest mismatch",
          code = "qpp_step_digest_mismatch",
          .subclass = "lna_error_qpp_integrity",
          location = "invert_step.qpp:step_digest"
        )
      }
    }
  }
  full_steps_verified <- identical(unique_parcels, seq_len(parcel_count))
  if (full_steps_verified) {
    checksum <- .parcel_mean_numeric_checksum(delta_unique, "lna_qpp_steps_v1")
    if (!identical(checksum, p$parcel_steps_checksum)) {
      abort_lna("qpp parcel step checksum mismatch",
                .subclass = "lna_error_validation",
                location = "invert_step.qpp:steps")
    }
  }

  q <- if (length(requested_time) && length(requested_parcels)) {
    q_unique[
      match(requested_time, unique_time),
      match(requested_parcels, unique_parcels),
      drop = FALSE
    ]
  } else {
    matrix(numeric(), nrow = length(requested_time),
           ncol = length(requested_parcels))
  }
  delta <- if (length(requested_parcels)) {
    delta_unique[match(requested_parcels, unique_parcels)]
  } else {
    numeric()
  }
  parcel_series <- if (ncol(q)) sweep(q, 2L, delta, "*") else q
  if (!is.null(handle$subset$time_idx %||% handle$subset$time)) {
    handle$subset_applied$time_idx <- TRUE
  }
  meta <- handle$meta
  meta$qpp_decoded_q_checksum <- .parcel_mean_numeric_checksum(q, "lna_qpp_q_v1")
  meta$qpp_units <- p$units
  meta$qpp_selected_parcels <- as.integer(requested_parcels)
  meta$qpp_roi_dense_labels <- parcel_request$roi_dense_labels
  meta$qpp_tile_cache <- cache
  meta$qpp_io_receipt <- .qpp_io_receipt(
    stats, cache, full_steps_verified,
    time_count, parcel_count, requested_time, requested_parcels
  )
  handle <- handle$with(meta = meta)
  output_key <- desc$outputs[[1]] %||% "qpp_packed"
  input_key <- desc$inputs[[1]] %||% "parcel_series"
  handle$update_stash(output_key, setNames(list(parcel_series), input_key))
}

.qpp_new_io_stats <- function() {
  state <- new.env(parent = emptyenv())
  state$index_rows_total <- 0
  state$index_rows_read <- 0
  state$payload_total_bytes <- 0
  state$payload_bytes_read <- 0
  state$payload_read_calls <- 0L
  state$step_values_total <- 0
  state$step_values_read <- 0L
  state$tile_digest_rows_total <- 0
  state$tile_digest_rows_read <- 0L
  state$step_digest_rows_total <- 0
  state$step_digest_rows_read <- 0L
  state$label_values_total <- 0
  state$label_values_read <- 0L
  state$tiles_decoded <- 0L
  state$tile_cache_hits <- 0L
  state$touched_tile_rows <- integer()
  state$read_calls <- list()
  state
}

.qpp_new_tile_cache <- function(max_bytes = 64 * 1024^2) {
  max_bytes <- as.numeric(max_bytes)
  if (length(max_bytes) != 1L || !is.finite(max_bytes) || max_bytes < 0) {
    abort_lna("qpp tile cache max_bytes must be finite and non-negative",
              .subclass = "lna_error_validation",
              location = "qpp_tile_cache:max_bytes")
  }
  state <- new.env(parent = emptyenv())
  state$entries <- new.env(parent = emptyenv())
  state$order <- character()
  state$entry_bytes <- numeric()
  state$bytes <- 0
  state$max_bytes <- max_bytes
  state$evictions <- 0L
  state
}

.qpp_tile_cache_clear <- function(cache) {
  if (!inherits(cache, "environment")) return(invisible(NULL))
  if (inherits(cache$entries, "environment")) {
    keys <- ls(cache$entries, all.names = TRUE)
    if (length(keys)) rm(list = keys, envir = cache$entries)
  }
  cache$order <- character()
  cache$entry_bytes <- numeric()
  cache$bytes <- 0
  invisible(cache)
}

.qpp_cache_get <- function(cache, key) {
  if (!inherits(cache, "environment") ||
      !inherits(cache$entries, "environment") ||
      !exists(key, envir = cache$entries, inherits = FALSE)) return(NULL)
  value <- get(key, envir = cache$entries, inherits = FALSE)
  cache$order <- c(setdiff(cache$order, key), key)
  value
}

.qpp_cache_put <- function(cache, key, value) {
  if (!inherits(cache, "environment") || cache$max_bytes <= 0) {
    return(invisible(cache))
  }
  bytes <- as.numeric(utils::object.size(value))
  if (bytes > cache$max_bytes) return(invisible(cache))
  if (exists(key, envir = cache$entries, inherits = FALSE)) {
    cache$bytes <- cache$bytes - unname(cache$entry_bytes[[key]])
  }
  assign(key, value, envir = cache$entries)
  cache$entry_bytes[[key]] <- bytes
  cache$order <- c(setdiff(cache$order, key), key)
  cache$bytes <- cache$bytes + bytes
  while (cache$bytes > cache$max_bytes && length(cache$order)) {
    victim <- cache$order[[1L]]
    cache$order <- cache$order[-1L]
    cache$bytes <- cache$bytes - unname(cache$entry_bytes[[victim]])
    cache$entry_bytes <- cache$entry_bytes[names(cache$entry_bytes) != victim]
    if (exists(victim, envir = cache$entries, inherits = FALSE)) {
      rm(list = victim, envir = cache$entries)
    }
    cache$evictions <- cache$evictions + 1L
  }
  invisible(cache)
}

.qpp_dataset_length <- function(root, path) {
  if (!path_exists_safely(root, path)) {
    abort_lna(paste0("qpp dataset is missing: ", path),
              .subclass = "lna_error_missing_path",
              location = "invert_step.qpp:dataset")
  }
  dset <- root[[path]]
  on.exit(safe_h5_close(dset), add = TRUE)
  dims <- as.numeric(dset$dims)
  if (length(dims) != 1L || !is.finite(dims) || dims < 0) {
    abort_lna(paste0("qpp dataset must be one-dimensional: ", path),
              .subclass = "lna_error_validation",
              location = "invert_step.qpp:dataset_extent")
  }
  dims
}

.qpp_dataset_dims <- function(root, path) {
  if (!path_exists_safely(root, path)) {
    abort_lna(
      paste0("qpp dataset is missing: ", path),
      code = "qpp_digest_missing",
      .subclass = "lna_error_qpp_integrity",
      location = "invert_step.qpp:digest_dataset"
    )
  }
  dset <- root[[path]]
  on.exit(safe_h5_close(dset), add = TRUE)
  as.numeric(dset$dims)
}

.qpp_read_digest_rows <- function(root, path, rows, stats, kind) {
  if (!length(rows)) return(matrix(integer(), nrow = 0L, ncol = 32L))
  values <- .qpp_read_subset(
    root, path, list(as.integer(rows), seq_len(32L)), stats, kind
  )
  values <- as.integer(values)
  if (length(values) != length(rows) * 32L ||
      anyNA(values) || any(values < 0L) || any(values > 255L)) {
    abort_lna(
      "qpp digest rows are malformed",
      code = "qpp_digest_shape",
      .subclass = "lna_error_qpp_integrity",
      location = "invert_step.qpp:digest_rows"
    )
  }
  matrix(values, nrow = length(rows), ncol = 32L)
}

.qpp_read_subset <- function(root, path, index, stats, kind) {
  selected <- as.numeric(index[[1L]])
  selection <- list(
    count = length(selected),
    first = if (length(selected)) selected[[1L]] else NULL,
    last = if (length(selected)) selected[[length(selected)]] else NULL,
    contiguous = length(selected) < 2L || all(diff(selected) == 1)
  )
  if (length(selected) <= 64L) selection$values <- selected
  stats$read_calls[[length(stats$read_calls) + 1L]] <- list(
    kind = kind, path = path, selection = selection,
    access_method = "H5D$read(args=...)"
  )
  h5_read_subset(root, path, index)
}

.qpp_normalize_subset_index <- function(index, size, field) {
  if (is.null(index)) return(seq_len(size))
  if (!is.numeric(index) || anyNA(index) || any(!is.finite(index)) ||
      any(index != trunc(index)) || any(index < 1) || any(index > size)) {
    abort_lna(paste0(field, " is outside qpp data or is not an exact integer"),
              .subclass = "lna_error_validation",
              location = paste0("invert_step.qpp:", field))
  }
  as.integer(index)
}

.qpp_requested_parcels <- function(handle, root, parcel_count, stats, desc) {
  roi <- handle$subset$roi_mask %||% handle$subset$roi
  if (is.null(roi)) {
    requested <- .qpp_normalize_subset_index(
      handle$subset$parcel_idx, parcel_count, "parcel_idx"
    )
    return(list(requested = requested, roi_dense_labels = NULL))
  }
  binding <- desc$params$parcel_binding %||% NULL
  if (is.null(binding)) {
    abort_lna(
      "ROI reads require a qpp descriptor bound to a spatial parcel asset",
      code = "qpp_roi_without_parcel_binding",
      .subclass = "lna_error_qpp_binding",
      location = "invert_step.qpp:roi_mask"
    )
  }
  .qpp_validate_parcel_binding_fields(
    binding, run_id = desc$run_id %||% handle$current_run_id,
    location = "invert_step.qpp:parcel_binding"
  )
  label_path <- binding$labels_path
  label_count <- .qpp_dataset_length(root, label_path)
  stats$label_values_total <- label_count
  roi <- as.logical(roi)
  mask <- handle$mask_info$mask %||% NULL
  if (!is.null(mask) && length(roi) == length(mask)) {
    active_roi <- roi[as.logical(mask)]
  } else if (length(roi) == label_count) {
    active_roi <- roi
  } else {
    abort_lna("roi_mask does not match the archive mask or active voxel order",
              .subclass = "lna_error_validation",
              location = "invert_step.qpp:roi_mask")
  }
  positions <- which(active_roi)
  labels <- if (length(positions)) {
    as.integer(.qpp_read_subset(
      root, label_path, list(positions), stats, "parcel_labels"
    ))
  } else {
    integer()
  }
  stats$label_values_read <- length(positions)
  if (length(labels) != length(positions) || anyNA(labels) ||
      any(labels < 1L) || any(labels > parcel_count)) {
    abort_lna("selected parcel labels are invalid",
              .subclass = "lna_error_validation",
              location = "invert_step.qpp:parcel_labels")
  }
  list(requested = sort(unique(labels)), roi_dense_labels = labels)
}

.qpp_resolve_bound_parcel <- function(desc, handle) {
  binding <- desc$params$parcel_binding %||% NULL
  run_id <- desc$run_id %||% handle$current_run_id %||% NULL
  .qpp_validate_parcel_binding_fields(
    binding, run_id = run_id, location = "invert_step.qpp:parcel_binding"
  )
  descriptors <- handle$meta$transform_descriptors %||% list()
  matches <- which(names(descriptors) == binding$descriptor)
  if (length(matches) != 1L) {
    abort_lna(
      "analysis-safe qpp parcel binding does not resolve exactly one descriptor",
      code = if (length(matches)) "qpp_parcel_duplicate" else "qpp_parcel_missing",
      .subclass = "lna_error_qpp_binding",
      location = "invert_step.qpp:parcel_binding"
    )
  }
  parcel <- descriptors[[matches]]
  p <- parcel$params %||% list()
  labels <- Filter(
    function(dataset) identical(dataset$role, "parcel_labels"),
    parcel$datasets %||% list()
  )
  valid <- identical(parcel$type, "parcel.mean") &&
    identical(parcel$version, "1.1.0") &&
    identical(parcel$run_id, binding$run_id) &&
    identical(as.integer(p$step_index), as.integer(binding$step_index)) &&
    identical(p$asset_id, binding$asset_id) &&
    identical(p$map_checksum, binding$map_checksum) &&
    identical(p$mask_checksum, binding$mask_checksum) &&
    length(labels) == 1L && identical(labels[[1L]]$path, binding$labels_path)
  if (!isTRUE(valid)) {
    abort_lna(
      "analysis-safe qpp parcel descriptor contradicts its recorded binding",
      code = "qpp_parcel_mismatch",
      .subclass = "lna_error_qpp_binding",
      location = "invert_step.qpp:parcel_binding"
    )
  }
  invisible(parcel)
}

.qpp_required_tile_rows <- function(time_index, parcel_index, time_chunk,
                                    parcel_chunk, n_parcel_tiles) {
  if (!length(time_index) || !length(parcel_index)) return(integer())
  time_tiles <- sort(unique((time_index - 1L) %/% time_chunk + 1L))
  parcel_tiles <- sort(unique((parcel_index - 1L) %/% parcel_chunk + 1L))
  rows <- as.vector(outer(
    time_tiles, parcel_tiles,
    function(t, p) (t - 1L) * n_parcel_tiles + p
  ))
  sort(unique(as.integer(rows)))
}

.qpp_validate_index_row <- function(row, row_number, time_count, parcel_count,
                                    time_chunk, parcel_chunk, n_parcel_tiles,
                                    payload_length) {
  time_ordinal <- (row_number - 1L) %/% n_parcel_tiles
  parcel_ordinal <- (row_number - 1L) %% n_parcel_tiles
  expected_time_start <- time_ordinal * time_chunk
  expected_parcel_start <- parcel_ordinal * parcel_chunk
  expected_time_count <- min(time_chunk, time_count - expected_time_start)
  expected_parcel_count <- min(parcel_chunk, parcel_count - expected_parcel_start)
  offset <- as.numeric(row$byte_offset)
  byte_length <- as.numeric(row$byte_length)
  valid <-
    as.numeric(row$tile_id) == row_number - 1L &&
    as.numeric(row$time_start) == expected_time_start &&
    as.numeric(row$parcel_start) == expected_parcel_start &&
    as.numeric(row$time_count) == expected_time_count &&
    as.numeric(row$parcel_count) == expected_parcel_count &&
    is.finite(offset) && is.finite(byte_length) && offset >= 0 &&
    byte_length >= 1 && offset + byte_length <= payload_length &&
    as.numeric(row$outer_codec) == 0 &&
    as.numeric(row$uncompressed_length) == byte_length
  if (!isTRUE(valid)) {
    abort_lna("qpp index row contradicts its deterministic tile or byte range",
              .subclass = "lna_error_validation",
              location = "invert_step.qpp:index")
  }
  invisible(TRUE)
}

.qpp_validate_selected_offsets <- function(index, tile_rows, index_total,
                                           payload_length) {
  offsets <- as.numeric(index$byte_offset)
  lengths <- as.numeric(index$byte_length)
  ends <- offsets + lengths
  if (length(offsets) > 1L &&
      (any(diff(offsets) <= 0) || any(offsets[-1L] < ends[-length(ends)]))) {
    abort_lna("selected qpp payload ranges overlap or are not monotone",
              .subclass = "lna_error_validation",
              location = "invert_step.qpp:index_ranges")
  }
  if (identical(as.integer(tile_rows), seq_len(index_total))) {
    exact <- offsets[[1L]] == 0 &&
      (length(offsets) == 1L || all(offsets[-1L] == ends[-length(ends)])) &&
      ends[[length(ends)]] == payload_length
    if (!isTRUE(exact)) {
      abort_lna("full qpp index does not exactly cover its payload",
                .subclass = "lna_error_validation",
                location = "invert_step.qpp:index_coverage")
    }
  }
  invisible(TRUE)
}

.qpp_io_receipt <- function(stats, cache, full_steps_verified, time_count,
                            parcel_count, requested_time, requested_parcels) {
  list(
    access_method = "H5D$read(args=...)",
    index_rows_read = stats$index_rows_read,
    index_rows_total = stats$index_rows_total,
    payload_bytes_read = stats$payload_bytes_read,
    payload_total_bytes = stats$payload_total_bytes,
    payload_read_calls = stats$payload_read_calls,
    step_values_read = stats$step_values_read,
    step_values_total = stats$step_values_total,
    tile_digest_rows_read = stats$tile_digest_rows_read,
    tile_digest_rows_total = stats$tile_digest_rows_total,
    step_digest_rows_read = stats$step_digest_rows_read,
    step_digest_rows_total = stats$step_digest_rows_total,
    selected_tile_digests_verified = stats$tile_digest_rows_read == stats$index_rows_read,
    selected_step_digests_verified = stats$step_digest_rows_read == stats$step_values_read,
    label_values_read = stats$label_values_read,
    label_values_total = stats$label_values_total,
    tiles_decoded = stats$tiles_decoded,
    tile_cache_hits = stats$tile_cache_hits,
    touched_tile_rows = stats$touched_tile_rows,
    requested_time_values = length(requested_time),
    requested_parcel_values = length(requested_parcels),
    full_step_checksum_verified = isTRUE(full_steps_verified),
    physical_time_partial = length(unique(requested_time)) < time_count,
    physical_parcel_partial = length(unique(requested_parcels)) < parcel_count,
    cache_bytes = cache$bytes,
    cache_max_bytes = cache$max_bytes,
    cache_entries = length(cache$order),
    cache_evictions = cache$evictions,
    read_calls = stats$read_calls
  )
}

.qpp_validate_delta <- function(delta, parcel_count) {
  if (!is.numeric(delta) || !length(delta) || any(!is.finite(delta)) || any(delta <= 0)) {
    abort_lna("fixed/external qpp requires finite positive delta",
              .subclass = "lna_error_validation",
              location = "forward_step.qpp:delta")
  }
  if (length(delta) == 1L) delta <- rep(as.numeric(delta), parcel_count)
  if (length(delta) != parcel_count) {
    abort_lna("qpp delta must be scalar or one value per parcel",
              .subclass = "lna_error_validation",
              location = "forward_step.qpp:delta")
  }
  as.numeric(delta)
}

.qpp_quantize_pointwise <- function(y, delta) {
  q <- round(sweep(y, 2L, delta, "/"))
  # QPP's integer domain has one zero. R can preserve a negative IEEE-754
  # sign bit when a small negative value rounds to zero, so canonicalize it
  # before q checksums are recorded.
  q[q == 0] <- 0
  if (any(!is.finite(q)) || any(q != trunc(q)) ||
      any(q < -2147483648) || any(q > 2147483647)) {
    abort_lna("qpp quantized values exceed signed-int32; increase delta",
              .subclass = "lna_error_validation",
              location = "forward_step.qpp:integer_domain")
  }
  q
}

.qpp_chunk_size <- function(value, field) {
  value <- as.numeric(value)
  if (length(value) != 1L || !is.finite(value) || value != trunc(value) ||
      value < 1 || value > 65535) {
    abort_lna(paste0(field, " must be an integer in 1..65535"),
              .subclass = "lna_error_validation",
              location = paste0("forward_step.qpp:", field))
  }
  as.integer(value)
}

.qpp_add_payload <- function(plan, path, value, dtype, role, run_id,
                             step_index, params_json) {
  plan$add_payload(path, value)
  plan$add_dataset_def(path, role, "qpp", run_id, as.integer(step_index),
                       params_json, path, "eager", dtype = dtype)
  invisible(plan)
}

.qpp_dataset <- function(desc, role) {
  entries <- Filter(function(x) identical(x$role, role), desc$datasets %||% list())
  if (length(entries) != 1L) {
    abort_lna(paste0("qpp descriptor requires one ", role, " dataset"),
              .subclass = "lna_error_descriptor",
              location = "invert_step.qpp:descriptor")
  }
  entries[[1L]]
}

.qpp_payload_as_raw <- function(payload) {
  if (is.raw(payload)) return(payload)
  if (!is.numeric(payload) || any(!is.finite(payload)) ||
      any(payload != trunc(payload)) || any(payload < 0) || any(payload > 255)) {
    abort_lna("qpp uint8 payload contains invalid byte values",
              .subclass = "lna_error_validation",
              location = "invert_step.qpp:payload")
  }
  as.raw(payload)
}

.qpp_raw_checksum <- function(value, prefix) {
  value <- .qpp_payload_as_raw(value)
  paste0(
    "sha256:",
    digest::digest(c(charToRaw(prefix), value), algo = "sha256", serialize = FALSE)
  )
}

.qpp_sha256_raw <- function(value) {
  value <- .qpp_payload_as_raw(value)
  digest::digest(
    value, algo = "sha256", serialize = FALSE, raw = TRUE
  )
}

.qpp_digest_matrix <- function(rows) {
  if (!is.list(rows) || !length(rows) ||
      any(vapply(rows, function(x) !is.raw(x) || length(x) != 32L, logical(1)))) {
    abort_lna(
      "qpp digest rows must contain exactly 32 raw SHA-256 bytes",
      code = "qpp_digest_shape",
      .subclass = "lna_error_qpp_integrity",
      location = ".qpp_digest_matrix"
    )
  }
  matrix(
    as.integer(do.call(c, rows)), nrow = length(rows), ncol = 32L,
    byrow = TRUE
  )
}

.qpp_step_digest_raw <- function(value, parcel_ordinal) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value <= 0 || !is.numeric(parcel_ordinal) || length(parcel_ordinal) != 1L ||
      !is.finite(parcel_ordinal) || parcel_ordinal != trunc(parcel_ordinal) ||
      parcel_ordinal < 0 || parcel_ordinal > 2^32 - 1) {
    abort_lna(
      "qpp step digest input is invalid",
      code = "qpp_step_digest",
      .subclass = "lna_error_qpp_integrity",
      location = ".qpp_step_digest_raw"
    )
  }
  connection <- rawConnection(raw(), "wb")
  on.exit(close(connection), add = TRUE)
  writeBin(charToRaw("lna_qpp_step_digest_v1"), connection)
  writeBin(as.integer(parcel_ordinal), connection, size = 4L, endian = "little")
  writeBin(as.double(value), connection, size = 8L, endian = "little")
  .qpp_sha256_raw(rawConnectionValue(connection))
}

.qpp_validate_parcel_binding_fields <- function(binding, run_id = NULL,
                                                location = "qpp:parcel_binding") {
  required <- c(
    "descriptor", "transform", "run_id", "step_index", "asset_id",
    "map_checksum", "mask_checksum", "labels_path"
  )
  if (!is.list(binding) || any(vapply(required, function(field) {
    is.null(binding[[field]])
  }, logical(1)))) {
    abort_lna(
      "analysis-safe qpp requires a complete parcel binding",
      code = "qpp_parcel_binding",
      .subclass = "lna_error_qpp_binding",
      location = location
    )
  }
  if (!is.character(binding$descriptor) || length(binding$descriptor) != 1L ||
      !grepl("^[0-9]+_parcel[.]mean[.]json$", binding$descriptor) ||
      !identical(binding$transform, "parcel.mean@1.1.0") ||
      !is.character(binding$run_id) || length(binding$run_id) != 1L ||
      !nzchar(binding$run_id) ||
      !is.numeric(binding$step_index) || length(binding$step_index) != 1L ||
      !is.finite(binding$step_index) || binding$step_index != trunc(binding$step_index) ||
      binding$step_index < 0) {
    abort_lna(
      "analysis-safe qpp parcel binding has invalid descriptor, version, run, or step",
      code = "qpp_parcel_binding",
      .subclass = "lna_error_qpp_binding",
      location = location
    )
  }
  if (!is.null(run_id) && !identical(binding$run_id, run_id)) {
    abort_lna(
      "analysis-safe qpp parcel binding belongs to a different run",
      code = "qpp_parcel_wrong_run",
      .subclass = "lna_error_qpp_binding",
      location = location
    )
  }
  for (field in c("asset_id", "map_checksum", "mask_checksum")) {
    .lna_canonical_sha256(binding[[field]], field, location = location)
  }
  expected_labels <- paste0(
    "/spatial/qpp_assets/",
    sub("^sha256:", "", binding$asset_id),
    "/parcel_labels"
  )
  if (!is.character(binding$labels_path) || length(binding$labels_path) != 1L ||
      !identical(binding$labels_path, expected_labels)) {
    abort_lna(
      "analysis-safe qpp parcel binding label path contradicts asset identity",
      code = "qpp_parcel_wrong_path",
      .subclass = "lna_error_qpp_binding",
      location = location
    )
  }
  invisible(TRUE)
}

.qpp_canonical_checksum <- function(value, prefix) {
  json <- jsonlite::toJSON(
    value, auto_unbox = TRUE, null = "null", digits = NA,
    dataframe = "columns", POSIXt = "ISO8601"
  )
  .qpp_raw_checksum(c(charToRaw(prefix), charToRaw(as.character(json))),
                    "lna_qpp_canonical_json_v1")
}

#' Default qpp parameters
#' @return QPP defaults; delta and scientific units are intentionally required.
#' @keywords internal
lna_default.qpp <- function() {
  list(
    delta = NULL,
    quantizer_policy = "fixed",
    units = NULL,
    calibration_id = NULL,
    calibration_checksum = NULL,
    normalization_provenance = NULL,
    fit_id = NULL,
    target_q = 32767L,
    zero_step = 1,
    allow_clip = FALSE,
    time_chunk = 128L,
    parcel_chunk = 128L,
    max_tile_values = .qpp_default_max_values,
    max_tile_bytes = .qpp_default_max_bytes
  )
}
