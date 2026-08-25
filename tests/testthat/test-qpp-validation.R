make_qpp_validation_file <- function(
    y = outer(1:12, 1:4, function(t, k) sin(t / 3) + k),
    time_chunk = 4L, parcel_chunk = 2L,
    policy = "fixed", with_parcel = FALSE) {
  path <- tempfile(fileext = ".h5")
  if (isTRUE(with_parcel)) {
    mask <- array(TRUE, dim = c(8, 1, 1))
    x <- array(seq_len(8 * 12), dim = c(8, 1, 1, 12))
    labels <- array(rep(1:4, each = 2), dim = dim(mask))
    affine <- diag(c(2, 2, 2, 1))
    header <- list(
      dims = dim(mask), spacing = c(2, 2, 2), origin = c(0, 0, 0),
      transform = affine, spatial_units = "mm",
      voxel_ordering = "first-axis-fastest"
    )
    result <- core_write(
      x, transforms = c("parcel.mean", "qpp"), mask = mask, header = header,
      transform_params = list(
        "parcel.mean" = list(
          labels = labels, fit_scope = "fixed", source_id = "map-v1",
          signal_units = "psc", spatial_units = "mm",
          voxel_ordering = "first-axis-fastest",
          map_grid_id = .lna_grid_id(
            dim(mask), affine, "mm", "first-axis-fastest"
          )
        ),
        qpp = list(
          delta = 0.01, quantizer_policy = "fixed", units = "psc",
          time_chunk = time_chunk, parcel_chunk = parcel_chunk
        )
      )
    )
  } else {
    params <- list(
      delta = 0.01, quantizer_policy = policy, units = "psc",
      time_chunk = time_chunk, parcel_chunk = parcel_chunk
    )
    if (identical(policy, "external-calibration")) {
      params$calibration_id <- "external-cal-v1"
      params$calibration_checksum <- paste0("sha256:", strrep("a", 64))
    }
    handle <- DataHandle$new(
      initial_stash = list(parcel_series = as.matrix(y)),
      initial_meta = list(parcel_mean_assets = list(signal_units = "psc")),
      plan = Plan$new(), run_ids = "run-01", current_run_id = "run-01"
    )
    desc <- list(
      type = "qpp", params = params,
      inputs = list("parcel_series"), outputs = list("qpp_packed")
    )
    transformed <- run_transform_step("forward", "qpp", desc, handle, 0L)
    result <- list(plan = transformed$plan)
    header <- NULL
  }
  h5 <- neuroarchive:::open_h5(path, "w")
  materialise_plan(h5, result$plan, checksum = "none", header = header,
                   run_ids = "run-01")
  neuroarchive:::close_h5_safely(h5)
  path
}

with_qpp_mutation <- function(file, code) {
  h5 <- neuroarchive:::open_h5(file, "r+")
  on.exit(neuroarchive:::close_h5_safely(h5), add = TRUE)
  group <- h5[["transforms"]]
  on.exit(safe_h5_close(group), add = TRUE)
  transforms <- discover_transforms(group)
  name <- transforms$name[transforms$type == "qpp"][[1L]]
  desc <- read_json_descriptor(group, name)
  root <- h5[["/"]]
  on.exit(safe_h5_close(root), add = TRUE)
  code(h5, root, group, name, desc)
  h5$flush()
  invisible(file)
}

mutate_qpp_descriptor <- function(file, mutate) {
  with_qpp_mutation(file, function(h5, root, group, name, desc) {
    write_json_descriptor(group, name, mutate(desc))
  })
  invisible(file)
}

qpp_validation_path <- function(desc, role) {
  Filter(function(x) identical(x$role, role), desc$datasets)[[1L]]$path
}

mutate_qpp_index <- function(file, mutate) {
  with_qpp_mutation(file, function(h5, root, group, name, desc) {
    path <- qpp_validation_path(desc, "qpp_tile_index")
    dataset <- h5[[path]]
    on.exit(safe_h5_close(dataset), add = TRUE)
    index <- dataset$read()
    index <- mutate(index)
    dataset$write(args = list(seq_len(nrow(index))), value = index)
  })
}

mutate_qpp_payload <- function(file, mutate) {
  with_qpp_mutation(file, function(h5, root, group, name, desc) {
    path <- qpp_validation_path(desc, "qpp_payload")
    dataset <- h5[[path]]
    on.exit(safe_h5_close(dataset), add = TRUE)
    payload <- .qpp_payload_as_raw(dataset$read())
    payload <- mutate(payload)
    dataset$write(
      args = list(seq_along(payload)), value = as.integer(payload)
    )
  })
}

refresh_qpp_tile_digests <- function(file) {
  with_qpp_mutation(file, function(h5, root, group, name, desc) {
    payload <- .qpp_payload_as_raw(h5_read(
      root, qpp_validation_path(desc, "qpp_payload")
    ))
    index <- h5_read(root, qpp_validation_path(desc, "qpp_tile_index"))
    rows <- lapply(seq_len(nrow(index)), function(row) {
      start <- as.numeric(index$byte_offset[[row]]) + 1L
      length <- as.numeric(index$byte_length[[row]])
      .qpp_sha256_raw(payload[seq.int(start, length.out = length)])
    })
    digests <- .qpp_digest_matrix(rows)
    dataset <- h5[[qpp_validation_path(desc, "qpp_tile_sha256")]]
    on.exit(safe_h5_close(dataset), add = TRUE)
    dataset$write(
      args = list(seq_len(nrow(digests)), seq_len(ncol(digests))),
      value = digests
    )
  })
  invisible(file)
}

expect_qpp_validation_error <- function(file, location, pattern = NULL,
                                        qpp = "all",
                                        class = "lna_error_qpp_corrupt", ...) {
  error <- rlang::catch_cnd(validate_lna(file, checksum = FALSE, qpp = qpp, ...))
  expect_s3_class(error, class)
  expect_identical(error$location, location)
  if (!is.null(pattern)) expect_match(conditionMessage(error), pattern)
  invisible(error)
}

test_that("QPP validation returns reproducible structural, sample, and full receipts", {
  file <- make_qpp_validation_file()
  on.exit(unlink(file), add = TRUE)
  structural <- validate_lna(
    file, checksum = FALSE, qpp = "structural", receipt = TRUE
  )
  sampled <- validate_lna(
    file, checksum = FALSE, qpp = "sample", qpp_sample_tiles = 2L,
    receipt = TRUE
  )
  full <- validate_lna(file, checksum = FALSE, qpp = "all", receipt = TRUE)
  expect_s3_class(full, "lna_validation_receipt")
  expect_identical(full, validate_lna(
    file, checksum = FALSE, qpp = "all", receipt = TRUE
  ))
  expect_identical(structural$qpp[[1L]]$tiles_checked, 0L)
  expect_identical(sampled$qpp[[1L]]$tiles_checked, 2L)
  expect_identical(full$qpp[[1L]]$tiles_checked,
                   full$qpp[[1L]]$tile_count)
  expect_true(full$qpp[[1L]]$q_roundtrip_verified)
  expect_length(full$qpp[[1L]]$parcelwise_RMSE, 4L)
  expect_named(full$qpp[[1L]]$checksums,
               c("receipt", "payload", "index", "steps", "settings"))
  expect_true(quick_validate_lna(file))
})

test_that("QPP index overlap, out-of-range entries, and truncation fail specifically", {
  overlap <- make_qpp_validation_file()
  on.exit(unlink(overlap), add = TRUE)
  mutate_qpp_index(overlap, function(index) {
    index$byte_offset[2] <- index$byte_offset[1]
    index
  })
  expect_qpp_validation_error(
    overlap, "validate_lna:qpp:index_ranges", "overlap|monotone",
    qpp = "structural"
  )

  outside <- make_qpp_validation_file()
  on.exit(unlink(outside), add = TRUE)
  mutate_qpp_index(outside, function(index) {
    index$byte_length[1] <- sum(index$byte_length) + 1
    index$uncompressed_length[1] <- index$byte_length[1]
    index
  })
  expect_qpp_validation_error(
    outside, "validate_lna:qpp:index_row", "byte range|deterministic",
    qpp = "structural"
  )

  truncated <- make_qpp_validation_file()
  on.exit(unlink(truncated), add = TRUE)
  mutate_qpp_index(truncated, function(index) {
    index$byte_length[nrow(index)] <- index$byte_length[nrow(index)] - 1
    index$uncompressed_length[nrow(index)] <-
      index$byte_length[nrow(index)]
    index
  })
  expect_qpp_validation_error(
    truncated, "validate_lna:qpp:index_ranges", "cover|payload",
    qpp = "structural"
  )
})

test_that("unsupported codecs and decompression bombs fail before tile decode", {
  codec <- make_qpp_validation_file()
  on.exit(unlink(codec), add = TRUE)
  mutate_qpp_index(codec, function(index) {
    index$outer_codec[1] <- 1
    index
  })
  expect_qpp_validation_error(
    codec, "validate_lna:qpp:outer_codec", "unsupported outer codec",
    qpp = "structural"
  )

  bomb <- make_qpp_validation_file()
  on.exit(unlink(bomb), add = TRUE)
  mutate_qpp_index(bomb, function(index) {
    index$uncompressed_length[1] <- 1000000
    index
  })
  expect_qpp_validation_error(
    bomb, "validate_lna:qpp:uncompressed_length", "exceeds its cap",
    qpp = "structural", qpp_limits = list(max_uncompressed_tile_bytes = 1024)
  )
})

test_that("bad predictor modes, widths, and anchors fail in bounded tile decode", {
  cases <- list(
    mode = list(byte = 9L, value = 9L, pattern = "predictor|header"),
    width = list(byte = 10L, value = 35L, pattern = "width|header"),
    anchors = list(
      byte = 11L, value = 0L, force_mode = 1L,
      pattern = "anchor|predictor lengths"
    )
  )
  for (case in cases) {
    file <- make_qpp_validation_file()
    on.exit(unlink(file), add = TRUE)
    mutate_qpp_payload(file, function(payload) {
      if (!is.null(case$force_mode)) {
        payload[[9L]] <- as.raw(case$force_mode)
      }
      payload[[case$byte]] <- as.raw(case$value)
      payload
    })
    refresh_qpp_tile_digests(file)
    expect_qpp_validation_error(
      file, "validate_lna:qpp:tile_decode", case$pattern
    )
  }
})

test_that("integer overflow and payload checksum changes are detected", {
  overflow <- make_qpp_validation_file(
    y = matrix(c(1000000, 999999.99), ncol = 1),
    time_chunk = 2L, parcel_chunk = 1L
  )
  on.exit(unlink(overflow), add = TRUE)
  overflow_frame <- c(
    charToRaw("QPP1"), .qpp_ref_u16le(2), .qpp_ref_u16le(1),
    as.raw(c(1, 2, 1, 0)), .qpp_ref_u16le(1), .qpp_ref_u32le(1),
    .qpp_ref_i32le(2147483647), .qpp_ref_pack(2, 2)
  )
  mutate_qpp_payload(overflow, function(payload) {
    expect_identical(length(payload), length(overflow_frame))
    overflow_frame
  })
  refresh_qpp_tile_digests(overflow)
  expect_qpp_validation_error(
    overflow, "validate_lna:qpp:tile_decode", "overflow"
  )

  changed <- make_qpp_validation_file()
  on.exit(unlink(changed), add = TRUE)
  mutate_qpp_payload(changed, function(payload) {
    payload[[length(payload)]] <- as.raw(bitwXor(as.integer(payload[[length(payload)]]), 1L))
    payload
  })
  error <- rlang::catch_cnd(validate_lna(changed, checksum = FALSE, qpp = "all"))
  expect_s3_class(error, "lna_error_qpp_integrity")
  expect_identical(error$code, "qpp_tile_digest_mismatch")
  expect_identical(error$location, "validate_lna:qpp:tile_digest")
})

test_that("receipt, steps, parcel maps, and calibration metadata are protected", {
  receipt_file <- make_qpp_validation_file()
  on.exit(unlink(receipt_file), add = TRUE)
  with_qpp_mutation(receipt_file, function(h5, root, group, name, desc) {
    path <- qpp_validation_path(desc, "qpp_quality_receipt")
    receipt <- jsonlite::fromJSON(h5_read(root, path), simplifyVector = FALSE)
    checksum <- receipt$receipt_checksum
    last <- substr(checksum, nchar(checksum), nchar(checksum))
    replacement <- if (identical(last, "0")) "1" else "0"
    receipt$receipt_checksum <- paste0(
      substr(checksum, 1L, nchar(checksum) - 1L), replacement
    )
    dataset <- h5[[path]]
    on.exit(safe_h5_close(dataset), add = TRUE)
    dataset$write(args = list(1L), value = as.character(jsonlite::toJSON(
      receipt, auto_unbox = TRUE, null = "null", digits = NA
    )))
  })
  expect_qpp_validation_error(
    receipt_file, "validate_lna:qpp:quality_receipt_checksum",
    "receipt checksum", qpp = "structural"
  )

  steps_file <- make_qpp_validation_file()
  on.exit(unlink(steps_file), add = TRUE)
  with_qpp_mutation(steps_file, function(h5, root, group, name, desc) {
    path <- qpp_validation_path(desc, "qpp_parcel_steps")
    dataset <- h5[[path]]
    on.exit(safe_h5_close(dataset), add = TRUE)
    values <- as.numeric(dataset$read())
    values[[1L]] <- values[[1L]] * 2
    dataset$write(args = list(seq_along(values)), value = values)
  })
  expect_qpp_validation_error(
    steps_file, "validate_lna:qpp:step_digest", "step",
    qpp = "structural", class = "lna_error_qpp_integrity"
  )

  map_file <- make_qpp_validation_file(with_parcel = TRUE)
  on.exit(unlink(map_file), add = TRUE)
  with_qpp_mutation(map_file, function(h5, root, group, name, desc) {
    path <- desc$params$parcel_binding$labels_path
    dataset <- h5[[path]]
    on.exit(safe_h5_close(dataset), add = TRUE)
    labels <- as.integer(dataset$read())
    labels[[1L]] <- 9L
    dataset$write(args = list(seq_along(labels)), value = labels)
  })
  expect_qpp_validation_error(
    map_file, "validate_lna:qpp:parcel_map_labels", "parcel map",
    qpp = "structural"
  )

  calibration <- make_qpp_validation_file(policy = "external-calibration")
  on.exit(unlink(calibration), add = TRUE)
  with_qpp_mutation(calibration, function(h5, root, group, name, desc) {
    desc$params$calibration_id <- ""
    write_json_descriptor(group, name, desc)
  })
  expect_qpp_validation_error(
    calibration, "analysis_contract:calibration_id", "calibration",
    qpp = "structural", class = "lna_error_analysis_contract"
  )
})

test_that("quick validation is bounded and full validation catches payload corruption", {
  file <- make_qpp_validation_file()
  on.exit(unlink(file), add = TRUE)
  mutate_qpp_payload(file, function(payload) {
    payload[[9L]] <- as.raw(9L)
    payload
  })
  refresh_qpp_tile_digests(file)
  expect_true(quick_validate_lna(file))
  expect_qpp_validation_error(
    file, "validate_lna:qpp:tile_decode", "predictor|header"
  )

  capped <- make_qpp_validation_file()
  on.exit(unlink(capped), add = TRUE)
  error <- rlang::catch_cnd(validate_lna(
    capped, checksum = FALSE,
    qpp_limits = list(max_tiles = 1)
  ))
  expect_s3_class(error, "lna_error_qpp_corrupt")
  expect_identical(error$code, "qpp_resource_limit")
  expect_identical(error$location, "validate_lna:qpp:tile_count")
})

test_that("unknown required decoders are never silently skipped", {
  file <- make_qpp_validation_file()
  on.exit(unlink(file), add = TRUE)
  with_qpp_mutation(file, function(h5, root, group, name, desc) {
    current <- as.character(h5_attr_read(root, "required_transforms"))
    h5_attr_write(root, "required_transforms", c(current, "unknown@1.0.0"))
  })
  error <- rlang::catch_cnd(validate_lna(file, checksum = FALSE))
  expect_s3_class(error, "lna_error_required_transform")
  expect_identical(error$location, "validate_lna:required_transforms")
  expect_match(conditionMessage(error), "unavailable")
  expect_error(quick_validate_lna(file), class = "lna_error_required_transform")
})

test_that("analysis-safe qpp binding rejects missing and contradictory provenance", {
  missing <- make_qpp_validation_file(with_parcel = TRUE)
  on.exit(unlink(missing), add = TRUE)
  mutate_qpp_descriptor(missing, function(desc) {
    desc$params$parcel_binding <- NULL
    desc
  })
  error <- rlang::catch_cnd(validate_lna(missing, checksum = FALSE))
  expect_s3_class(error, "lna_error_qpp_binding")
  expect_identical(error$code, "qpp_parcel_binding")
  expect_identical(error$location, "analysis_contract:qpp_parcel_binding")

  wrong_run <- make_qpp_validation_file(with_parcel = TRUE)
  on.exit(unlink(wrong_run), add = TRUE)
  mutate_qpp_descriptor(wrong_run, function(desc) {
    desc$params$parcel_binding$run_id <- "run-02"
    desc
  })
  error <- rlang::catch_cnd(validate_lna(wrong_run, checksum = FALSE))
  expect_s3_class(error, "lna_error_qpp_binding")
  expect_identical(error$code, "qpp_parcel_wrong_run")
  expect_identical(error$location, "validate_lna:qpp:parcel_binding")

  wrong_step <- make_qpp_validation_file(with_parcel = TRUE)
  on.exit(unlink(wrong_step), add = TRUE)
  mutate_qpp_descriptor(wrong_step, function(desc) {
    desc$params$parcel_binding$step_index <-
      desc$params$parcel_binding$step_index + 1L
    desc
  })
  error <- rlang::catch_cnd(validate_lna(wrong_step, checksum = FALSE))
  expect_s3_class(error, "lna_error_qpp_binding")
  expect_identical(error$code, "qpp_parcel_missing")

  wrong_map <- make_qpp_validation_file(with_parcel = TRUE)
  on.exit(unlink(wrong_map), add = TRUE)
  mutate_qpp_descriptor(wrong_map, function(desc) {
    desc$params$parcel_binding$map_checksum <-
      paste0("sha256:", strrep("0", 64))
    desc
  })
  error <- rlang::catch_cnd(validate_lna(wrong_map, checksum = FALSE))
  expect_s3_class(error, "lna_error_qpp_binding")
  expect_identical(error$code, "qpp_parcel_missing")
})

test_that("analysis-safe qpp rejects duplicate bindings but ignores link order", {
  reordered <- make_qpp_validation_file(with_parcel = TRUE)
  on.exit(unlink(reordered), add = TRUE)
  mutate_qpp_descriptor(reordered, identity)
  expect_true(validate_lna(reordered, checksum = FALSE, qpp = "structural"))

  duplicate <- make_qpp_validation_file(with_parcel = TRUE)
  on.exit(unlink(duplicate), add = TRUE)
  with_qpp_mutation(duplicate, function(h5, root, group, name, desc) {
    parcel <- read_json_descriptor(group, desc$params$parcel_binding$descriptor)
    write_json_descriptor(group, "02_parcel.mean.json", parcel)
  })
  error <- rlang::catch_cnd(validate_lna(duplicate, checksum = FALSE))
  expect_s3_class(error, "lna_error_qpp_binding")
  expect_identical(error$code, "qpp_parcel_duplicate")
  expect_identical(error$location, "validate_lna:qpp:parcel_binding")
})

test_that("standalone qpp cannot be relabelled analysis-safe", {
  file <- make_qpp_validation_file(with_parcel = FALSE)
  on.exit(unlink(file), add = TRUE)
  mutate_qpp_descriptor(file, function(desc) {
    desc$profile_id <- "analysis-safe"
    desc$params$profile_id <- "analysis-safe"
    desc
  })
  error <- rlang::catch_cnd(validate_lna(file, checksum = FALSE))
  expect_s3_class(error, "lna_error_qpp_binding")
  expect_identical(error$code, "qpp_parcel_binding")
  expect_identical(error$location, "analysis_contract:qpp_parcel_binding")
})

test_that("unknown and legacy qpp transform versions fail explicitly", {
  for (version in c("1.0.0", "9.9.9")) {
    file <- make_qpp_validation_file()
    on.exit(unlink(file), add = TRUE)
    mutate_qpp_descriptor(file, function(desc) {
      desc$version <- version
      desc
    })
    error <- rlang::catch_cnd(validate_lna(file, checksum = FALSE))
    expect_s3_class(error, "lna_error_required_transform")
    expect_identical(error$location, "validate_lna:qpp:decoder_version")
    expect_identical(
      error$code,
      if (identical(version, "1.0.0")) "qpp_legacy_profile" else "required_transform"
    )
    quick_error <- rlang::catch_cnd(quick_validate_lna(file))
    expect_s3_class(quick_error, "lna_error_required_transform")
    expect_identical(quick_error$code, "required_transform")
    expect_identical(quick_error$location, "validate_lna:required_transforms")
  }
})

test_that("legacy qpp is decode-only when named validation is disabled", {
  file <- make_qpp_validation_file()
  on.exit(unlink(file), add = TRUE)
  current <- read_lna(
    file, validate = FALSE, output_dtype = "float64", time_idx = 1:4
  )
  mutate_qpp_descriptor(file, function(desc) {
    desc$version <- "1.0.0"
    desc$codec_version <- NULL
    desc
  })

  expect_warning(
    legacy <- read_lna(
      file, validate = FALSE, output_dtype = "float64", time_idx = 1:4
    ),
    class = "lna_warning_qpp_legacy_integrity"
  )
  expect_equal(legacy$stash$input, current$stash$input, tolerance = 1e-12)
  expect_identical(legacy$meta$qpp_io_receipt$tile_digest_rows_read, 0L)
  expect_identical(legacy$meta$qpp_io_receipt$step_digest_rows_read, 0L)
  expect_false(legacy$meta$qpp_io_receipt$selected_tile_digests_verified)
  expect_false(legacy$meta$qpp_io_receipt$selected_step_digests_verified)
})

test_that("normal reads reject selected payload and parcel-step mutations", {
  payload_file <- make_qpp_validation_file()
  on.exit(unlink(payload_file), add = TRUE)
  mutate_qpp_payload(payload_file, function(payload) {
    payload[[1L]] <- as.raw(bitwXor(as.integer(payload[[1L]]), 1L))
    payload
  })
  error <- rlang::catch_cnd(read_lna(
    payload_file, time_idx = 1L, validate = FALSE
  ))
  expect_s3_class(error, "lna_error_qpp_integrity")
  expect_identical(error$code, "qpp_tile_digest_mismatch")
  expect_identical(error$location, "invert_step.qpp:tile_digest")

  step_file <- make_qpp_validation_file()
  on.exit(unlink(step_file), add = TRUE)
  with_qpp_mutation(step_file, function(h5, root, group, name, desc) {
    path <- qpp_validation_path(desc, "qpp_parcel_steps")
    dataset <- h5[[path]]
    on.exit(safe_h5_close(dataset), add = TRUE)
    values <- as.numeric(dataset$read())
    values[[1L]] <- values[[1L]] * 2
    dataset$write(args = list(seq_along(values)), value = values)
  })
  error <- rlang::catch_cnd(read_lna(
    step_file, time_idx = 1L, validate = FALSE
  ))
  expect_s3_class(error, "lna_error_qpp_integrity")
  expect_identical(error$code, "qpp_step_digest_mismatch")
  expect_identical(error$location, "invert_step.qpp:step_digest")
})

test_that("unselected corruption remains bounded while full validation detects it", {
  file <- make_qpp_validation_file(time_chunk = 4L, parcel_chunk = 2L)
  on.exit(unlink(file), add = TRUE)
  mutate_qpp_payload(file, function(payload) {
    payload[[length(payload)]] <- as.raw(
      bitwXor(as.integer(payload[[length(payload)]]), 1L)
    )
    payload
  })
  selected <- read_lna(
    file, time_idx = 1L, validate = FALSE, output_dtype = "float64"
  )
  expect_identical(selected$meta$qpp_io_receipt$index_rows_read, 2L)
  error <- rlang::catch_cnd(validate_lna(file, checksum = FALSE, qpp = "all"))
  expect_s3_class(error, "lna_error_qpp_integrity")
  expect_identical(error$code, "qpp_tile_digest_mismatch")
})
