qpp_transform_handle <- function(y, units = "percent-signal-change") {
  DataHandle$new(
    initial_stash = list(parcel_series = as.matrix(y)),
    initial_meta = list(parcel_mean_assets = list(signal_units = units)),
    plan = Plan$new(), run_ids = "run-01", current_run_id = "run-01"
  )
}

run_qpp_transform <- function(y, params, output = "qpp_packed") {
  handle <- qpp_transform_handle(y)
  desc <- list(type = "qpp", params = params, inputs = list("parcel_series"),
               outputs = list(output))
  run_transform_step("forward", "qpp", desc, handle, 0L)
}

materialize_qpp_plan <- function(out) {
  path <- withr::local_tempfile(fileext = ".h5", .local_envir = parent.frame())
  h5 <- neuroarchive:::open_h5(path, "w")
  materialise_plan(h5, out$plan, run_ids = "run-01")
  list(path = path, h5 = h5)
}

test_that("qpp fixed scalar delta stores only terminal representation", {
  y <- cbind(c(0, -1, 2, -3, 4), c(1000, -1000, 0, 0.5, -0.5))
  out <- run_qpp_transform(y, list(
    delta = 0.25, quantizer_policy = "fixed",
    units = "percent-signal-change", time_chunk = 3L, parcel_chunk = 1L
  ))
  roles <- unlist(out$plan$datasets$role, use.names = FALSE)
  expect_identical(roles, c("qpp_payload", "qpp_tile_index",
                            "qpp_parcel_steps", "qpp_tile_sha256",
                            "qpp_parcel_step_sha256",
                            "qpp_quality_receipt"))
  expect_false(any(grepl("uncompressed|parcel_series|quantized", roles)))
  expect_false(out$has_key("parcel_series"))
  expect_true(out$has_key("qpp_packed"))
  expect_identical(out$meta$qpp_receipt$fit_scope, "fixed")
  expect_identical(out$meta$qpp_receipt$tile_count, 4L)
  expect_equal(out$plan$payloads[[grep("parcel_steps", names(out$plan$payloads),
                                      value = TRUE)]], c(0.25, 0.25))
})

test_that("qpp supports per-parcel external calibration and exact q decode", {
  y <- cbind(c(-1.2, 0, 1.2, 2.4), c(-10, 5, 0, 10))
  delta <- c(0.1, 2.5)
  out <- run_qpp_transform(y, list(
    delta = delta,
    quantizer_policy = "external-calibration",
    units = "percent-signal-change",
    calibration_id = "cal-set-v1",
    calibration_checksum = paste0("sha256:", strrep("a", 64)),
    time_chunk = 64L, parcel_chunk = 64L
  ))
  stored <- materialize_qpp_plan(out)
  root <- stored$h5[["/"]]
  expect_identical(h5_attr_read(root, "required_transforms"), "qpp@1.1.0")
  desc <- read_json_descriptor(stored$h5[["transforms"]], "00_qpp.json")
  expect_identical(desc$profile_id, "parcel-series-qpp")
  expect_null(desc$params$parcel_binding)
  expect_identical(desc$params$fit_scope, "external-calibration")
  payload_path <- Filter(function(x) x$role == "qpp_payload", desc$datasets)[[1]]$path
  payload <- h5_read(root, payload_path)
  expected_q <- round(sweep(y, 2L, delta, "/"))
  expect_identical(qpp_decode_tile(.qpp_payload_as_raw(payload)), expected_q)

  read_handle <- DataHandle$new(h5 = stored$h5, current_run_id = "run-01")
  decoded <- run_transform_step("invert", "qpp", desc, read_handle, 0L)
  expect_equal(decoded$stash$parcel_series, sweep(expected_q, 2L, delta, "*"))
  safe_h5_close(root)
  neuroarchive:::close_h5_safely(stored$h5)
  expect_true(validate_lna(stored$path))
})

test_that("analysis-safe qpp rejects fitted or ambiguous parameter policies", {
  y <- matrix(c(-1, 0, 1), ncol = 1)
  for (policy in c("range", "sd", "mad", "maximum", "subject")) {
    expect_error(run_qpp_transform(y, list(
      delta = 0.1, quantizer_policy = policy,
      units = "percent-signal-change"
    )), "unsupported qpp quantizer_policy", info = policy)
  }
  expect_error(run_qpp_transform(y, list(
    delta = 0.1, quantizer_policy = "external-calibration",
    units = "percent-signal-change"
  )), "calibration")
  expect_error(run_qpp_transform(y, list(
    delta = 0.1, quantizer_policy = "fixed",
    units = "percent-signal-change", allow_clip = TRUE
  )), "never permits clipping")
  expect_error(run_qpp_transform(y, list(
    delta = 0.1, quantizer_policy = "fixed", units = "raw"
  )), "do not match")

  bad <- rlang::catch_cnd(run_qpp_transform(y, list(
    delta = 0.1, quantizer_policy = "external-calibration",
    units = "percent-signal-change", calibration_id = "external-v1",
    calibration_checksum = "not-a-sha256"
  )))
  expect_s3_class(bad, "lna_error_checksum")
  expect_identical(bad$code, "noncanonical_sha256")
  expect_identical(bad$location, "forward_step.qpp:calibration_checksum")
})

test_that("fitted-run-maxabs is explicit and cannot masquerade as analysis-safe", {
  y <- cbind(c(-10, 0, 5), c(0, 0, 0))
  params <- list(
    quantizer_policy = "fitted-run-maxabs",
    units = "percent-signal-change", fit_id = "run-01-fit",
    target_q = 100L, zero_step = 0.5
  )
  out <- run_qpp_transform(y, params)
  desc <- out$plan$descriptors[[1L]]
  expect_identical(desc$profile_id, "fitted-run")
  expect_identical(desc$params$fit_scope, "run")
  expect_identical(desc$params$quantizer_policy, "fitted-run-maxabs")
  expect_equal(out$plan$payloads[[grep("parcel_steps", names(out$plan$payloads),
                                      value = TRUE)]], c(0.1, 0.5))
  params$delta <- 0.1
  expect_error(run_qpp_transform(y, params), "forbids an explicit delta")
  params$delta <- NULL
  params$fit_id <- NULL
  expect_error(run_qpp_transform(y, params), "fit_id")
})

test_that("qpp detects invalid deltas and signed-int32 overflow", {
  y <- matrix(c(0, -1, 1, 1e12), ncol = 1)
  for (delta in list(NULL, 0, -1, NA_real_, c(1, 2))) {
    expect_error(run_qpp_transform(y, list(
      delta = delta, quantizer_policy = "fixed",
      units = "percent-signal-change"
    )), "delta")
  }
  expect_error(run_qpp_transform(y, list(
    delta = 1e-12, quantizer_policy = "fixed",
    units = "percent-signal-change"
  )), "signed-int32")
})

test_that("qpp tiles and run instances reset anchors independently", {
  first_y <- matrix(0:8, ncol = 1)
  second_y <- matrix(100:108, ncol = 1)
  params <- list(delta = 1, quantizer_policy = "fixed",
                 units = "percent-signal-change", time_chunk = 4L,
                 parcel_chunk = 1L)
  first <- run_qpp_transform(first_y, params)
  second <- run_qpp_transform(second_y, params)
  first_index <- first$plan$payloads[[grep("tile_index", names(first$plan$payloads),
                                           value = TRUE)]]
  second_index <- second$plan$payloads[[grep("tile_index", names(second$plan$payloads),
                                             value = TRUE)]]
  expect_identical(first_index$time_start, c(0L, 4L, 8L))
  expect_identical(second_index$time_start, c(0L, 4L, 8L))
  expect_true(all(first_index$unit_ordinal == 0))
  expect_true(all(second_index$unit_ordinal == 0))
  expect_false(identical(first$stash$qpp_packed, second$stash$qpp_packed))

  list_handle <- qpp_transform_handle(first_y)
  list_handle$stash$parcel_series <- list(run1 = first_y, run2 = second_y)
  desc <- list(type = "qpp", params = params, inputs = list("parcel_series"),
               outputs = list("qpp_packed"))
  expect_error(run_transform_step("forward", "qpp", desc, list_handle, 0L),
               "one run")
})

test_that("parcel.mean to qpp archive roundtrips end to end without intermediates", {
  mask <- array(c(TRUE, FALSE, TRUE, TRUE), dim = c(2, 2, 1))
  x <- array(seq_len(20), dim = c(2, 2, 1, 5))
  labels <- array(c(10, 0, 30, 10), dim = dim(mask))
  affine <- diag(c(2, 2, 2, 1))
  header <- list(
    dims = dim(mask), spacing = c(2, 2, 2), origin = c(0, 0, 0),
    transform = affine, spatial_units = "mm",
    voxel_ordering = "first-axis-fastest"
  )
  parcel_params <- list(
    labels = labels, fit_scope = "fixed", source_id = "fixture-atlas-v1",
    signal_units = "percent-signal-change", spatial_units = "mm",
    voxel_ordering = "first-axis-fastest",
    map_grid_id = .lna_grid_id(
      dim(mask), affine, "mm", "first-axis-fastest"
    )
  )
  qpp_params <- list(
    delta = c(0.5, 1), quantizer_policy = "fixed",
    units = "percent-signal-change", time_chunk = 3L, parcel_chunk = 2L
  )
  result <- core_write(
    x, transforms = c("parcel.mean", "qpp"),
    transform_params = list("parcel.mean" = parcel_params, qpp = qpp_params),
    mask = mask, header = header
  )
  roles <- unlist(result$plan$datasets$role, use.names = FALSE)
  expect_false(any(roles %in% c("parcel_series", "quantized", "q")))
  expect_true(all(c("parcel_labels", "qpp_payload", "qpp_tile_index") %in% roles))

  tmp <- local_tempfile(fileext = ".h5")
  h5 <- neuroarchive:::open_h5(tmp, "w")
  materialise_plan(
    h5, result$plan, run_ids = "run-01",
    header = header
  )
  root <- h5[["/"]]
  expect_identical(h5_attr_read(root, "required_transforms"),
                   c("parcel.mean@1.1.0", "qpp@1.1.0"))
  safe_h5_close(root)
  neuroarchive:::close_h5_safely(h5)
  expect_true(validate_lna(tmp))

  decoded <- core_read(tmp, validate = TRUE, output_dtype = "float64")
  vox_time <- convert_to_masked_vox_time_matrix(x, mask)
  parcel <- cbind(colMeans(vox_time[c(1, 3), , drop = FALSE]),
                  vox_time[2, ])
  q <- round(sweep(parcel, 2L, qpp_params$delta, "/"))
  reconstructed_parcel <- sweep(q, 2L, qpp_params$delta, "*")
  expected <- reconstructed_parcel[, c(1, 2, 1), drop = FALSE]
  expect_equal(decoded$stash$input, expected)
})

test_that("identical qpp inputs produce identical payload index steps and receipt", {
  y <- cbind(sin(seq_len(17)), cos(seq_len(17)))
  params <- list(delta = c(0.01, 0.02), quantizer_policy = "fixed",
                 units = "percent-signal-change", time_chunk = 5L,
                 parcel_chunk = 1L)
  first <- run_qpp_transform(y, params)
  second <- run_qpp_transform(y, params)
  qpp_payloads <- function(out) out$plan$payloads[
    vapply(out$plan$datasets$role, function(role) {
      role %in% c("qpp_payload", "qpp_tile_index", "qpp_parcel_steps",
                  "qpp_tile_sha256", "qpp_parcel_step_sha256",
                  "qpp_quality_receipt")
    }, logical(1))
  ]
  expect_identical(qpp_payloads(first), qpp_payloads(second))
})

test_that("qpp canonicalizes signed zero before recording q evidence", {
  y <- matrix(c(-0.001, 0.001, -0.004, 0.004), ncol = 1L)
  q <- .qpp_quantize_pointwise(y, 0.01)
  expect_identical(q, matrix(0, nrow(y), 1L))
  decoded <- qpp_decode_tile(qpp_encode_tile(q))
  expect_identical(
    .parcel_mean_numeric_checksum(q, "lna_qpp_tile_q_v1"),
    .parcel_mean_numeric_checksum(decoded, "lna_qpp_tile_q_v1")
  )
})

test_that("qpp schema is valid JSON", {
  schema <- test_inst_file("schemas", "qpp.schema.json")
  expect_true(jsonlite::validate(paste(readLines(schema, warn = FALSE),
                                       collapse = "\n")))
})
