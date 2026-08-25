qpp_remediation_normalization <- function() {
  list(
    method = "none", fit_scope = "fixed", source_id = "input-values-v1",
    calibration_identity = NULL, fit_id = NULL
  )
}

qpp_remediation_fixture <- function() {
  mask <- array(TRUE, dim = c(8, 1, 1))
  labels <- array(rep(1:4, each = 2), dim = dim(mask))
  affine <- diag(c(2, 2, 2, 1))
  header <- list(
    dims = dim(mask), spacing = c(2, 2, 2), origin = c(0, 0, 0),
    transform = affine, spatial_units = "mm",
    voxel_ordering = "first-axis-fastest"
  )
  base <- outer(
    seq_len(8), seq_len(12),
    function(voxel, time) voxel + sin(time / 3)
  )
  runs <- list(
    "run-01" = array(base, dim = c(8, 1, 1, 12)),
    "run-02" = array(base + outer(seq_len(8), seq_len(12), function(v, t) {
      10 + v / 10 + cos(t / 4)
    }), dim = c(8, 1, 1, 12))
  )
  profile <- lna_profile_analysis_safe(
    labels = labels, delta = c(0.01, 0.02, 0.03, 0.04),
    signal_units = "percent-signal-change", source_id = "atlas-v1",
    normalization_provenance = qpp_remediation_normalization(),
    map_grid_id = .lna_grid_id(
      dim(mask), affine, "mm", "first-axis-fastest"
    ),
    time_chunk = 4L, parcel_chunk = 2L
  )
  list(
    mask = mask, labels = labels, affine = affine, header = header,
    runs = runs, profile = profile
  )
}

write_qpp_remediation_fixture <- function(path, fixture = qpp_remediation_fixture()) {
  write_lna(
    fixture$runs, path,
    transforms = fixture$profile$transforms,
    transform_params = fixture$profile$transform_params,
    mask = fixture$mask,
    run_geometry = list(
      "run-01" = fixture$header,
      "run-02" = fixture$header
    ),
    checksum = "none"
  )
  invisible(fixture)
}

read_qpp_remediation_descriptors <- function(file) {
  h5 <- neuroarchive:::open_h5(file, "r")
  on.exit(neuroarchive:::close_h5_safely(h5), add = TRUE)
  group <- h5[["transforms"]]
  on.exit(safe_h5_close(group), add = TRUE)
  discovered <- discover_transforms(group)
  descriptors <- lapply(discovered$name, function(name) {
    read_json_descriptor(group, name)
  })
  names(descriptors) <- discovered$name
  descriptors
}

test_that("public two-run profile records exact run-matched parcel bindings", {
  file <- local_tempfile(fileext = ".h5")
  fixture <- write_qpp_remediation_fixture(file)
  expect_true(quick_validate_lna(file))
  receipt <- validate_lna(file, checksum = FALSE, qpp = "all", receipt = TRUE)
  expect_identical(
    receipt$required_transforms,
    c("parcel.mean@1.1.0", "qpp@1.1.0")
  )
  expect_length(receipt$qpp, 2L)

  descriptors <- read_qpp_remediation_descriptors(file)
  parcels <- descriptors[vapply(
    descriptors, function(desc) identical(desc$type, "parcel.mean"), logical(1)
  )]
  qpps <- descriptors[vapply(
    descriptors, function(desc) identical(desc$type, "qpp"), logical(1)
  )]
  expect_length(parcels, 2L)
  expect_length(qpps, 2L)
  expect_identical(
    unname(vapply(parcels, function(desc) desc$run_id, character(1))),
    c("run-01", "run-02")
  )
  expect_identical(
    unname(vapply(qpps, function(desc) desc$run_id, character(1))),
    c("run-01", "run-02")
  )
  for (qpp_name in names(qpps)) {
    qpp <- qpps[[qpp_name]]
    binding <- qpp$params$parcel_binding
    expect_identical(binding$run_id, qpp$run_id)
    expect_true(binding$descriptor %in% names(parcels))
    parcel <- parcels[[binding$descriptor]]
    expect_identical(binding$transform, "parcel.mean@1.1.0")
    expect_identical(binding$step_index, parcel$params$step_index)
    expect_identical(binding$asset_id, parcel$params$asset_id)
    expect_identical(binding$map_checksum, parcel$params$map_checksum)
    expect_identical(binding$mask_checksum, parcel$params$mask_checksum)
    expect_identical(
      binding$labels_path,
      Filter(function(x) x$role == "parcel_labels", parcel$datasets)[[1L]]$path
    )
  }
  expect_identical(
    length(unique(vapply(parcels, function(desc) desc$params$asset_id, ""))),
    1L
  )

  eager <- read_lna(
    file, run_id = "run-*", validate = FALSE, output_dtype = "float64"
  )
  expect_identical(names(eager), c("run-01", "run-02"))
  expect_false(identical(eager[[1L]]$stash$input, eager[[2L]]$stash$input))

  lazy <- read_lna(
    file, run_id = "run-*", validate = FALSE, output_dtype = "float64",
    lazy = TRUE
  )
  on.exit(lapply(lazy, function(reader) reader$close()), add = TRUE)
  expect_identical(names(lazy), c("run-01", "run-02"))
  expect_true(all(vapply(lazy, inherits, logical(1), "lna_reader")))
  expect_identical(
    vapply(lazy, function(reader) reader$current_run_id, character(1)),
    c("run-01" = "run-01", "run-02" = "run-02")
  )
  for (rid in names(lazy)) {
    multi_subset <- lazy[[rid]]$data(time_idx = 3:8)
    expect_equal(
      multi_subset$stash$input,
      eager[[rid]]$stash$input[3:8, , drop = FALSE],
      tolerance = 1e-12
    )

    single_file <- local_tempfile(fileext = ".h5")
    write_lna(
      fixture$runs[[rid]], single_file,
      transforms = fixture$profile$transforms,
      transform_params = fixture$profile$transform_params,
      mask = fixture$mask, header = fixture$header, run_id = rid,
      checksum = "none"
    )
    single_subset <- read_lna(
      single_file, run_id = rid, validate = FALSE, output_dtype = "float64",
      time_idx = 3:8
    )
    expect_equal(
      multi_subset$stash$input,
      single_subset$stash$input,
      tolerance = 1e-12
    )
    expect_identical(
      multi_subset$meta$qpp_decoded_q_checksum,
      single_subset$meta$qpp_decoded_q_checksum
    )

    roi <- array(FALSE, dim = dim(fixture$mask))
    roi[c(1L, 8L)] <- TRUE
    multi_roi <- lazy[[rid]]$data(roi_mask = roi, time_idx = c(12L, 1L, 5L))
    single_roi <- read_lna(
      single_file, run_id = rid, validate = FALSE, output_dtype = "float64",
      roi_mask = roi, time_idx = c(12L, 1L, 5L)
    )
    expect_equal(multi_roi$stash$input, single_roi$stash$input, tolerance = 1e-12)
    expect_identical(multi_roi$meta$qpp_selected_parcels, c(1L, 4L))
    expect_identical(multi_roi$meta$qpp_io_receipt$label_values_read, 2L)

    single_receipt <- validate_lna(
      single_file, checksum = FALSE, qpp = "all", receipt = TRUE
    )
    multi_qpp_name <- names(qpps)[vapply(
      qpps, function(desc) identical(desc$run_id, rid), logical(1)
    )][[1L]]
    multi_validation <- receipt$qpp[[which(vapply(
      receipt$qpp,
      function(entry) identical(entry$descriptor, multi_qpp_name),
      logical(1)
    ))]]
    expect_true(multi_validation$q_roundtrip_verified)
    expect_true(single_receipt$qpp[[1L]]$q_roundtrip_verified)
    expect_identical(
      multi_validation$checksums[c("payload", "steps")],
      single_receipt$qpp[[1L]]$checksums[c("payload", "steps")]
    )
    expect_identical(
      multi_validation$parcelwise_RMSE,
      single_receipt$qpp[[1L]]$parcelwise_RMSE
    )
  }
  lazy[[1L]]$close()
  expect_null(lazy[[1L]]$h5)
  expect_true(lazy[[2L]]$h5$is_valid)

  h5 <- neuroarchive:::open_h5(file, "r")
  on.exit(neuroarchive:::close_h5_safely(h5), add = TRUE)
  assets <- h5[["spatial/qpp_assets"]]
  on.exit(safe_h5_close(assets), add = TRUE)
  expect_length(names(assets), 1L)
  for (i in seq_along(qpps)) {
    index_path <- Filter(
      function(x) x$role == "qpp_tile_index", qpps[[i]]$datasets
    )[[1L]]$path
    index <- h5_read(h5[["/"]], index_path)
    expect_true(all(index$run_ordinal == i - 1L))
    expect_true(all(index$unit_ordinal == i - 1L))
  }
})

test_that("cross-run parcel substitution fails in quick and full validation", {
  file <- local_tempfile(fileext = ".h5")
  write_qpp_remediation_fixture(file)
  h5 <- neuroarchive:::open_h5(file, "r+")
  group <- h5[["transforms"]]
  run_one <- read_json_descriptor(group, "01_qpp.json")
  run_two <- read_json_descriptor(group, "03_qpp.json")
  run_one$params$parcel_binding <- run_two$params$parcel_binding
  write_json_descriptor(group, "01_qpp.json", run_one)
  h5$flush()
  safe_h5_close(group)
  neuroarchive:::close_h5_safely(h5)

  quick_error <- rlang::catch_cnd(quick_validate_lna(file))
  expect_s3_class(quick_error, "lna_error_qpp_binding")
  expect_identical(quick_error$code, "qpp_parcel_wrong_run")
  expect_identical(quick_error$location, "validate_lna:qpp:parcel_binding")

  full_error <- rlang::catch_cnd(validate_lna(file, checksum = FALSE))
  expect_s3_class(full_error, "lna_error_qpp_binding")
  expect_identical(full_error$code, "qpp_parcel_wrong_run")
  expect_identical(full_error$location, "validate_lna:qpp:parcel_binding")
})

test_that("identity-addressed parcel assets deduplicate or reject collisions", {
  plan <- Plan$new("run-01")
  args <- list(
    plan = plan,
    path = "/spatial/qpp_assets/fixture/parcel_labels",
    value = c(1L, 1L, 2L, 2L),
    dtype = "uint32", role = "parcel_labels", run_id = "run-01",
    step_index = 0L, params_json = "{}"
  )
  expect_invisible(do.call(.parcel_mean_add_asset, args))
  expect_invisible(do.call(.parcel_mean_add_asset, args))
  expect_length(plan$payloads, 1L)

  args$value <- c(1L, 2L, 1L, 2L)
  collision <- rlang::catch_cnd(do.call(.parcel_mean_add_asset, args))
  expect_s3_class(collision, "lna_error_spatial_identity")
  expect_identical(collision$code, "parcel_asset_collision")
  expect_identical(collision$location, ".parcel_mean_add_asset")
})
