qpp_fixed_normalization <- function() {
  list(
    method = "none", fit_scope = "fixed", source_id = "input-values-v1",
    calibration_identity = NULL, fit_id = NULL
  )
}

qpp_profile_fixture <- function() {
  mask <- array(TRUE, dim = c(8, 1, 1))
  labels <- array(rep(1:4, each = 2), dim = dim(mask))
  x <- array(
    outer(seq_len(8), seq_len(12), function(v, t) v + sin(t / 3)),
    dim = c(8, 1, 1, 12)
  )
  affine <- diag(c(2, 2, 2, 1))
  list(
    mask = mask,
    labels = labels,
    x = x,
    affine = affine,
    grid_id = .lna_grid_id(
      dim(mask), affine, "mm", "first-axis-fastest"
    )
  )
}

test_that("analysis-safe profile produces explicit experimental write arguments", {
  fixture <- qpp_profile_fixture()
  profile <- lna_profile_analysis_safe(
    labels = fixture$labels, delta = 0.01,
    signal_units = "percent-signal-change", source_id = "atlas-v1",
    normalization_provenance = qpp_fixed_normalization(),
    map_grid_id = fixture$grid_id,
    time_chunk = 4L, parcel_chunk = 2L
  )

  expect_s3_class(profile, "lna_profile")
  expect_identical(profile$profile_id, "analysis-safe")
  expect_identical(profile$lifecycle, "experimental")
  expect_true(profile$available)
  expect_identical(profile$allowed_split_unit, "pointwise")
  expect_identical(
    profile$allowed_split_description,
    "individual TR observations may be assigned independently"
  )
  expect_identical(profile$transforms, c("parcel.mean", "qpp"))
  expect_identical(profile$transform_params$qpp$quantizer_policy, "fixed")
  expect_identical(profile$transform_params$qpp$time_chunk, 4L)
  expect_identical(
    profile$transform_params[["parcel.mean"]]$normalization_provenance,
    qpp_fixed_normalization()
  )
  expect_output(print(profile), "LNA profile: analysis-safe", fixed = TRUE)
  expect_output(print(profile), "lifecycle: experimental", fixed = TRUE)
  expect_output(
    print(profile), "allowed split unit (machine): pointwise", fixed = TRUE
  )
  expect_output(print(profile), "individual TR observations", fixed = TRUE)
})

test_that("analysis-safe profile fails closed on provenance", {
  fixture <- qpp_profile_fixture()
  common <- list(
    labels = fixture$labels, delta = 0.01, signal_units = "psc",
    source_id = "atlas-v1", map_grid_id = fixture$grid_id
  )
  expect_error(
    do.call(lna_profile_analysis_safe, common),
    "normalization_provenance is required",
    class = "lna_error_analysis_contract"
  )
  without_grid <- common
  without_grid$map_grid_id <- NULL
  expect_error(
    do.call(lna_profile_analysis_safe, c(without_grid, list(
      normalization_provenance = qpp_fixed_normalization()
    ))),
    "map_grid_id",
    class = "lna_error_spatial_identity"
  )
  expect_error(
    do.call(lna_profile_analysis_safe, c(common, list(
      normalization_provenance = qpp_fixed_normalization(),
      parcellation_fit_scope = "subject"
    ))),
    "should be one of"
  )
  run_normalized <- qpp_fixed_normalization()
  run_normalized$method <- "percent-signal-change"
  run_normalized$fit_scope <- "run"
  run_normalized$source_id <- "whole-run-mean"
  run_normalized$fit_id <- "run-01-mean-fit"
  expect_error(
    do.call(lna_profile_analysis_safe, c(common, list(
      normalization_provenance = run_normalized
    ))),
    "normalization fit_scope must be fixed or external-calibration",
    class = "lna_error_analysis_contract"
  )
  expect_error(
    do.call(lna_profile_analysis_safe, c(common, list(
      normalization_provenance = qpp_fixed_normalization(),
      quantizer_policy = "external-calibration"
    ))),
    "quantizer_calibration_id"
  )

  sha <- paste0("sha256:", strrep("A", 64))
  external <- do.call(lna_profile_analysis_safe, c(common, list(
    normalization_provenance = qpp_fixed_normalization(),
    parcellation_fit_scope = "external-calibration",
    parcellation_calibration_id = "atlas-calibration-v1",
    parcellation_calibration_checksum = sha,
    quantizer_policy = "external-calibration",
    quantizer_calibration_id = "step-calibration-v1",
    quantizer_calibration_checksum = sha
  )))
  expect_identical(
    external$transform_params[["parcel.mean"]]$calibration_id,
    "atlas-calibration-v1"
  )
  expect_identical(
    external$transform_params$qpp$calibration_id,
    "step-calibration-v1"
  )
  expect_identical(
    external$transform_params$qpp$calibration_checksum, tolower(sha)
  )

  bad_normalization <- list(
    method = "external-psc", fit_scope = "external-calibration",
    source_id = "normalization-v1",
    calibration_identity = list(
      id = "normalization-cal-v1", checksum = "not-a-sha256",
      relationship = "external-disjoint"
    ),
    fit_id = NULL
  )
  error <- rlang::catch_cnd(do.call(
    lna_profile_analysis_safe,
    c(common, list(normalization_provenance = bad_normalization))
  ))
  expect_s3_class(error, "lna_error_checksum")
  expect_identical(error$code, "noncanonical_sha256")

  contradictory <- qpp_fixed_normalization()
  contradictory$calibration_identity <- list(
    id = "unexpected-calibration",
    checksum = paste0("sha256:", strrep("f", 64)),
    relationship = "external-disjoint"
  )
  expect_error(
    do.call(lna_profile_analysis_safe, c(common, list(
      normalization_provenance = contradictory
    ))),
    "calibration_identity requires external-calibration",
    class = "lna_error_analysis_contract"
  )
})

test_that("public profile writes, validates, checks splits, and reads subsets", {
  fixture <- qpp_profile_fixture()
  profile <- lna_profile_analysis_safe(
    labels = fixture$labels, delta = 0.01,
    signal_units = "percent-signal-change", source_id = "atlas-v1",
    normalization_provenance = qpp_fixed_normalization(),
    map_grid_id = fixture$grid_id,
    time_chunk = 4L, parcel_chunk = 2L
  )
  file <- local_tempfile(fileext = ".h5")
  write_lna(
    fixture$x, file, transforms = profile$transforms,
    transform_params = profile$transform_params, mask = fixture$mask,
    header = list(
      dims = c(8L, 1L, 1L), spacing = c(2, 2, 2), origin = c(0, 0, 0),
      transform = fixture$affine, spatial_units = "mm",
      voxel_ordering = "first-axis-fastest"
    ),
    run_id = "run-01", checksum = "none"
  )

  receipt <- validate_lna(
    file, checksum = FALSE, qpp = "all", receipt = TRUE
  )
  expect_s3_class(receipt, "lna_validation_receipt")
  expect_true(receipt$qpp[[1L]]$q_roundtrip_verified)
  contract <- analysis_contract(file)
  expect_identical(contract$allowed_split_unit, "pointwise")
  folds <- data.frame(observation_id = seq_len(12), fold = rep(c("A", "B"), 6))
  expect_true(assert_split_compatible(file, folds)$compatible)

  roi <- array(FALSE, dim = dim(fixture$mask))
  roi[c(1, 8)] <- TRUE
  time_idx <- c(11L, 2L, 5L)
  out <- read_lna(
    file, roi_mask = roi, time_idx = time_idx, validate = FALSE
  )
  expect_identical(dim(out$stash$input), c(3L, 2L))
  expect_true(out$meta$qpp_io_receipt$physical_time_partial)
  expect_true(out$meta$qpp_io_receipt$physical_parcel_partial)
})

test_that("multiresolution profile is an executable no-go decision", {
  profile <- lna_profile_multiresolution_run()
  expect_s3_class(profile, "lna_profile_unavailable")
  expect_false(profile$available)
  expect_identical(profile$lifecycle, "prototype-no-go")
  expect_identical(profile$allowed_split_unit, "run")
  expect_identical(
    profile$allowed_split_description,
    "each transformed run must remain in one fold"
  )
  expect_null(profile$transforms)
  expect_match(profile$lifecycle_reason, "95.317%")
})
