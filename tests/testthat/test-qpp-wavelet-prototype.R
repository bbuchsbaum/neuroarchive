source(test_inst_file(
  "prototypes", "temporal-lift53-reference.R"
), local = TRUE)
source(test_inst_file(
  "benchmarks", "qpp-rate-distortion.R"
), local = TRUE)
source(test_inst_file(
  "benchmarks", "qpp-wavelet-prototype.R"
), local = TRUE)

qpp_small_wavelet_court <- local({
  receipt <- NULL
  function() {
    if (is.null(receipt)) {
      receipt <<- qpp_run_wavelet_prototype_court(
        n_runs = 1L, n_time = 72L, n_parcel = 8L,
        levels = 3L, repetitions = 1L
      )
    }
    receipt
  }
})

test_that("wavelet prototype compares CDF 5/3 and Haar lossless references", {
  receipt <- qpp_small_wavelet_court()
  candidates <- receipt$calibration_candidates
  ids <- vapply(candidates, function(x) x$candidate$id, "")
  expect_true(all(c("cdf53-lossless", "haar-lossless") %in% ids))
  expect_true(receipt$lossless_exact_for_both_wavelets)
  exact <- candidates[ids %in% c("cdf53-lossless", "haar-lossless")]
  expect_true(all(vapply(
    exact, `[[`, logical(1), "full_coefficient_q_roundtrip_exact"
  )))
})

test_that("sparse band codec measures actual bitmap and bit-packed bytes", {
  payload <- .qpp_wavelet_band_bytes(c(0, -1, 2), 1, 1L, "detail")
  expect_identical(payload$bitmap_bytes, 1L)
  expect_identical(payload$bit_width, 3L)
  expect_identical(payload$packed_value_bytes, 1L)
  expect_identical(payload$total_bytes, 26L)
  expect_equal(payload$zero_fraction, 1 / 3)
})

test_that("wavelet candidate selection uses only calibration-pass bytes", {
  receipt <- qpp_small_wavelet_court()
  candidates <- receipt$calibration_candidates
  passing <- which(vapply(candidates, function(x) x$gate$pass, logical(1)))
  bytes <- vapply(candidates[passing], function(x) x$rate$total_bytes, numeric(1))
  expect_identical(receipt$selected_candidate_index,
                   passing[[which.min(bytes)]])
  expect_match(receipt$selection_rule, "calibration")
  expect_true(receipt$corpus_manifest$disjoint)
  expect_identical(receipt$corpus_manifest$delta_fit_scope,
                   "external-calibration")
})

test_that("wavelet admission is fail-closed against scientific and rate gates", {
  receipt <- qpp_small_wavelet_court()
  expected <- isTRUE(receipt$evaluation_gate$pass) &&
    receipt$wavelet_rate_fraction_of_qpp <=
      receipt$gate_definition$admission_rate_fraction_of_qpp
  expect_identical(receipt$admit_to_qpp12, expected)
  expect_false(receipt$admit_to_qpp12)
  expect_gt(receipt$wavelet_rate_fraction_of_qpp, 0.90)
  expect_match(receipt$admission_reason, "do not integrate")
  expect_match(receipt$qpp_comparator$byte_boundary, "payload only")
})

test_that("package boundary keeps stable decoder in neuroarchive core", {
  receipt <- qpp_small_wavelet_court()
  expect_identical(
    receipt$package_boundary$production_decoder_owner, "neuroarchive-core"
  )
  expect_identical(
    receipt$package_boundary$fmrilatent_role,
    "fitting, calibration, and experimental transforms"
  )
  expect_false(receipt$package_boundary$fmrilatent_source_edited)
  expect_identical(
    receipt$package_boundary$fmrilatent_companion_mote,
    "bd-01M0WJX0QVS0ZE07EF3VCN3AG4"
  )
  expect_false(
    receipt$selected_calibration$partition_contract$cross_unit_operations
  )
  expect_identical(
    receipt$selected_calibration$partition_contract$lossy_support,
    "complete-run"
  )
})

test_that("committed wavelet receipt records the admission decision", {
  path <- testthat::test_path(
    "..", "..", "inst", "benchmarks", "qpp-wavelet-prototype-v1.json"
  )
  skip_if_not(file.exists(path), "wavelet prototype receipt is not generated")
  receipt <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  expect_identical(receipt$profile_status, "prototype")
  expect_true(receipt$lossless_exact_for_both_wavelets)
  expect_false(receipt$admit_to_qpp12)
  expect_false(receipt$package_boundary$fmrilatent_source_edited)
})
