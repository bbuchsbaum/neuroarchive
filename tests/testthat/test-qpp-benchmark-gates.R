qpp_benchmark_script <- test_inst_file(
  "benchmarks", "qpp-rate-distortion.R"
)
source(qpp_benchmark_script, local = TRUE)

qpp_small_court <- local({
  receipt <- NULL
  function() {
    if (is.null(receipt)) {
      receipt <<- qpp_run_rate_distortion_court(
        n_runs = 1L, n_time = 72L, n_parcel = 8L,
        alpha = c(0.00025, 0.001),
        time_chunk = c(32L, 64L), parcel_chunk = c(4L, 8L),
        repetitions = 1L
      )
    }
    receipt
  }
})

test_that("QPP benchmark corpus separates calibration from evaluation", {
  calibration <- qpp_benchmark_corpus("calibration", 2L, 72L, 8L)
  evaluation <- qpp_benchmark_corpus("evaluation", 2L, 72L, 8L)
  calibration_ids <- vapply(calibration, `[[`, "", "id")
  evaluation_ids <- vapply(evaluation, `[[`, "", "id")
  expect_length(intersect(calibration_ids, evaluation_ids), 0L)
  expect_true(all(vapply(calibration, `[[`, "", "split") == "calibration"))
  expect_true(all(vapply(evaluation, `[[`, "", "split") == "evaluation"))
  expect_false(identical(
    .qpp_bench_hash_corpus(calibration), .qpp_bench_hash_corpus(evaluation)
  ))
})

test_that("QPP court reports explicit rate runtime fidelity and uncertainty", {
  receipt <- qpp_small_court()
  expect_identical(receipt$receipt_version, "1.0.0")
  expect_identical(receipt$profile_status, "experimental")
  expect_false(receipt$stable_profile_claim)
  expect_true(receipt$corpus_manifest$disjoint)
  expect_true(receipt$calibration_gate$pass)
  expect_true(receipt$evaluation_gate$pass)
  expect_true(all(unlist(receipt$edge_case_q_roundtrip)))
  expect_identical(receipt$outer_codec_decision$selected, "raw")

  rate_names <- c(
    "payload_bytes", "index_bytes", "step_bytes", "metadata_bytes",
    "total_transform_bytes", "effective_payload_bits_per_parcel_TR",
    "effective_total_bits_per_parcel_TR",
    "measured_parcel_float32_compression_ratio", "tile_count"
  )
  expect_true(all(rate_names %in% names(receipt$selected_calibration$rate)))
  expect_named(
    receipt$selected_calibration$runtime,
    c("encode", "sequential_decode", "random_tile_decode",
      "exploratory_gzip_outer", "allocation_profile", "peak_memory")
  )
  fidelity_names <- c(
    "maximum_absolute_error", "RMSE", "spectrum_relative_RMSE",
    "autocorrelation_max_absolute_error", "glm_beta_relative_RMSE",
    "glm_contrast_max_absolute_error", "glm_t_max_absolute_error",
    "connectivity_edge_RMSE", "connectivity_edge_max_absolute_error",
    "single_trial_beta_RMSE", "mvpa_prediction_agreement",
    "mvpa_decision_value_RMSE"
  )
  expect_named(receipt$selected_calibration$fidelity, fidelity_names)
  expect_true(nzchar(receipt$gate_definition$estimand))
  expect_true(nzchar(receipt$gate_definition$aggregation))
  expect_true(nzchar(receipt$execution$uncertainty))
  expect_gte(receipt$execution$repetitions, 1L)
  expect_true(nzchar(receipt$execution$platform))
})

test_that("QPP candidate is selected only from calibration-pass frontier", {
  receipt <- qpp_small_court()
  candidates <- receipt$calibration_candidates
  passing <- which(vapply(candidates, `[[`, logical(1), "pass"))
  ratios <- vapply(candidates[passing], `[[`, numeric(1),
                   "parcel_float32_compression_ratio")
  expected <- passing[[which.max(ratios)]]
  expect_identical(receipt$selected_candidate_index, expected)
  expect_match(receipt$selection_rule, "calibration")
  expect_match(receipt$corpus_manifest$evaluation$selection_role, "untouched")
  expect_true(receipt$corpus_manifest$external_data_promotion_required)
})

test_that("QPP baselines name actual byte boundaries and loss contracts", {
  receipt <- qpp_small_court()
  baselines <- receipt$baselines
  expect_identical(
    vapply(baselines, `[[`, "", "name"),
    c("uncompressed parcel float32", "parcel float32 plus gzip",
      "current quant uint8", "N-bit 12", "delta plus RLE", "QPP selected")
  )
  expect_true(all(vapply(baselines, function(x) x$bytes > 0, logical(1))))
  expect_true(all(vapply(baselines, function(x) nzchar(x$boundary), logical(1))))
  expect_true(all(vapply(baselines, function(x) nzchar(x$loss), logical(1))))
  expect_match(receipt$accounting$derived_identity, "V/K")
  expect_false(any(grepl("voxel.*ratio", names(receipt$selected_calibration$rate))))
})

test_that("QPP archive-level measurement validates every installed tile", {
  receipt <- qpp_small_court()
  archive <- receipt$archive_measurement
  expect_gt(archive$file_bytes, 0)
  expect_identical(archive$installed_decoder, "qpp@1.1.0")
  expect_identical(archive$tiles_checked, archive$tile_count)
  expect_match(archive$qpp_validation_receipt_checksum, "^sha256:[0-9a-f]{64}$")
})

test_that("committed QPP gate receipt stays explicitly non-promotional", {
  path <- testthat::test_path(
    "..", "..", "inst", "benchmarks", "qpp-gates-v1.json"
  )
  skip_if_not(file.exists(path), "benchmark receipt has not been generated yet")
  receipt <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  expect_identical(receipt$profile_status, "experimental")
  expect_false(receipt$stable_profile_claim)
  expect_true(receipt$corpus_manifest$disjoint)
  expect_true(receipt$corpus_manifest$external_data_promotion_required)
  expect_true(receipt$calibration_gate$pass)
  expect_false(receipt$evaluation_gate$pass)
  failed <- names(Filter(
    function(value) !isTRUE(value), receipt$evaluation_gate$checks
  ))
  expect_identical(failed, "connectivity_max")
})
