library(testthat)

test_that("pipeline subset profile distinguishes basis semantics", {
  q_caps <- neuroarchive:::.pipeline_supports_subset("quant")
  expect_true(q_caps$roi_mask)
  expect_equal(q_caps$roi_strategy, "direct")
  expect_false(q_caps$coefficient_pipeline)

  bq_caps <- neuroarchive:::.pipeline_supports_subset(c("basis", "quant"))
  expect_true(bq_caps$roi_mask)
  expect_equal(bq_caps$roi_strategy, "dense-pass-through")
  expect_false(bq_caps$coefficient_pipeline)
  expect_true(bq_caps$dense_passthrough_pipeline)

  beq_caps <- neuroarchive:::.pipeline_supports_subset(c("basis", "embed", "quant"))
  expect_true(beq_caps$roi_mask)
  expect_equal(beq_caps$roi_strategy, "coefficient")
  expect_true(beq_caps$coefficient_pipeline)
  expect_false(beq_caps$dense_passthrough_pipeline)

  sb_caps <- neuroarchive:::.pipeline_supports_subset(c("embed.shared_basis", "quant"))
  expect_true(sb_caps$roi_mask)
  expect_equal(sb_caps$roi_strategy, "shared-coefficient")
  expect_true(sb_caps$coefficient_pipeline)
})

test_that("pipeline subset profile respects unsupported transforms", {
  caps <- neuroarchive:::.pipeline_supports_subset(c("aggregate_runs", "quant"))
  expect_false(caps$roi_mask)
  expect_false(caps$time_idx)
})

test_that("transform contract classifies lifecycle and subset support", {
  expect_equal(neuroarchive:::.transform_lifecycle_status("quant"), "stable")
  expect_equal(neuroarchive:::.transform_lifecycle_status("spat.hrbf"), "experimental")
  expect_equal(neuroarchive:::.transform_lifecycle_status("aggregate_runs"), "developer")
  expect_equal(neuroarchive:::.transform_lifecycle_status("my_custom"), "custom")

  expect_equal(
    neuroarchive:::.transform_lifecycle_status("quant", list(lifecycle = "experimental")),
    "experimental"
  )

  q_caps <- neuroarchive:::.transform_supports_subset("quant")
  expect_true(q_caps$roi_mask)
  expect_true(q_caps$time_idx)

  aggregate_caps <- neuroarchive:::.transform_supports_subset("aggregate_runs")
  expect_false(aggregate_caps$roi_mask)
  expect_false(aggregate_caps$time_idx)
})

test_that("QPP capabilities distinguish semantics from physical tile I/O", {
  qpp <- neuroarchive:::.transform_supports_subset("qpp")
  expect_true(qpp$semantic$time_idx)
  expect_false(qpp$semantic$roi_mask)
  expect_true(qpp$semantic$parcel_idx)
  expect_true(qpp$physical$time_idx)
  expect_false(qpp$physical$roi_mask)
  expect_true(qpp$physical$parcel_idx)
  expect_identical(qpp$temporal_granularity, "qpp-time-tile")
  expect_identical(qpp$parcel_granularity, "qpp-parcel-tile")

  pipeline <- neuroarchive:::.pipeline_supports_subset(c("parcel.mean", "qpp"))
  expect_true(pipeline$roi_mask)
  expect_true(pipeline$time_idx)
  expect_true(pipeline$parcel_idx)
  expect_true(pipeline$physical_roi_mask)
  expect_true(pipeline$physical_time_idx)
  expect_true(pipeline$physical_parcel_idx)
  expect_identical(pipeline$roi_strategy, "parcel-gather")
  expect_match(pipeline$temporal_granularity, "qpp-time-tile")
  expect_match(pipeline$parcel_granularity, "qpp-parcel-tile")

  suppressWarnings(register_transform("qpp", validate = FALSE))
  info <- get_transform_info("qpp")
  expect_true(info$physical_time_idx)
  expect_false(info$physical_roi_mask)
  expect_true(info$physical_parcel_idx)
  expect_identical(info$temporal_granularity, "qpp-time-tile")
  expect_identical(info$roi_strategy, "parcel-index")

  quant <- neuroarchive:::.transform_supports_subset("quant")
  expect_true(quant$semantic$time_idx)
  expect_false(quant$physical$time_idx)
})
