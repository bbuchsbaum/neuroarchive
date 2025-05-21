library(testthat)
library(neuroarchive)


test_that("sanitize_run_id accepts valid id", {
  expect_equal(neuroarchive:::sanitize_run_id("run-01"), "run-01")
})

test_that("sanitize_run_id rejects invalid patterns", {
  expect_error(neuroarchive:::sanitize_run_id("run1"), class = "lna_error_validation")
  expect_error(neuroarchive:::sanitize_run_id("run-1"), class = "lna_error_validation")
  expect_error(neuroarchive:::sanitize_run_id("run/01"), class = "lna_error_validation")
  expect_error(neuroarchive:::sanitize_run_id("run\\01"), class = "lna_error_validation")
})
