library(testthat)

# Tests for convenience alias functions

test_that("compress_fmri is a direct passthrough wrapper", {
  expect_true(is.pairlist(formals(neuroarchive::compress_fmri)))
  expect_identical(names(formals(neuroarchive::compress_fmri)), "...")
  # Under covr, function bodies are instrumented so exact body match fails
  if (!(requireNamespace("covr", quietly = TRUE) && covr::in_covr())) {
    expect_identical(body(neuroarchive::compress_fmri), quote(write_lna(...)))
  }
})

test_that("open_lna is an alias of read_lna", {
  # Explicitly reference functions from the namespace to ensure they are found
  # This assumes devtools::load_all() has correctly loaded the package.
  expect_identical(neuroarchive::open_lna, neuroarchive::read_lna)
})
