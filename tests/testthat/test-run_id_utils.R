library(testthat)
library(hdf5r)

# Tests for discover_run_ids and resolve_run_ids

test_that("discover_run_ids returns sorted run identifiers", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- H5File$new(tmp, mode = "w")
  scans <- h5$create_group("scans")
  scans$create_group("run-03")
  scans$create_group("run-01")
  scans$create_group("run-02")
  scans$create_group("misc")

  runs <- neuroarchive:::discover_run_ids(h5)
  expect_equal(runs, c("run-01", "run-02", "run-03"))
  h5$close_all()
})

test_that("resolve_run_ids matches names and globs uniquely", {
  available <- c("run-01", "run-02", "run-03")
  patterns <- c("run-*", "run-02")
  res <- neuroarchive:::resolve_run_ids(patterns, available)
  expect_equal(res, c("run-01", "run-02", "run-03"))
})
