library(testthat)
library(hdf5r)
library(withr)

# Tests for open_h5 and close_h5_safely helpers

test_that("open_h5 creates and closes files", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- neuroarchive:::open_h5(tmp, mode = "w")
  expect_s3_class(h5, "H5File")
  expect_true(h5$is_valid())
  neuroarchive:::close_h5_safely(h5)
  expect_false(h5$is_valid())
})

test_that("open_h5 errors for missing file", {
  missing <- file.path(tempdir(), "does_not_exist.h5")
  expect_error(neuroarchive:::open_h5(missing, mode = "r"), "Failed to open HDF5 file")
})

test_that("close_h5_safely tolerates invalid objects", {
  expect_silent(neuroarchive:::close_h5_safely(NULL))
  expect_silent(neuroarchive:::close_h5_safely("not a handle"))
})
