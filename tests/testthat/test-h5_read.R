library(testthat)
library(hdf5r)
library(withr)
library(neuroarchive)
# Tests for h5_read and h5_read_subset

test_that("h5_read returns dataset contents", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- H5File$new(tmp, mode = "w")
  root <- h5[["/"]]
  mat <- matrix(1:9, nrow = 3)
  root$create_dataset("mat", mat)

  read_back <- h5_read(root, "mat")
  expect_equal(read_back, mat)

  h5$close_all()
})

test_that("h5_read errors when dataset is missing", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- H5File$new(tmp, mode = "w")
  root <- h5[["/"]]

  expect_error(h5_read(root, "missing"), "Dataset 'missing' not found")
  h5$close_all()
})

test_that("h5_read_subset returns correct subset", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- H5File$new(tmp, mode = "w")
  root <- h5[["/"]]
  mat <- matrix(1:16, nrow = 4)
  root$create_dataset("mat", mat)

  sub <- h5_read_subset(root, "mat", list(1:2, 2:3))
  expect_equal(sub, mat[1:2, 2:3])
  h5$close_all()
})
