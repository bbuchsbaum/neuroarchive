library(testthat)

test_that("h5_write_dataset creates dataset at path", {
  skip_if_not_installed("hdf5r")
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)

  neuroarchive:::h5_write_dataset(h5[["/"]], "data/values", 1:10)
  expect_true(h5$exists("data/values"))
  expect_equal(as.integer(h5[["data/values"]]$read()), 1:10)
  h5$close_all()
})

test_that("h5_write_dataset rejects non-string path", {
  skip_if_not_installed("hdf5r")
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)

  expect_error(neuroarchive:::h5_write_dataset(h5[["/"]], 42, 1:5), "single character string")
  h5$close_all()
})

test_that("h5_write_dataset handles leading slash", {
  skip_if_not_installed("hdf5r")
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)

  neuroarchive:::h5_write_dataset(h5[["/"]], "/group/dataset", c(1.5, 2.5, 3.5))
  expect_true(h5$exists("group/dataset"))
  h5$close_all()
})
