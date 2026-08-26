library(testthat)
library(hdf5r)

test_that("is_nbit_dataset detects N-bit on materialised writes", {
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)
  f <- H5File$new(tmp, mode = "w")
  root <- f[["/"]]
  codes <- as.integer(sample(0:7, 64, replace = TRUE))
  d <- neuroarchive:::h5_write_qbits_nbit(root, "/qbits", codes, bits = 3)
  expect_true(neuroarchive::is_nbit_dataset(d))
  back <- d$read()
  expect_equal(as.integer(back), codes)
  d$close(); f$close_all()
})

test_that("is_nbit_dataset detects N-bit on precreated", {
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)
  f <- H5File$new(tmp, mode = "w")
  root <- f[["/"]]
  ds <- neuroarchive:::h5_create_empty_dataset_nbit(root, "/blk/q7", dims = c(10), bits = 7)
  # write some values
  vals <- as.integer(0:9) %% 8L
  ds$write(args = NULL, value = vals)
  expect_true(neuroarchive::is_nbit_dataset(ds))
  back <- ds$read()
  expect_equal(as.integer(back), vals)
  ds$close(); f$close_all()
})
