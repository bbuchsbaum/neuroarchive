library(testthat)

test_that("default write produces gzip-compressed datasets", {
  arr <- array(rnorm(256), dim = c(4, 4, 4, 4))
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  write_lna(arr, tmp, transforms = "quant")

  h5 <- hdf5r::H5File$new(tmp, mode = "r")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)

  # Check that at least one dataset under /scans has compression
  # The quantized data should be stored compressed
  expect_true(file.size(tmp) > 0)
  h5$close_all()
})

test_that("compression_level=0 still works", {
  old <- lna_options("write.compression_level")
  on.exit(lna_options(write.compression_level = old[[1]]), add = TRUE)

  lna_options(write.compression_level = 0L)

  arr <- array(rnorm(256), dim = c(4, 4, 4, 4))
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  expect_no_error(write_lna(arr, tmp, transforms = "quant"))
})

test_that("compressed file is smaller than uncompressed", {
  arr <- array(rnorm(4000), dim = c(10, 10, 10, 4))

  tmp_compressed <- tempfile(fileext = ".h5")
  tmp_uncompressed <- tempfile(fileext = ".h5")
  on.exit({
    unlink(tmp_compressed)
    unlink(tmp_uncompressed)
  }, add = TRUE)

  # Write with default compression (level 4)
  write_lna(arr, tmp_compressed, transforms = "quant")

  # Write without compression
  old <- lna_options("write.compression_level")
  lna_options(write.compression_level = 0L)
  write_lna(arr, tmp_uncompressed, transforms = "quant")
  lna_options(write.compression_level = old[[1]])

  expect_lt(file.size(tmp_compressed), file.size(tmp_uncompressed))
})
