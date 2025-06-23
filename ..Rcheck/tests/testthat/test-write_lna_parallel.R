library(testthat)
library(withr)

# Test recommended atomic rename pattern for parallel writes

test_that("write_lna temp file can be atomically renamed", {
  dir <- local_tempdir()
  final <- file.path(dir, "dest.lna.h5")
  tmp <- tempfile(tmpdir = dir, fileext = ".h5")

  res <- write_lna(x = array(1, dim = c(1, 1, 1)), file = tmp, transforms = character(0))
  expect_true(file.exists(tmp))
  expect_true(file.rename(tmp, final))
  expect_false(file.exists(tmp))
  expect_true(file.exists(final))

  h5 <- neuroarchive:::open_h5(final, mode = "r")
  root <- h5[["/"]]
  expect_identical(h5_attr_read(root, "lna_spec"), "LNA R v2.0")
  neuroarchive:::close_h5_safely(h5)
})
