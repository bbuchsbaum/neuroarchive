library(testthat)
library(hdf5r)
library(neuroarchive)

# Tests for inverse, sparsification and full roundtrip of spat.haar_octwave

# Basic inverse reconstruction without quant

test_that("invert_step.spat.haar_octwave reconstructs input", {
  set.seed(1)
  arr <- array(rnorm(16), dim = c(2,2,2,2))
  mask <- array(TRUE, dim = c(2,2,2))
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(arr, file = tmp, transforms = "spat.haar_octwave",
            mask = mask,
            transform_params = list(spat.haar_octwave = list(levels = 1)))
  h <- read_lna(tmp)
  out <- h$stash$input
  expect_equal(dim(out), c(dim(arr)[4], sum(mask)))
  expect_equal(out, as_dense_mat(arr), tolerance = 1e-6)
})

# Sparsification threshold introduces zeros but reconstruction still succeeds

test_that("detail sparsification zeros coefficients", {
  arr <- array(runif(16), dim = c(2,2,2,2))
  mask <- array(TRUE, dim = c(2,2,2))
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(arr, file = tmp, transforms = "spat.haar_octwave",
            mask = mask,
            transform_params = list(spat.haar_octwave = list(levels = 1,
                                                               detail_threshold_type = "absolute",
                                                               detail_threshold_value = 100)))
  h5 <- H5File$new(tmp, mode = "r")
  dset <- h5[["/scans/run-01/wavelet/level_0/detail_coefficients"]]
  coeff <- dset$read()
  expect_true(all(coeff == 0))
  dset$close(); h5$close_all()
  h <- read_lna(tmp)
  out <- h$stash$input
  expect_equal(dim(out), c(dim(arr)[4], sum(mask)))
})

# Full roundtrip with quant following spat.haar_octwave

test_that("spat.haar_octwave roundtrip with quant", {
  set.seed(2)
  arr <- array(runif(16), dim = c(2,2,2,2))
  mask <- array(TRUE, dim = c(2,2,2))
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(arr, file = tmp, mask = mask,
            transforms = c("spat.haar_octwave", "quant"),
            transform_params = list(quant = list(bits = 8L)))
  h <- read_lna(tmp)
  out <- h$stash$input
  expect_equal(dim(out), c(dim(arr)[4], sum(mask)))
  expect_lt(mean(abs(out - as_dense_mat(arr))), 1)
})

# ROI streaming returns subset of voxels

test_that("read_lna applies roi_mask after spat.haar_octwave", {
  arr <- array(runif(16), dim = c(2,2,2,2))
  mask <- array(TRUE, dim = c(2,2,2))
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(arr, file = tmp, mask = mask, transforms = "spat.haar_octwave")
  roi <- mask
  roi[1,1,1] <- FALSE
  h <- read_lna(tmp, roi_mask = roi)
  out <- h$stash$input
  expect_equal(dim(out), c(dim(arr)[4], sum(roi)))
})
