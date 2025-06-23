library(testthat)
library(hdf5r)
library(withr)


test_that("invert_step.delta applies roi_mask subset", {
  arr <- matrix(seq_len(20), nrow = 4, ncol = 5)
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(arr, file = tmp, transforms = "delta",
            transform_params = list(delta = list(axis = 2)))
  roi <- c(TRUE, FALSE, TRUE, FALSE)
  h <- read_lna(tmp, roi_mask = roi)
  out <- h$stash$input
  expect_equal(dim(out), c(sum(roi), ncol(arr)))
  expect_equal(out, arr[roi, , drop = FALSE])
})


test_that("invert_step.delta applies time_idx subset", {
  arr <- matrix(seq_len(20), nrow = 4, ncol = 5)
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(arr, file = tmp, transforms = "delta",
            transform_params = list(delta = list(axis = 2)))
  idx <- c(1, 5)
  h <- read_lna(tmp, time_idx = idx)
  out <- h$stash$input
  expect_equal(dim(out), c(nrow(arr), length(idx)))
  expect_equal(out, arr[, idx, drop = FALSE])
})


test_that("invert_step.delta applies roi_mask and time_idx subset", {
  arr <- matrix(seq_len(20), nrow = 4, ncol = 5)
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(arr, file = tmp, transforms = "delta",
            transform_params = list(delta = list(axis = 2)))
  roi <- c(TRUE, FALSE, TRUE, FALSE)
  idx <- c(2, 5)
  h <- read_lna(tmp, roi_mask = roi, time_idx = idx)
  out <- h$stash$input
  expect_equal(dim(out), c(sum(roi), length(idx)))
  expect_equal(out, arr[roi, idx, drop = FALSE])
})

