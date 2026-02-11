library(testthat)

test_that("successful write creates target file", {
  arr <- array(rnorm(64), dim = c(4, 4, 4, 1))
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  write_lna(arr, tmp, transforms = "quant")
  expect_true(file.exists(tmp))

  # No temp files should remain
  tmp_files <- list.files(dirname(tmp), pattern = "\\.lna\\.tmp$", full.names = TRUE)
  expect_length(tmp_files, 0)
})

test_that("failed write preserves existing file", {
  arr <- array(rnorm(64), dim = c(4, 4, 4, 1))
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  # Write a valid file first
  write_lna(arr, tmp, transforms = "quant")
  original_size <- file.size(tmp)
  expect_true(original_size > 0)

  # Now attempt a write that will fail (unknown transform)
  tryCatch(
    write_lna(arr, tmp, transforms = "nonexistent_transform_xyz"),
    error = function(e) NULL
  )

  # Original file should be preserved (same size)
  expect_true(file.exists(tmp))
  expect_equal(file.size(tmp), original_size)
})

test_that("temp files cleaned up after failure", {
  arr <- array(rnorm(64), dim = c(4, 4, 4, 1))
  tmp <- tempfile(fileext = ".h5")
  tmp_dir <- dirname(tmp)

  # Count existing .lna.tmp files before
  before <- length(list.files(tmp_dir, pattern = "\\.lna\\.tmp$"))

  tryCatch(
    write_lna(arr, tmp, transforms = "nonexistent_transform_xyz"),
    error = function(e) NULL
  )

  # No new .lna.tmp files should remain
  after <- length(list.files(tmp_dir, pattern = "\\.lna\\.tmp$"))
  expect_equal(after, before)
})

test_that("in-memory path still works", {
  arr <- array(rnorm(64), dim = c(4, 4, 4, 1))

  # file=NULL should still work (in-memory mode)
  expect_no_error(write_lna(arr, NULL, transforms = "quant"))
})
