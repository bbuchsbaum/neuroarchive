library(testthat)

test_that("core_write rejects unknown transform types early", {
  arr <- array(rnorm(64), dim = c(4, 4, 4, 1))
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  expect_error(
    write_lna(arr, tmp, transforms = "nonexistent_transform_xyz"),
    "Unknown transform"
  )
})

test_that("core_write accepts known transform types", {
  arr <- array(rnorm(64), dim = c(4, 4, 4, 1))
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  expect_no_error(
    write_lna(arr, tmp, transforms = "quant")
  )
})
