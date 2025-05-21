library(testthat)
#library(neuroarchive)

library(hdf5r)
library(withr)


test_that("default_params for delta loads schema", {
  neuroarchive:::default_param_cache_clear()
  p <- neuroarchive:::default_params("delta")
  expect_equal(p$order, 1)

  expect_true(is.numeric(p$axis))
  expect_equal(p$reference_value_storage, "first_value_verbatim")
})


test_that("delta transform forward and inverse roundtrip", {
  arr <- matrix(1:10, nrow = 5, ncol = 2)
  tmp <- local_tempfile(fileext = ".h5")

  res <- write_lna(arr, file = tmp, transforms = "delta")
  expect_true(file.exists(tmp))

  h <- read_lna(tmp)
  out <- h$stash$input
  expect_equal(out, arr)
})


test_that("delta transform with rle coding works", {
  arr <- matrix(rep(1:5, each = 2), nrow = 5, ncol = 2)
  tmp <- local_tempfile(fileext = ".h5")

  res <- write_lna(arr, file = tmp, transforms = "delta",
                   transform_params = list(delta = list(coding_method = "rle")))
  expect_true(file.exists(tmp))

  h <- read_lna(tmp)
  out <- h$stash$input
  expect_equal(out, arr)

  p <- neuroarchive:::default_params("delta")

  expect_equal(p$axis, 4)
  expect_equal(p$reference_value_storage, "first_value_verbatim")
  expect_equal(p$coding_method, "none")

})

test_that("delta transform rejects unsupported coding_method", {
  arr <- matrix(1:4, nrow = 2)
  tmp <- local_tempfile(fileext = ".h5")

  expect_error(
    write_lna(arr, file = tmp, transforms = "delta",
              transform_params = list(delta = list(coding_method = "bogus"))),
    class = "lna_error_validation",
    regexp = "coding_method"
  )
})
