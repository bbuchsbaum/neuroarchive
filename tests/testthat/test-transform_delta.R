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
  expect_equal(drop(out), arr)
})

test_that("forward_step.delta uses custom desc$outputs for stash", {
  plan <- Plan$new()
  handle <- DataHandle$new(initial_stash = list(input = matrix(1:4, nrow = 2)),
                           plan = plan, run_ids = "run-01",
                           current_run_id = "run-01")
  desc <- list(type = "delta", params = list(), inputs = c("input"),
               outputs = c("my_delta"))

  h <- neuroarchive:::forward_step.delta("delta", desc, handle)

  expect_true(h$has_key("my_delta"))
  expect_false(h$has_key("delta_stream"))
})


test_that("delta transform with rle coding works", {
  arr <- matrix(rep(1:5, each = 2), nrow = 5, ncol = 2)
  tmp <- local_tempfile(fileext = ".h5")

  res <- write_lna(arr, file = tmp, transforms = "delta",
                   transform_params = list(delta = list(coding_method = "rle")))
  expect_true(file.exists(tmp))

  h <- read_lna(tmp)
  out <- h$stash$input
  expect_equal(drop(out), arr)

  p <- neuroarchive:::default_params("delta")

  expect_equal(p$axis, -1L)
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

test_that("rle coding compresses delta stream for 1D input", {
  arr <- rep(0, 10)
  tmp <- local_tempfile(fileext = ".h5")

  write_lna(arr, file = tmp, transforms = "delta",
           transform_params = list(delta = list(coding_method = "rle")))
  expect_true(file.exists(tmp))

  h5 <- H5File$new(tmp, mode = "r")
  ds_path <- "/scans/run-01/deltas/00_delta/delta_stream"
  delta_stream <- h5[[ds_path]][]
  h5$close_all()

  expect_true(nrow(delta_stream) < length(arr) - 1)

  h <- read_lna(tmp)
  out <- h$stash$input
  expect_equal(drop(out), arr)
})

test_that("rle coding compresses delta stream for matrix input", {
  arr <- matrix(0, nrow = 4, ncol = 5)
  tmp <- local_tempfile(fileext = ".h5")

  write_lna(arr, file = tmp, transforms = "delta",
           transform_params = list(delta = list(coding_method = "rle")))
  expect_true(file.exists(tmp))

  h5 <- H5File$new(tmp, mode = "r")
  ds_path <- "/scans/run-01/deltas/00_delta/delta_stream"
  delta_stream <- h5[[ds_path]][]
  h5$close_all()

  expect_true(nrow(delta_stream) < prod(dim(arr)) - nrow(arr))

  h <- read_lna(tmp)
  out <- h$stash$input
  expect_equal(drop(out), arr)
})
