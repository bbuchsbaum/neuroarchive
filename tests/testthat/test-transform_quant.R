library(testthat)
#library(neuroarchive)
library(hdf5r)
library(withr)


test_that("default_params for quant loads schema", {
  neuroarchive:::default_param_cache_clear()
  p <- neuroarchive:::default_params("quant")
  expect_equal(p$bits, 8)
  expect_equal(p$method, "range")
  expect_true(p$center)
})


test_that("quant transform forward and inverse roundtrip", {
  arr <- array(runif(12), dim = c(3,4))
  tmp <- local_tempfile(fileext = ".h5")

  res <- write_lna(arr, file = tmp, transforms = "quant")
  expect_equal(nrow(res$plan$datasets), 3)

  h <- read_lna(tmp)
  out <- h$stash$input
  expect_equal(dim(out), dim(arr))
  expect_lt(mean(abs(out - arr)), 1)
})

test_that("quant transform supports sd method and voxel scope", {
  arr <- array(runif(40), dim = c(2,2,2,5))
  tmp <- local_tempfile(fileext = ".h5")

  res <- write_lna(arr, file = tmp, transforms = "quant",
                   transform_params = list(quant = list(method = "sd",
                                                         scale_scope = "voxel")))
  expect_equal(nrow(res$plan$datasets), 3)

  h <- read_lna(tmp)
  out <- h$stash$input
  expect_equal(dim(out), dim(arr))
  expect_lt(mean(abs(out - arr)), 1)
})

test_that("invert_step.quant applies roi_mask and time_idx", {
  arr <- array(seq_len(40), dim = c(2,2,2,5))
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(arr, file = tmp, transforms = "quant")
  roi <- array(c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE), dim = c(2,2,2))
  h <- read_lna(tmp, roi_mask = roi, time_idx = c(2,5))
  out <- h$stash$input
  expect_equal(dim(out), c(sum(roi), 2))
})


test_that("forward_step.quant validates parameters", {
  plan <- Plan$new()
  h <- DataHandle$new(initial_stash = list(input = array(1:8, dim = c(2,4))),
                      plan = plan)

  desc <- list(type = "quant", params = list(bits = 0L), inputs = c("input"))
  expect_error(
    neuroarchive:::forward_step.quant("quant", desc, h),
    class = "lna_error_validation",
    regexp = "bits"
  )

  desc$params$bits <- 17L
  expect_error(
    neuroarchive:::forward_step.quant("quant", desc, h),
    class = "lna_error_validation"
  )

  desc$params <- list(method = "bad")
  expect_error(
    neuroarchive:::forward_step.quant("quant", desc, h),
    class = "lna_error_validation",
    regexp = "method"
  )

  desc$params <- list(center = c(TRUE, FALSE))
  expect_error(
    neuroarchive:::forward_step.quant("quant", desc, h),
    class = "lna_error_validation",
    regexp = "center"
  )

  desc$params <- list(scale_scope = "nonsense")
  expect_error(
    neuroarchive:::forward_step.quant("quant", desc, h),
    class = "lna_error_validation",
    regexp = "scale_scope"
  )
})

test_that(".quantize_global handles constant arrays", {
  x <- rep(5, 10)
  res <- neuroarchive:::.quantize_global(x, bits = 8, method = "range", center = TRUE)
  expect_equal(res$scale, 1)
  expect_true(all(res$q == 0))

  res2 <- neuroarchive:::.quantize_global(x, bits = 8, method = "range", center = FALSE)
  expect_equal(res2$scale, 1)
  expect_true(all(res2$q == 0))
})

test_that(".quantize_voxel handles constant arrays", {
  arr <- array(5, dim = c(2,2,2,3))
  res <- neuroarchive:::.quantize_voxel(arr, bits = 8, method = "range", center = TRUE)
  expect_true(all(as.numeric(res$scale) == 1))
  expect_true(all(res$q == 0))

  res2 <- neuroarchive:::.quantize_voxel(arr, bits = 8, method = "range", center = FALSE)
  expect_true(all(as.numeric(res2$scale) == 1))
  expect_true(all(res2$q == 0))
})
          
test_that("quant transform errors on non-finite input", {
  arr <- c(1, NA, 3)
  tmp <- local_tempfile(fileext = ".h5")
  expect_error(
    write_lna(arr, file = tmp, transforms = "quant"),
    class = "lna_error_validation",
    regexp = "non-finite"
  )

  arr_nan <- c(1, NaN, 3)
  tmp2 <- local_tempfile(fileext = ".h5")
  expect_error(
    write_lna(arr_nan, file = tmp2, transforms = "quant"),
    class = "lna_error_validation",
    regexp = "non-finite"
  )

})
