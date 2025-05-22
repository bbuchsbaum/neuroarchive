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
  expect_false(p$allow_clip)
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
  expect_equal(res$n_clipped_total, 0L)
  expect_equal(res$clip_pct, 0)

  res2 <- neuroarchive:::.quantize_global(x, bits = 8, method = "range", center = FALSE)
  expect_equal(res2$scale, 1)
  expect_true(all(res2$q == 0))
})

test_that(".quantize_global counts clipping", {
  set.seed(1)
  x <- c(rep(0, 98), 5, -5)
  res <- neuroarchive:::.quantize_global(x, bits = 8, method = "sd", center = TRUE)
  expect_equal(res$n_clipped_total, 2L)
  expect_equal(res$clip_pct, 2)
  expect_true(all(res$q >= 0 & res$q <= 255))
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

test_that("forward_step.quant stores clipping stats in handle meta", {
  arr <- c(rep(0, 98), 5, -5)
  tmp <- local_tempfile(fileext = ".h5")
  res <- write_lna(arr, file = tmp, transforms = "quant",
                   transform_params = list(quant = list(method = "sd")))
  expect_equal(res$handle$meta$quant_stats$n_clipped_total, 2L)
  expect_equal(res$handle$meta$quant_stats$clip_pct, 2)
})

test_that("invert_step.quant warns when quant_bits attribute missing", {
  arr <- array(runif(6), dim = c(2,3))
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(arr, file = tmp, transforms = "quant")
  expect_warning(read_lna(tmp), regexp = "quant_bits")
})

test_that("quantized dataset uses uint8 or uint16 storage", {
  arr <- array(runif(6), dim = c(2,3))

  tmp8 <- local_tempfile(fileext = ".h5")
  write_lna(arr, file = tmp8, transforms = "quant",
           transform_params = list(quant = list(bits = 7)))
  h5 <- H5File$new(tmp8, mode = "r")
  dset <- h5[["scans/run-01/quantized"]]
  dt <- hdf5r:::datatype_to_char(dset$get_type())
  expect_equal(dt, "H5T_STD_U8LE")
  expect_equal(h5_attr_read(dset, "quant_bits"), 7L)
  st <- h5[["scans/run-01/quant_scale"]]
  ot <- h5[["scans/run-01/quant_offset"]]
  expect_equal(hdf5r:::datatype_to_char(st$get_type()), "H5T_IEEE_F32LE")
  expect_equal(hdf5r:::datatype_to_char(ot$get_type()), "H5T_IEEE_F32LE")
  dset$close(); st$close(); ot$close(); h5$close_all()

  tmp16 <- local_tempfile(fileext = ".h5")
  write_lna(arr, file = tmp16, transforms = "quant",
           transform_params = list(quant = list(bits = 12)))
  h52 <- H5File$new(tmp16, mode = "r")
  dset2 <- h52[["scans/run-01/quantized"]]
  dt2 <- hdf5r:::datatype_to_char(dset2$get_type())
  expect_equal(dt2, "H5T_STD_U16LE")
  expect_equal(h5_attr_read(dset2, "quant_bits"), 12L)
  st2 <- h52[["scans/run-01/quant_scale"]]
  ot2 <- h52[["scans/run-01/quant_offset"]]
  expect_equal(hdf5r:::datatype_to_char(st2$get_type()), "H5T_IEEE_F32LE")
  expect_equal(hdf5r:::datatype_to_char(ot2$get_type()), "H5T_IEEE_F32LE")
  dset2$close(); st2$close(); ot2$close(); h52$close_all()
})

test_that("quant_bits attribute is validated against descriptor", {
  arr <- array(runif(6), dim = c(2,3))
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(arr, file = tmp, transforms = "quant",
           transform_params = list(quant = list(bits = 9)))
  expect_silent(read_lna(tmp))

  h5 <- H5File$new(tmp, mode = "a")
  dset <- h5[["scans/run-01/quantized"]]
  h5_attr_write(dset, "quant_bits", 8L)
  dset$close(); h5$close_all()

  expect_error(read_lna(tmp), class = "lna_error_validation")
})

test_that("missing quant_bits attribute triggers warning", {
  arr <- array(runif(6), dim = c(2,3))
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(arr, file = tmp, transforms = "quant")
  h5 <- H5File$new(tmp, mode = "a")
  dset <- h5[["scans/run-01/quantized"]]
  if (h5_attr_exists(dset, "quant_bits")) h5_attr_delete(dset, "quant_bits")
  dset$close(); h5$close_all()
  expect_warning(read_lna(tmp), regexp = "quant_bits")
})

test_that("forward_step.quant warns or errors based on clipping thresholds", {
  opts_env <- get(".lna_opts", envir = neuroarchive:::lna_options_env)
  withr::defer(rm(list = c("quant.clip_warn_pct", "quant.clip_abort_pct"),
                  envir = opts_env))
  lna_options(quant.clip_warn_pct = 0.5, quant.clip_abort_pct = 5)

  arr_warn <- c(rep(0, 98), 100, -100)
  tmp_warn <- local_tempfile(fileext = ".h5")
  expect_warning(
    write_lna(arr_warn, file = tmp_warn, transforms = "quant",
              transform_params = list(quant = list(method = "sd"))),
    regexp = "Clipping"
  )

  arr_err <- c(rep(0, 94), rep(100, 6))
  tmp_err <- local_tempfile(fileext = ".h5")
  expect_error(
    write_lna(arr_err, file = tmp_err, transforms = "quant",
              transform_params = list(quant = list(method = "sd"))),
    class = "lna_error_validation"
  )

  tmp_allow <- local_tempfile(fileext = ".h5")
  expect_warning(
    write_lna(arr_err, file = tmp_allow, transforms = "quant",
              transform_params = list(quant = list(method = "sd",
                                                   allow_clip = TRUE))),
    regexp = "Clipping"
  )
})

test_that("forward_step.quant hard clips output range", {
  arr <- c(rep(0, 98), 10, -10)
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(arr, file = tmp, transforms = "quant",
            transform_params = list(quant = list(method = "sd")))
  h5 <- H5File$new(tmp, mode = "r")
  dset <- h5[["scans/run-01/quantized"]]
  qvals <- dset$read()
  expect_true(all(qvals >= 0 & qvals <= 255))
  dset$close(); h5$close_all()
})
