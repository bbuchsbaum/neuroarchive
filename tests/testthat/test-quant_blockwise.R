library(testthat)
library(hdf5r)

# Block-wise forward_step.quant should match in-memory quantization

test_that("block-wise voxel scope matches full array", {
  set.seed(1)
  arr <- array(runif(2*3*4*5), dim = c(20, 10, 6, 5))
  baseline <- neuroarchive:::.quantize_voxel(arr, bits = 8, method = "range", center = TRUE)

  tmp <- local_tempfile(fileext = ".h5")
  h5 <- neuroarchive:::open_h5(tmp, mode = "w")
  plan <- Plan$new()
  handle <- DataHandle$new(initial_stash = list(input = arr), plan = plan,
                           h5 = h5, run_ids = "run-01", current_run_id = "run-01")
  desc <- list(type = "quant", params = list(scale_scope = "voxel", bits = 8),
               inputs = c("input"))
  neuroarchive:::forward_step.quant("quant", desc, handle)

  root <- h5[["/"]]
  q_disk <- neuroarchive:::h5_read(root, "scans/run-01/quantized")
  sc_disk <- neuroarchive:::h5_read(root, "scans/run-01/quant_scale")
  off_disk <- neuroarchive:::h5_read(root, "scans/run-01/quant_offset")
  expect_equal(q_disk, baseline$q)
  expect_equal(sc_disk, baseline$scale)
  expect_equal(off_disk, baseline$offset)
  expect_equal(handle$meta$quant_stats$n_clipped_total, 0L)
  neuroarchive:::close_h5_safely(h5)
})
