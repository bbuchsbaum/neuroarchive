library(testthat)
library(hdf5r)

# forward_step.quant pre-creates datasets when handle has h5

test_that("forward_step.quant precreates datasets for voxel scope", {
  arr <- array(runif(24), dim = c(4,3,2,1))
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- open_h5(tmp, mode = "w")
  plan <- Plan$new()
  handle <- DataHandle$new(initial_stash = list(input = arr), plan = plan,
                           h5 = h5, run_ids = "run-01", current_run_id = "run-01")
  desc <- list(type = "quant", params = list(scale_scope = "voxel", bits = 8),
               inputs = c("input"))
  forward_step.quant("quant", desc, handle)
  root <- h5[["/"]]
  expect_true(root$exists("scans/run-01/quantized"))
  expect_true(root$exists("scans/run-01/quant_scale"))
  expect_true(root$exists("scans/run-01/quant_offset"))
  dset <- root[["scans/run-01/quantized"]]
  res <- auto_block_size(dim(arr)[1:3], element_size_bytes = 1)
  expected_chunk <- c(res$slab_dims, dim(arr)[4])
  dcpl <- dset$get_create_plist()
  expect_equal(dcpl$get_chunk(4), expected_chunk)
  dset$close(); dcpl$close(); close_h5_safely(h5)
})
