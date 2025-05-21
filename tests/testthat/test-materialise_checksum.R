library(testthat)
library(hdf5r)
library(withr)
library(digest)

# Test checksum writing

test_that("materialise_plan writes sha256 checksum attribute", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- neuroarchive:::open_h5(tmp, mode = "w")
  plan <- Plan$new()
  plan$add_descriptor("00_dummy.json", list(type = "dummy"))
  plan$add_payload("payload", matrix(1:4, nrow = 2))
  plan$add_dataset_def("/scans/run-01/data", "data", "dummy", "run-01", 0L, "{}", "payload", "eager")

  materialise_plan(h5, plan, checksum = "sha256")

  expect_false(h5$is_valid())
  h5r <- neuroarchive:::open_h5(tmp, mode = "r")
  root <- h5r[["/"]]
  expect_true(h5_attr_exists(root, "lna_checksum"))
  expected <- digest(file = tmp, algo = "sha256")
  expect_identical(h5_attr_read(root, "lna_checksum"), expected)
  neuroarchive:::close_h5_safely(h5r)
})
