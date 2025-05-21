library(testthat)
library(hdf5r)
library(withr)

# Test materialise_plan basic functionality

test_that("materialise_plan creates structure and updates plan", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- H5File$new(tmp, mode = "w")
  plan <- Plan$new()
  plan$add_descriptor("00_dummy.json", list(type = "dummy"))
  plan$add_payload("payload", matrix(1:4, nrow = 2))
  plan$add_dataset_def("/scans/run-01/data", "data", "dummy", "run-01", 0L, "{}", "payload", "eager")

  materialise_plan(h5, plan)

  expect_true(h5$exists("transforms"))
  expect_true(h5$exists("basis"))
  expect_true(h5$exists("scans"))

  root <- h5[["/"]]
  expect_identical(h5_attr_read(root, "lna_spec"), "LNA R v2.0")
  expect_identical(h5_attr_read(root, "creator"), "lna R package v0.0.1")
  expect_identical(h5_attr_read(root, "required_transforms"), character(0))

  desc <- read_json_descriptor(h5[["transforms"]], "00_dummy.json")
  expect_identical(desc, list(type = "dummy"))

  expect_equal(plan$datasets$write_mode_effective, "eager")
  expect_true(is.null(plan$payloads$payload))
  expect_true(h5$exists("scans/run-01/data"))
  expect_equal(h5[["scans/run-01/data"]]$read(), matrix(1:4, nrow = 2))
  expect_true(h5$is_valid())
  h5$close_all()
})
