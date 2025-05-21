library(testthat)
library(hdf5r)
library(withr)

# Test materialise_plan basic functionality

test_that("materialise_plan creates structure and updates plan", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- neuroarchive:::open_h5(tmp, mode = "w")
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
  neuroarchive:::close_h5_safely(h5)
})

test_that("materialise_plan writes header attributes", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- neuroarchive:::open_h5(tmp, mode = "w")
  plan <- Plan$new()

  materialise_plan(h5, plan, header = list(vox = 1L, note = "hi"))

  expect_true(h5$exists("header/global"))
  grp <- h5[["header/global"]]
  expect_identical(h5_attr_read(grp, "vox"), 1L)
  expect_identical(h5_attr_read(grp, "note"), "hi")
  neuroarchive:::close_h5_safely(h5)
})

test_that("materialise_plan respects progress handlers", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- neuroarchive:::open_h5(tmp, mode = "w")
  plan <- Plan$new()
  for (i in 1:3) {
    key <- paste0("p", i)
    path <- paste0("/scans/run-01/d", i)
    plan$add_payload(key, 1:5)
    plan$add_dataset_def(path, "data", "dummy", "run-01", 0L, "{}", key, "eager")
  }
  progressr::handlers(progressr::handler_void)
  expect_silent(progressr::with_progress(materialise_plan(h5, plan)))
  progressr::handlers(NULL)
  neuroarchive:::close_h5_safely(h5)
})
