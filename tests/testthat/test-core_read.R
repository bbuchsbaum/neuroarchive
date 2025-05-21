library(testthat)
library(hdf5r)
library(withr)

# Ensure core_read and helpers are loaded
# source("../R/core_read.R")

# Helper to create empty transforms group
create_empty_lna <- function(path) {
  h5 <- neuroarchive:::open_h5(path, mode = "w")
  h5$create_group("transforms")
  neuroarchive:::close_h5_safely(h5)
}

# Helper to create lna with one dummy descriptor
create_dummy_lna <- function(path) {
  h5 <- neuroarchive:::open_h5(path, mode = "w")
  tf <- h5$create_group("transforms")
  write_json_descriptor(tf, "00_dummy.json", list(type = "dummy"))
  neuroarchive:::close_h5_safely(h5)
}

test_that("core_read handles empty /transforms group", {
  tmp <- local_tempfile(fileext = ".h5")
  create_empty_lna(tmp)

  handle <- core_read(tmp)
  expect_true(inherits(handle, "DataHandle"))
  expect_false(handle$h5$is_valid())
})

test_that("core_read closes file if invert_step errors", {
  tmp <- local_tempfile(fileext = ".h5")
  create_dummy_lna(tmp)

  captured_h5 <- NULL
  with_mocked_bindings(
    invert_step.default = function(type, desc, handle) {
      captured_h5 <<- handle$h5
      stop("mock error")
    },
    expect_error(core_read(tmp), "mock error")
  )
  expect_true(inherits(captured_h5, "H5File"))
  expect_false(captured_h5$is_valid())
})

test_that("core_read lazy=TRUE keeps file open", {
  tmp <- local_tempfile(fileext = ".h5")
  create_empty_lna(tmp)

  handle <- core_read(tmp, lazy = TRUE)
  expect_true(handle$h5$is_valid())
  neuroarchive:::close_h5_safely(handle$h5)
})

test_that("core_read output_dtype stored and float16 check", {
  tmp <- local_tempfile(fileext = ".h5")
  create_empty_lna(tmp)

  h <- core_read(tmp, output_dtype = "float64")
  expect_equal(h$meta$output_dtype, "float64")
  err <- expect_error(
    core_read(tmp, output_dtype = "float16"),
    class = "lna_error_float16_unsupported"
  )
  expect_true(grepl("core_read", err$location))
})

test_that("core_read allows float16 when support present", {
  tmp <- local_tempfile(fileext = ".h5")
  create_empty_lna(tmp)

  with_mocked_bindings(
    has_float16_support = function() TRUE,
    {
      h <- core_read(tmp, output_dtype = "float16")
      expect_equal(h$meta$output_dtype, "float16")
      expect_true(h$h5$is_valid())
      neuroarchive:::close_h5_safely(h$h5)
    }
  )
})

test_that("core_read works with progress handlers", {
  tmp <- local_tempfile(fileext = ".h5")
  create_dummy_lna(tmp)
  progressr::handlers(progressr::handler_void)
  expect_silent(progressr::with_progress(core_read(tmp)))
  progressr::handlers(NULL)
})

test_that("core_read validate=TRUE calls validate_lna", {
  tmp <- local_tempfile(fileext = ".h5")
  create_empty_lna(tmp)
  called <- FALSE
  with_mocked_bindings(
    validate_lna = function(file) { called <<- TRUE },
    { core_read(tmp, validate = TRUE) }
  )
  expect_true(called)
})

test_that("core_read allow_plugins='none' errors on unknown transform", {
  tmp <- local_tempfile(fileext = ".h5")
  create_dummy_lna(tmp)
  expect_error(core_read(tmp, allow_plugins = "none"),
               class = "lna_error_no_method")
})

test_that("core_read allow_plugins='prompt' falls back when non-interactive", {
  tmp <- local_tempfile(fileext = ".h5")
  create_dummy_lna(tmp)
  with_mocked_bindings(
    rlang::is_interactive = function() FALSE,
    { expect_warning(core_read(tmp, allow_plugins = "prompt"),
                    "Missing invert_step") }
  )
})

test_that("core_read allow_plugins='prompt' interactive respects user choice", {
  tmp <- local_tempfile(fileext = ".h5")
  create_dummy_lna(tmp)
  with_mocked_bindings(
    rlang::is_interactive = function() TRUE,
    readline = function(prompt = "") "n",
    { expect_error(core_read(tmp, allow_plugins = "prompt"),
                  class = "lna_error_no_method") }
  )

  with_mocked_bindings(
    rlang::is_interactive = function() TRUE,
    readline = function(prompt = "") "y",
    { expect_warning(core_read(tmp, allow_plugins = "prompt"),
                    "Missing invert_step") }
  )
})


test_that("core_read stores subset parameters", {
  tmp <- local_tempfile(fileext = ".h5")
  create_empty_lna(tmp)
  roi <- array(TRUE, dim = c(1,1,1))
  h <- core_read(tmp, roi_mask = roi, time_idx = 1:2)
  expect_identical(h$subset$roi_mask, roi)
  expect_identical(h$subset$time_idx, 1:2)
  expect_false(h$h5$is_valid())
})

test_that("core_read run_id globbing returns handles", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- neuroarchive:::open_h5(tmp, mode = "w")
  plan <- Plan$new()
  plan$add_descriptor("00_dummy.json", list(type = "dummy"))
  plan$add_payload("p1", 1L)
  plan$add_dataset_def("/scans/run-01/data", "data", "dummy", "run-01", 0L, "{}", "p1", "eager")
  plan$add_payload("p2", 2L)
  plan$add_dataset_def("/scans/run-02/data", "data", "dummy", "run-02", 0L, "{}", "p2", "eager")
  materialise_plan(h5, plan)
  neuroarchive:::close_h5_safely(h5)

  with_mocked_bindings(
    invert_step.dummy = function(type, desc, handle) {
      path <- paste0("/scans/", handle$current_run_id, "/data")
      root <- handle$h5[["/"]]
      val <- h5_read(root, path)
      handle$update_stash(keys = character(), new_values = list(input = val))
    }, {
      res <- core_read(tmp, run_id = "run-0*")
    }
  )
  expect_true(is.list(res))
  expect_equal(names(res), c("run-01", "run-02"))
  expect_equal(res[["run-01"]]$stash$input, 1)
  expect_equal(res[["run-02"]]$stash$input, 2)
})

test_that("core_read run_id globbing lazy returns first", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- neuroarchive:::open_h5(tmp, mode = "w")
  plan <- Plan$new()
  plan$add_descriptor("00_dummy.json", list(type = "dummy"))
  plan$add_payload("p1", 1L)
  plan$add_dataset_def("/scans/run-01/data", "data", "dummy", "run-01", 0L, "{}", "p1", "eager")
  plan$add_payload("p2", 2L)
  plan$add_dataset_def("/scans/run-02/data", "data", "dummy", "run-02", 0L, "{}", "p2", "eager")
  materialise_plan(h5, plan)
  neuroarchive:::close_h5_safely(h5)

  with_mocked_bindings(
    invert_step.dummy = function(type, desc, handle) {
      path <- paste0("/scans/", handle$current_run_id, "/data")
      root <- handle$h5[["/"]]
      val <- h5_read(root, path)
      handle$update_stash(keys = character(), new_values = list(input = val))
    }, {
      expect_warning(h <- core_read(tmp, run_id = "run-*", lazy = TRUE), "first match")
    }
  )
  expect_s3_class(h, "DataHandle")
  expect_equal(h$current_run_id, "run-01")
  neuroarchive:::close_h5_safely(h$h5)

})
