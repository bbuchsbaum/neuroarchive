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
  expect_error(core_read(tmp, output_dtype = "float16"),
               class = "lna_error_float16_unsupported")
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
