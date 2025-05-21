library(testthat)
library(hdf5r)
library(withr)

# Ensure core_read and helpers are loaded
# source("../R/core_read.R")

# Helper to create empty transforms group
create_empty_lna <- function(path) {
  h5 <- H5File$new(path, mode = "w")
  h5$create_group("transforms")
  h5$close_all()
}

# Helper to create lna with one dummy descriptor
create_dummy_lna <- function(path) {
  h5 <- H5File$new(path, mode = "w")
  tf <- h5$create_group("transforms")
  write_json_descriptor(tf, "00_dummy.json", list(type = "dummy"))
  h5$close_all()
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
