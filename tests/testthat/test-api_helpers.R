library(testthat)
library(hdf5r)
library(withr)
library(neuroarchive)

# Tests for open_output_h5 and close_output_h5

test_that("open_output_h5 creates in-memory file and cleans up", {
  info <- neuroarchive:::open_output_h5(NULL)
  expect_true(info$in_memory)
  expect_true(inherits(info$h5, "H5File"))
  temp_path <- info$file
  expect_true(file.exists(temp_path))
  expect_true(info$h5$is_valid)
  neuroarchive:::close_output_h5(info)
  expect_false(file.exists(temp_path))
})


test_that("open_output_h5 persists disk file after close", {
  tmp <- local_tempfile(fileext = ".h5")
  info <- neuroarchive:::open_output_h5(tmp)
  expect_false(info$in_memory)
  expect_identical(info$file, tmp)
  expect_true(file.exists(tmp))
  expect_true(info$h5$is_valid)
  neuroarchive:::close_output_h5(info)
  expect_true(file.exists(tmp))
  h5 <- neuroarchive:::open_h5(tmp, mode = "r")
  expect_true(h5$is_valid)
  neuroarchive:::close_h5_safely(h5)
})

# Test derive_header_from_input directly

test_that("derive_header_from_input uses space() on NeuroObj", {
  FakeSpace <- function(dim, spacing_v, origin_v, trans_m) {
    structure(list(dim = dim, spacing = spacing_v,
                   origin = origin_v, trans = trans_m),
              class = "FakeSpace")
  }
  dim.FakeSpace <- function(x) x$dim
  spacing <- function(x, ...) UseMethod("spacing")
  spacing.FakeSpace <- function(x, ...) x$spacing
  origin <- function(x, ...) UseMethod("origin")
  origin.FakeSpace <- function(x, ...) x$origin
  trans <- function(x, ...) UseMethod("trans")
  trans.FakeSpace <- function(x, ...) x$trans
  space <- function(x, ...) UseMethod("space")
  space.DenseNeuroVec <- function(x, ...) attr(x, "space")

  assign("FakeSpace", FakeSpace, envir = .GlobalEnv)
  assign("dim.FakeSpace", dim.FakeSpace, envir = .GlobalEnv)
  assign("spacing", spacing, envir = .GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir = .GlobalEnv)
  assign("origin", origin, envir = .GlobalEnv)
  assign("origin.FakeSpace", origin.FakeSpace, envir = .GlobalEnv)
  assign("trans", trans, envir = .GlobalEnv)
  assign("trans.FakeSpace", trans.FakeSpace, envir = .GlobalEnv)
  assign("space", space, envir = .GlobalEnv)
  assign("space.DenseNeuroVec", space.DenseNeuroVec, envir = .GlobalEnv)

  withr::defer({
    rm(list = c("FakeSpace", "dim.FakeSpace", "spacing", "spacing.FakeSpace",
                 "origin", "origin.FakeSpace", "trans", "trans.FakeSpace",
                 "space", "space.DenseNeuroVec"), envir = .GlobalEnv)
  }, envir = parent.frame())

  sp <- FakeSpace(c(3,3,3,1), c(1,1,1), c(0,0,0), diag(4))
  x <- array(1, dim = c(3,3,3,1))
  neuro_obj <- structure(x, class = c("DenseNeuroVec", "NeuroObj"))
  attr(neuro_obj, "space") <- sp

  hdr <- neuroarchive:::derive_header_from_input(list(neuro_obj))
  expect_equal(hdr$dims, c(3,3,3))
  expect_equal(hdr$spacing, c(1,1,1))
  expect_equal(hdr$origin, c(0,0,0))
  expect_equal(hdr$transform, diag(4))
})

