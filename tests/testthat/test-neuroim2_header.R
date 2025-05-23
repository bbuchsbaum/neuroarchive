library(testthat)
library(neuroarchive)

# Helper to remove globals after test
cleanup_fake <- function() {
  objects_to_rm <- c("FakeSpace", "dim.FakeSpace", "spacing", "spacing.FakeSpace",
                     "origin", "origin.FakeSpace", "trans.FakeSpace", "ndim",
                     "space", "space.DenseNeuroVec", "trans")
  
  for (obj in objects_to_rm) {
    if (exists(obj, envir = .GlobalEnv)) {
      rm(list = obj, envir = .GlobalEnv)
    }
  }
}


test_that("neuroim2_space_to_lna_header converts NeuroSpace", {
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
  ndim <- function(x) length(dim(x))

  assign("FakeSpace", FakeSpace, envir = .GlobalEnv)
  assign("dim.FakeSpace", dim.FakeSpace, envir = .GlobalEnv)
  assign("spacing", spacing, envir = .GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir = .GlobalEnv)
  assign("origin", origin, envir = .GlobalEnv)
  assign("origin.FakeSpace", origin.FakeSpace, envir = .GlobalEnv)
  assign("trans", trans, envir = .GlobalEnv)
  assign("trans.FakeSpace", trans.FakeSpace, envir = .GlobalEnv)
  assign("ndim", ndim, envir = .GlobalEnv)

  withr::defer(cleanup_fake(), envir = parent.frame())

  sp <- FakeSpace(c(10, 11, 12, 2), c(1, 1, 1), c(0, 0, 0), diag(4))
  hdr <- neuroim2_space_to_lna_header(sp)
  expect_equal(hdr$dims, c(10, 11, 12))
  expect_equal(hdr$spacing, c(1, 1, 1))
  expect_equal(hdr$origin, c(0, 0, 0))
  expect_equal(hdr$transform, diag(4))
})



test_that("write_lna auto-populates header from NeuroObj", {
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

  withr::defer(cleanup_fake(), envir = parent.frame())

  sp <- FakeSpace(c(2,2,2,1), c(1,1,1), c(0,0,0), diag(4))
  x <- array(1, dim = c(2,2,2,1))
  # DenseNeuroVec should inherit from NeuroObj in neuroim2
  neuro_obj <- structure(x, class = c("DenseNeuroVec", "NeuroObj"))
  attr(neuro_obj, "space") <- sp

  # Test that our fake objects work
  expect_true(methods::is(neuro_obj, "NeuroObj"))
  expect_equal(space(neuro_obj), sp)

  captured <- list()
  local_mocked_bindings(
    core_write = function(x, transforms, transform_params, mask = NULL,
                          header = NULL, plugins = NULL, run_id = NULL) {
      captured$header <<- header
      captured$x <<- x
      list(handle = DataHandle$new(), plan = Plan$new())
    },
    materialise_plan = function(...) list(),
    .env = asNamespace("neuroarchive")
  )

  write_lna(neuro_obj, file = tempfile(fileext = ".h5"), transforms = character(0))

  # Debug: check what we captured
  expect_true(methods::is(captured$x, "DenseNeuroVec"))
  expect_false(is.null(captured$header))

  expect_equal(captured$header$dims, c(2,2,2))
  expect_equal(captured$header$spacing, c(1,1,1))
  expect_equal(captured$header$origin, c(0,0,0))
  expect_equal(captured$header$transform, diag(4))
})


