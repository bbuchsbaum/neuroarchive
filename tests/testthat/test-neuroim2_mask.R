library(testthat)
library(neuroarchive)

test_that("LogicalNeuroVol mask warns on space mismatch", {
  FakeSpace <- function(dim, trans) structure(list(dim = dim, trans = trans), class = "FakeSpace")
  trans.FakeSpace <- function(x) x$trans
  dim.FakeSpace <- function(x) x$dim
  space <- function(x, ...) UseMethod("space")
  space.FakeLogicalNeuroVol <- function(x, ...) attr(x, "space")
  as.array.FakeLogicalNeuroVol <- function(x, ...) x$arr
  space.DenseNeuroVec <- function(x, ...) attr(x, "space")

  mask <- structure(list(arr = array(TRUE, dim = c(2,2,2))), class = c("FakeLogicalNeuroVol", "LogicalNeuroVol"))
  attr(mask, "space") <- FakeSpace(c(2,2,2), matrix(1:16, 4,4))

  # Create a simple object with the necessary attributes, not a list
  input_obj <- 1 
  attr(input_obj, "space") <- FakeSpace(c(2,2,2), diag(4))
  class(input_obj) <- "DenseNeuroVec"

  withr::defer({
    rm(FakeSpace, trans.FakeSpace, dim.FakeSpace, space, space.FakeLogicalNeuroVol,
       as.array.FakeLogicalNeuroVol, space.DenseNeuroVec, envir = .GlobalEnv)
  }, envir = parent.frame())

  assign("FakeSpace", FakeSpace, envir = .GlobalEnv)
  assign("trans.FakeSpace", trans.FakeSpace, envir = .GlobalEnv)
  assign("dim.FakeSpace", dim.FakeSpace, envir = .GlobalEnv)
  assign("space", space, envir = .GlobalEnv)
  assign("space.FakeLogicalNeuroVol", space.FakeLogicalNeuroVol, envir = .GlobalEnv)
  assign("as.array.FakeLogicalNeuroVol", as.array.FakeLogicalNeuroVol, envir = .GlobalEnv)
  assign("space.DenseNeuroVec", space.DenseNeuroVec, envir = .GlobalEnv)

  expect_warning(neuroarchive:::validate_mask(mask, input_obj),
                 "Mask orientation/space differs")
})

