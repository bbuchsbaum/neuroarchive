library(testthat)
library(neuroarchive)

# Helper utilities for fake neuroim2 objects
FakeSpace <- function(dim, spacing_v) {
  structure(list(dim = dim, spacing = spacing_v, trans = diag(4), origin = c(0,0,0)),
            class = "FakeSpace")
}
space.FakeLogicalNeuroVol <- function(x, ...) attr(x, "space")
spacing.FakeSpace <- function(x, ...) x$spacing
as.array.FakeLogicalNeuroVol <- function(x, ...) x$arr


test_that("poisson_disk_sample_neuroim2 deterministic", {
  mask <- array(TRUE, dim = c(5,5,5))
  vol <- structure(list(arr = mask), class = "LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(5,5,5), c(1,1,1))

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.FakeLogicalNeuroVol", space.FakeLogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("as.array.FakeLogicalNeuroVol", as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  withr::defer({
    rm(FakeSpace, space.FakeLogicalNeuroVol, spacing.FakeSpace,
       as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  c1 <- neuroarchive:::poisson_disk_sample_neuroim2(vol, radius_mm = 2, seed = 42)
  c2 <- neuroarchive:::poisson_disk_sample_neuroim2(vol, radius_mm = 2, seed = 42)
  expect_identical(c1, c2)
})


test_that("poisson_disk_sample_neuroim2 guard rail on tiny ROI", {
  mask <- array(TRUE, dim = c(2,2,2))
  vol <- structure(list(arr = mask), class = "LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(2,2,2), c(1,1,1))

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.FakeLogicalNeuroVol", space.FakeLogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("as.array.FakeLogicalNeuroVol", as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  withr::defer({
    rm(FakeSpace, space.FakeLogicalNeuroVol, spacing.FakeSpace,
       as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  pts <- neuroarchive:::poisson_disk_sample_neuroim2(vol, radius_mm = 5, seed = 1)
  expect_equal(nrow(pts), 1)
})


test_that("poisson_disk_sample_neuroim2 handles disconnected components", {
  mask <- array(FALSE, dim = c(4,4,4))
  mask[1:2,1:2,1:2] <- TRUE
  mask[3:4,3:4,3:4] <- TRUE
  vol <- structure(list(arr = mask), class = "LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(4,4,4), c(1,1,1))

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.FakeLogicalNeuroVol", space.FakeLogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("as.array.FakeLogicalNeuroVol", as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  withr::defer({
    rm(FakeSpace, space.FakeLogicalNeuroVol, spacing.FakeSpace,
       as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  pts <- neuroarchive:::poisson_disk_sample_neuroim2(vol, radius_mm = 1, seed = 99)
  expect_true(any(pts[,1] <= 2))
  expect_true(any(pts[,1] > 2))
})


test_that("generate_hrbf_atom gaussian normalisation", {
  coords <- matrix(rbind(c(0,0,0), c(1,0,0), c(0,1,0)), ncol = 3, byrow = TRUE)
  idx <- 1:3
  res <- neuroarchive:::generate_hrbf_atom(coords, idx, c(0,0,0), sigma_mm = 1,
                                           kernel_type = "gaussian")
  expect_equal(res$indices, idx)
  expect_equal(length(res$values), 3)
  expect_equal(sum(res$values^2), 1, tolerance = 1e-6)
})

test_that("generate_hrbf_atom wendland_c4 normalisation", {
  coords <- matrix(rbind(c(0,0,0), c(1,0,0), c(0,1,0)), ncol = 3, byrow = TRUE)
  idx <- 1:3
  res <- neuroarchive:::generate_hrbf_atom(coords, idx, c(0,0,0), sigma_mm = 2,
                                           kernel_type = "wendland_c4")
  expect_equal(res$indices, idx)
  expect_equal(sum(res$values^2), 1, tolerance = 1e-6)
  expect_true(all(res$values >= 0))
})
