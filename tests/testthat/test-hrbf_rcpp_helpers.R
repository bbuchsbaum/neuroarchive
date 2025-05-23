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


test_that("label_components Rcpp vs R fallback", {
  mask <- array(FALSE, dim = c(3,3,3))
  mask[1:2,1:2,1:2] <- TRUE
  mask[3,3,3] <- TRUE

  withr::defer({
    rm(FakeSpace, space.FakeLogicalNeuroVol, spacing.FakeSpace,
       as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())
  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.FakeLogicalNeuroVol", space.FakeLogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("as.array.FakeLogicalNeuroVol", as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)

  opts <- options(lna.hrbf.use_rcpp_helpers = TRUE)
  res_cpp <- neuroarchive:::label_components(mask)
  options(lna.hrbf.use_rcpp_helpers = FALSE)
  res_R <- neuroarchive:::label_components(mask)
  options(opts)

  expect_equal(res_cpp$count, res_R$count)
  expect_identical(res_cpp$labels, res_R$labels)
  expect_identical(res_cpp$labels > 0, mask)
})


test_that("poisson_disk_sample_neuroim2 Rcpp vs R fallback and radius", {
  mask <- array(TRUE, dim = c(5,5,5))
  vol <- structure(list(arr = mask), class = "LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(5,5,5), c(1,1,1))

  withr::defer({
    rm(FakeSpace, space.FakeLogicalNeuroVol, spacing.FakeSpace,
       as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())
  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.FakeLogicalNeuroVol", space.FakeLogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("as.array.FakeLogicalNeuroVol", as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)

  opts <- options(lna.hrbf.use_rcpp_helpers = TRUE)
  c_cpp <- neuroarchive:::poisson_disk_sample_neuroim2(vol, radius_mm = 2, seed = 42)
  c_cpp2 <- neuroarchive:::poisson_disk_sample_neuroim2(vol, radius_mm = 2, seed = 42)
  options(lna.hrbf.use_rcpp_helpers = FALSE)
  c_R <- neuroarchive:::poisson_disk_sample_neuroim2(vol, radius_mm = 2, seed = 42)
  options(opts)

  expect_identical(c_cpp, c_cpp2)
  expect_identical(c_cpp, c_R)
  if (nrow(c_cpp) > 1) {
    d2 <- as.matrix(dist(c_cpp))^2
    expect_true(all(d2[upper.tri(d2)] >= 4))
  }
})


test_that("rcpp helpers performance smoke test", {
  skip_on_cran()
  if (!exists("label_components_6N_rcpp") ||
      !exists("poisson_disk_sample_component_rcpp")) {
    skip("Rcpp helpers not available")
  }
  mask <- array(TRUE, dim = rep(128L,3))
  vol <- structure(list(arr = mask), class = "LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(rep(128L,3), c(1,1,1))

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.FakeLogicalNeuroVol", space.FakeLogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("as.array.FakeLogicalNeuroVol", as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  withr::defer({
    rm(FakeSpace, space.FakeLogicalNeuroVol, spacing.FakeSpace,
       as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  opts <- options(lna.hrbf.use_rcpp_helpers = TRUE)
  t1 <- system.time(res <- neuroarchive:::label_components(mask))["elapsed"]
  expect_lt(t1, 0.1)
  t2 <- system.time(neuroarchive:::poisson_disk_sample_neuroim2(vol, radius_mm = 2, seed = 1))["elapsed"]
  expect_lt(t2, 0.1)
  options(opts)
})


