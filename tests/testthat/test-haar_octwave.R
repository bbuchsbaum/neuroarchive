library(testthat)
library(neuroarchive)

# Basic tests for Haar octwave helpers

context("Haar lifting utilities")

set.seed(1)

test_that("precompute_haar_scalings handles simple masks", {
  mask_full <- array(TRUE, dim = c(2,2,2))
  s_full <- neuroarchive:::precompute_haar_scalings(mask_full, 1)
  expect_equal(length(s_full), 1L)
  expect_equal(s_full[[1]]$sqrt_nvalid, sqrt(8))

  mask_edge <- array(FALSE, dim = c(2,2,2))
  mask_edge[1,1,1] <- TRUE
  s_edge <- neuroarchive:::precompute_haar_scalings(mask_edge, 1)
  expect_equal(s_edge[[1]]$sqrt_nvalid, sqrt(1))

  mask_disc <- array(FALSE, dim = c(2,2,2))
  mask_disc[1,1,1] <- TRUE
  mask_disc[2,2,2] <- TRUE
  s_disc <- neuroarchive:::precompute_haar_scalings(mask_disc, 1)
  expect_equal(sum(round(s_disc[[1]]$sqrt_nvalid^2)), 2)
})


test_that("forward/inverse lifting reconstructs signal", {
  mask <- array(TRUE, dim = c(2,2,2))
  X <- matrix(rnorm(16), nrow = 2)

  opts <- options(lna.hwt.use_rcpp = FALSE)
  coeff <- neuroarchive:::perform_haar_lift_analysis(X, mask, levels = 1)
  reco  <- neuroarchive:::perform_haar_lift_synthesis(coeff, mask, levels = 1)
  options(opts)

  expect_equal(reco, X, tolerance = 1e-6)
})
