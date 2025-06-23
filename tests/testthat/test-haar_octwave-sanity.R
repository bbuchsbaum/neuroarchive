library(testthat)

# Additional sanity tests for Haar lifting

test_that("energy is preserved by forward/inverse", {
  mask <- make_ball_mask(radius = 1L)
  X <- array(rnorm(sum(mask) * 3), dim = c(3, sum(mask)))
  opts <- options(lna.hwt.use_rcpp = FALSE)
  coeff <- neuroarchive:::perform_haar_lift_analysis(t(X), mask, levels = 1)
  reco <- neuroarchive:::perform_haar_lift_synthesis(coeff, mask, levels = 1)
  options(opts)
  expect_equal(sum(X^2), sum(reco^2), tolerance = 1e-6)
})
