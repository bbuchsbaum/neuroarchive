library(testthat)

# Multi-level forward/inverse Haar lifting should reconstruct the input

test_that("multi-level Haar lifting roundtrip reconstructs X", {
  mask <- make_ball_mask(radius = 1L)
  arr  <- make_synthetic_fmri(time = 4L, mask = mask)
  X    <- t(neuroarchive:::convert_to_masked_vox_time_matrix(arr, mask))
  opts <- options(lna.hwt.use_rcpp = FALSE)
  coeff <- neuroarchive:::perform_haar_lift_analysis(X, mask, levels = 2)
  reco  <- neuroarchive:::perform_haar_lift_synthesis(coeff, mask, levels = 2)
  options(opts)
  expect_equal(reco, X, tolerance = 1e-6)
})

