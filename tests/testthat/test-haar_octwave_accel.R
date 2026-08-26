library(neuroarchive)

test_that("Morton indices ordering is deterministic", {
  mask <- array(FALSE, dim = c(8, 5, 4))
  mask[cbind(1:8, rep(3, 8), rep(2, 8))] <- TRUE

  set.seed(123)
  seed_before <- .Random.seed
  idx1 <- neuroarchive:::get_morton_ordered_indices(mask, z_order_seed = 42)
  seed_after <- .Random.seed
  expect_identical(seed_before, seed_after)

  idx2 <- neuroarchive:::get_morton_ordered_indices(mask, z_order_seed = 1234)
  expect_identical(idx1, idx2)
})

test_that("Detail widths shrink with level", {
  opts <- options(lna.hwt.use_rcpp = FALSE)
  on.exit(options(opts), add = TRUE)

  mask <- array(TRUE, dim = c(8, 8, 8))
  X <- matrix(rnorm(4 * sum(mask)), nrow = 4)
  levels <- 3L

  coeff <- neuroarchive:::perform_haar_lift_analysis(X, mask, levels)
  scalings <- neuroarchive:::precompute_haar_scalings(mask, levels)
  counts <- lapply(scalings, function(sc) as.integer(round(sc$sqrt_nvalid^2)))

  expect_equal(ncol(coeff$detail[[1]]), sum(counts[[1]]))
  expect_equal(ncol(coeff$detail[[2]]), length(counts[[1]]))
  expect_equal(ncol(coeff$detail[[3]]), length(counts[[2]]))
})

test_that("Forward and inverse lifting agree in R and Rcpp", {
  mask <- array(FALSE, dim = c(6, 6, 4))
  mask[sample.int(length(mask), 90)] <- TRUE
  X <- matrix(rnorm(10 * sum(mask)), nrow = 10)
  levels <- 3L

  opts <- options(lna.hwt.use_rcpp = FALSE)
  coeff_r <- neuroarchive:::perform_haar_lift_analysis(X, mask, levels)
  reco_r <- neuroarchive:::perform_haar_lift_synthesis(coeff_r, mask, levels)
  options(opts)

  coeff_cpp <- neuroarchive:::perform_haar_lift_analysis(X, mask, levels)
  reco_cpp <- neuroarchive:::perform_haar_lift_synthesis(coeff_cpp, mask, levels)

  expect_equal(reco_r, X, tolerance = 1e-8)
  expect_equal(reco_cpp, X, tolerance = 1e-8)
})

test_that("Non-singleton top level completes without error", {
  mask <- array(TRUE, dim = c(16, 8, 8))
  res <- neuroarchive:::perform_haar_lift_analysis(matrix(0, 2, sum(mask)), mask, levels = 2)
  expect_true(is.list(res))
  expect_true("root" %in% names(res))
})
