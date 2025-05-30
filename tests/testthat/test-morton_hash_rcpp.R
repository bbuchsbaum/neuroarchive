library(testthat)

# Test Rcpp morton hash helper

test_that("morton_indices_to_hash_rcpp matches R version", {
  idx <- as.integer(c(1,5,9,2,3))
  r_val <- neuroarchive:::morton_indices_to_hash(idx)
  opts <- options(lna.hwt.use_rcpp = TRUE)
  expect_true(exists("morton_indices_to_hash_rcpp"))
  rc_val <- neuroarchive:::morton_indices_to_hash(idx)
  options(opts)
  expect_identical(rc_val, r_val)
})
