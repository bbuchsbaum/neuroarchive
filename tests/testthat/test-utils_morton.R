library(testthat)
library(digest)

# Tests for Morton order helpers

test_that("get_morton_ordered_indices works on 2x2x2 full mask", {
  mask <- array(TRUE, dim = c(2, 2, 2))
  idx <- neuroarchive:::get_morton_ordered_indices(mask, z_order_seed = 1L)
  expect_equal(idx, 1:8)
})

test_that("get_morton_ordered_indices handles empty mask", {
  mask <- array(FALSE, dim = c(2, 2, 2))
  idx <- neuroarchive:::get_morton_ordered_indices(mask)
  expect_length(idx, 0)
})

test_that("morton_indices_to_hash computes sha1 hash", {
  vec <- c(1L, 2L, 3L)
  h <- neuroarchive:::morton_indices_to_hash(vec)
  # Check format is correct (sha1: prefix + 40 hex chars)
  expect_match(h, "^sha1:[a-f0-9]{40}$")
  # Check it's deterministic
  h2 <- neuroarchive:::morton_indices_to_hash(vec)
  expect_identical(h, h2)
  # Different input should give different hash
  h3 <- neuroarchive:::morton_indices_to_hash(c(1L, 2L, 4L))
  expect_false(identical(h, h3))
})
