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
  expected <- paste0("sha1:", digest::digest(vec, algo = "sha1", serialize = TRUE))
  expect_identical(h, expected)
  expect_match(h, "^sha1:[a-f0-9]{40}$")
})
