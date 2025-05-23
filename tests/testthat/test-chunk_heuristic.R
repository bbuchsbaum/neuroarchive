library(testthat)
library(neuroarchive)
# Tests for guess_chunk_dims heuristic

test_that("guess_chunk_dims targets ~1MiB", {
  dims <- c(100, 100, 10)
  res <- guess_chunk_dims(dims, 8L)
  expect_equal(length(res), length(dims))
  expect_true(all(res <= dims))
  expect_lt(prod(res) * 8L, 1.1 * 1024^2)
})

test_that("guess_chunk_dims limits chunks for large data", {
  dims <- c(30000, 20000) # >4 GiB for double
  res <- guess_chunk_dims(dims, 8L)
  expect_lt(prod(res) * 8L, 1024^3)
})
