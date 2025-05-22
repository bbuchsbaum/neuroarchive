library(testthat)

# Basic behavior of auto_block_size helper

test_that("auto_block_size reduces slab to target", {
  res <- auto_block_size(c(64, 64, 32), element_size_bytes = 4,
                         target_slab_bytes = 64 * 1024)
  expect_equal(length(res$slab_dims), 3)
  expect_true(all(res$slab_dims <= c(64, 64, 1)))
  expect_true(prod(res$slab_dims) * 4 <= 64 * 1024)
  expect_equal(res$iterate_slabs, ceiling(c(64, 64, 32) / res$slab_dims))
})

