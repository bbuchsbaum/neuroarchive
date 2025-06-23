test_that("transform_min_dims returns correct dimensions", {
  # Test known transforms
  expect_equal(transform_min_dims("delta"), 1L)
  expect_equal(transform_min_dims("quant"), 1L)
  expect_equal(transform_min_dims("basis"), 2L)
  expect_equal(transform_min_dims("embed"), 2L)
  
  # Test default case
  expect_equal(transform_min_dims("unknown_transform"), 3L)
  expect_equal(transform_min_dims("spat_hrbf"), 3L)
  expect_equal(transform_min_dims("temporal"), 3L)
  
  # Test with empty string
  expect_equal(transform_min_dims(""), 3L)
  
  # Test case sensitivity
  expect_equal(transform_min_dims("DELTA"), 3L)  # Should return default
  expect_equal(transform_min_dims("Delta"), 3L)  # Should return default
})

test_that("transform_min_dims handles edge cases", {
  # Test with NULL should return default
  expect_equal(transform_min_dims(NULL), 3L)
  
  # Test with numeric input should still work (coerced to character)
  expect_equal(transform_min_dims(123), 3L)
  
  # Test with multiple inputs (only first is used)
  expect_equal(transform_min_dims(c("delta", "quant")), 1L)
})