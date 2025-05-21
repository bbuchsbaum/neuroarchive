library(testthat)

# Basic behaviour of has_float16_support

test_that("has_float16_support returns logical scalar", {
  res <- has_float16_support()
  expect_type(res, "logical")
  expect_length(res, 1)
})

test_that("has_float16_support detects packages", {
  with_mocked_bindings(
    requireNamespace = function(pkg, quietly = TRUE) TRUE,
    {
      expect_true(has_float16_support())
    }
  )
})
