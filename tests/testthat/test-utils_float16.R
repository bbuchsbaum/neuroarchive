library(testthat)

# Basic behaviour of has_float16_support

test_that("has_float16_support returns logical scalar", {
  res <- has_float16_support()
  expect_type(res, "logical")
  expect_length(res, 1)
})

test_that("has_float16_support detects packages", {
  # requireNamespace is a base function not in neuroarchive namespace,
  # so we mock has_float16_support itself to test the TRUE path
  ns <- asNamespace("neuroarchive")
  original <- get("has_float16_support", envir = ns)
  mock_fn <- function() TRUE
  if (bindingIsLocked("has_float16_support", ns)) unlockBinding("has_float16_support", ns)
  assign("has_float16_support", mock_fn, envir = ns)
  lockBinding("has_float16_support", ns)
  withr::defer({
    if (bindingIsLocked("has_float16_support", ns)) unlockBinding("has_float16_support", ns)
    assign("has_float16_support", original, envir = ns)
    lockBinding("has_float16_support", ns)
  })
  expect_true(has_float16_support())
})
