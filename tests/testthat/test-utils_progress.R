library(testthat)
library(progressr)

# Tests for is_progress_globally_enabled

test_that("is_progress_globally_enabled reflects active handlers", {
  old_handlers <- progressr::handlers()
  withr::defer(progressr::handlers(old_handlers))

  progressr::handlers(progressr::handler_void())
  expect_false(neuroarchive:::is_progress_globally_enabled())

  progressr::handlers(progressr::handler_txtprogressbar())
  expect_true(neuroarchive:::is_progress_globally_enabled())
})
