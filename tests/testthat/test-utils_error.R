library(testthat)

# Basic tests for abort_lna helper and error classes

test_that("abort_lna creates lna_error_io", {
  expect_error(
    abort_lna("io fail", .subclass = "lna_error_io"),
    class = "lna_error_io"
  )
})

test_that("abort_lna creates lna_error_validation", {
  expect_error(
    abort_lna("bad input", .subclass = "lna_error_validation"),
    class = "lna_error_validation"
  )
})
