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

test_that("abort_lna creates all error subclasses correctly", {
  # Test all error subclasses
  expect_error(
    abort_lna("file not found", .subclass = "lna_error_missing_data"),
    class = "lna_error_missing_data"
  )
  
  expect_error(
    abort_lna("file exists", .subclass = "lna_error_file_exists"),
    class = "lna_error_file_exists"
  )
  
  expect_error(
    abort_lna("path missing", .subclass = "lna_error_missing_path"),
    class = "lna_error_missing_path"
  )
  
  expect_error(
    abort_lna("dependency missing", .subclass = "lna_error_dependency"),
    class = "lna_error_dependency"
  )
  
  expect_error(
    abort_lna("not implemented", .subclass = "lna_error_not_implemented"),
    class = "lna_error_not_implemented"
  )
  
  expect_error(
    abort_lna("duplicate key", .subclass = "lna_error_duplicate_key"),
    class = "lna_error_duplicate_key"
  )
})

test_that("abort_lna includes location information", {
  err <- expect_error(
    abort_lna("test error", .subclass = "lna_error_io", location = "test_function"),
    class = "lna_error_io"
  )
  
  # Check that location is included in the condition
  expect_true("location" %in% names(err))
  expect_equal(err$location, "test_function")
})

test_that("warn_lna creates warnings with correct subclasses", {
  # Test warning subclasses
  expect_warning(
    warn_lna("performance issue", .subclass = "lna_warning_performance"),
    class = "lna_warning_performance"
  )
  
  expect_warning(
    warn_lna("data quality issue", .subclass = "lna_warning_data_quality"),
    class = "lna_warning_data_quality"
  )
  
  expect_warning(
    warn_lna("overwriting file", .subclass = "lna_warning_overwrite"),
    class = "lna_warning_overwrite"
  )
  
  expect_warning(
    warn_lna("missing attribute", .subclass = "lna_warning_missing_attribute"),
    class = "lna_warning_missing_attribute"
  )
  
  expect_warning(
    warn_lna("parameter adjusted", .subclass = "lna_warning_parameter_adjustment"),
    class = "lna_warning_parameter_adjustment"
  )
  
  expect_warning(
    warn_lna("dependency missing", .subclass = "lna_warning_dependency"),
    class = "lna_warning_dependency"
  )
  
  expect_warning(
    warn_lna("data format issue", .subclass = "lna_warning_data_format"),
    class = "lna_warning_data_format"
  )
})

test_that("warn_lna includes location information", {
  warn <- expect_warning(
    warn_lna("test warning", .subclass = "lna_warning_performance", location = "test_location"),
    class = "lna_warning_performance"
  )
  
  # Check that location is included
  expect_true("location" %in% names(warn))
  expect_equal(warn$location, "test_location")
})

test_that("error handlers work without subclass", {
  # Should still create an lna_error
  expect_error(
    abort_lna("generic error"),
    class = "lna_error"
  )
  
  # Should still create an lna_warning
  expect_warning(
    warn_lna("generic warning"),
    class = "lna_warning"
  )
})
