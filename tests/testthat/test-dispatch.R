library(testthat)

# Source functions if not running via devtools::test()
# source("../R/dispatch.R")
# Need a dummy DataHandle class for testing signature
R6::R6Class("DataHandle", list(initialize = function(){}))

test_that("forward_step generic dispatches to default and errors", {
  dummy_type <- "non_existent_type"
  dummy_desc <- list(a = 1)
  dummy_handle <- DataHandle$new()

  # Check that the default method is listed
  expect_true("forward_step.default" %in% methods("forward_step"))

  # Check that calling with an undefined type dispatches to default and errors
  expect_error(
    forward_step(dummy_type, dummy_desc, dummy_handle),
    regexp = paste("No forward_step method implemented for transform type:", dummy_type),
     class = "lna_error_no_method"
  )
})

test_that("invert_step generic dispatches to default and errors", {
  dummy_type <- "another_undefined_type"
  dummy_desc <- list(b = 2)
  dummy_handle <- DataHandle$new()

  # Check that the default method is listed
  expect_true("invert_step.default" %in% methods("invert_step"))

  # Check that calling with an undefined type dispatches to default and errors
  expect_error(
    invert_step(dummy_type, dummy_desc, dummy_handle),

    regexp = paste("No invert_step method implemented for transform type:", dummy_type),
     class = "lna_error_no_method"
  )
}) 
