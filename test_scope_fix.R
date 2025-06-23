library(testthat)
library(withr)
library(neuroarchive)

# Demonstrate the problem
test_that("problem demonstration", {
  cat("\n--- PROBLEM DEMONSTRATION ---\n")
  
  # This will fail
  tmp <- local_tempfile(fileext = ".h5") 
  arr <- matrix(1:5, nrow = 1)
  
  # Function that writes and returns handle
  get_handle <- function(tmp, arr) {
    write_lna(arr, file = tmp, transforms = "delta",
              transform_params = list(delta = list(axis = 1, coding_method = "none")))
    cat("In function - file exists:", file.exists(tmp), "\n")
    read_lna(tmp)
  }
  
  h <- get_handle(tmp, arr)
  cat("After function - file exists:", file.exists(tmp), "\n")
  
  # File still exists because we're still in the test_that scope
  expect_true(file.exists(tmp))
})

# Show the fix - use regular tempfile
test_that("fix with regular tempfile", {
  cat("\n--- FIX WITH REGULAR TEMPFILE ---\n")
  
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)  # Clean up manually
  
  arr <- matrix(1:5, nrow = 1)
  
  get_handle <- function(tmp, arr) {
    write_lna(arr, file = tmp, transforms = "delta",
              transform_params = list(delta = list(axis = 1, coding_method = "none")))
    cat("In function - file exists:", file.exists(tmp), "\n")
    read_lna(tmp)
  }
  
  h <- get_handle(tmp, arr)
  cat("After function - file exists:", file.exists(tmp), "\n")
  
  expect_true(file.exists(tmp))
  expect_equal(h$stash$input, arr)
})

# Show another fix - defer in parent frame
test_that("fix with defer parent", {
  cat("\n--- FIX WITH DEFER PARENT ---\n")
  
  create_temp_in_parent <- function() {
    tmp <- tempfile(fileext = ".h5")
    withr::defer(unlink(tmp), envir = parent.frame())
    tmp
  }
  
  tmp <- create_temp_in_parent()
  arr <- matrix(1:5, nrow = 1)
  
  get_handle <- function(tmp, arr) {
    write_lna(arr, file = tmp, transforms = "delta",
              transform_params = list(delta = list(axis = 1, coding_method = "none")))
    cat("In function - file exists:", file.exists(tmp), "\n")
    read_lna(tmp)
  }
  
  h <- get_handle(tmp, arr)
  cat("After function - file exists:", file.exists(tmp), "\n")
  
  expect_true(file.exists(tmp))
  expect_equal(h$stash$input, arr)
})