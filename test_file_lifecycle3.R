library(testthat)
library(withr)
library(hdf5r)
library(neuroarchive)

# Simple test without tracing
test_that("file lifecycle simple", {
  # Test 1: Regular tempfile
  tmp1 <- tempfile(fileext = ".h5")
  arr <- matrix(1:5, nrow = 1)
  
  write_lna(arr, file = tmp1, transforms = "delta",
            transform_params = list(delta = list(axis = 1, coding_method = "none")))
  
  cat("\nRegular tempfile:\n")
  cat("  File exists after write:", file.exists(tmp1), "\n")
  
  h1 <- read_lna(tmp1)
  cat("  Read succeeded\n")
  
  # Test 2: local_tempfile
  tmp2 <- local_tempfile(fileext = ".h5")
  
  write_lna(arr, file = tmp2, transforms = "delta",
            transform_params = list(delta = list(axis = 1, coding_method = "none")))
  
  cat("\nlocal_tempfile:\n")
  cat("  File exists after write:", file.exists(tmp2), "\n")
  
  h2 <- read_lna(tmp2)
  cat("  Read succeeded\n")
  
  # Clean up
  unlink(tmp1)
})

# Test timing issues
test_that("file timing issues", {
  tmp <- local_tempfile(fileext = ".h5")
  arr <- matrix(1:5, nrow = 1)
  
  # Write the file
  write_lna(arr, file = tmp, transforms = "delta",
            transform_params = list(delta = list(axis = 1, coding_method = "none")))
  
  # Check file multiple times
  for (i in 1:5) {
    exists <- file.exists(tmp)
    cat("Check", i, "- File exists:", exists, "\n")
    if (!exists) {
      cat("  File disappeared!\n")
      break
    }
    Sys.sleep(0.1)
  }
  
  # Try to read
  if (file.exists(tmp)) {
    h <- read_lna(tmp)
    expect_equal(h$stash$input, arr)
  }
})

# Test environment and scope
test_that("environment scope test", {
  # Create a function that mimics the test scenario
  run_in_scope <- function() {
    tmp <- local_tempfile(fileext = ".h5")
    arr <- matrix(1:5, nrow = 1)
    
    write_lna(arr, file = tmp, transforms = "delta",
              transform_params = list(delta = list(axis = 1, coding_method = "none")))
    
    # Return both the path and whether it exists
    list(path = tmp, exists_after_write = file.exists(tmp))
  }
  
  result <- run_in_scope()
  cat("\nAfter function return:\n")
  cat("  Path:", result$path, "\n")
  cat("  Existed after write:", result$exists_after_write, "\n")
  cat("  Exists now:", file.exists(result$path), "\n")
})

# Test the actual test file pattern
check_roundtrip <- function(tmp, arr, axis = -1, coding = "none") {
  write_lna(arr, file = tmp, transforms = "delta",
            transform_params = list(delta = list(axis = axis,
                                               coding_method = coding)))
  h <- read_lna(tmp)
  h # Return the handle itself
}

test_that("delta handles 1D input", {
  vec <- 1:5
  for (coding in c("none", "rle")) {
    tmp <- local_tempfile(fileext = ".h5")
    h <- check_roundtrip(tmp, vec, axis = 1, coding = coding)
    expect_equal(h$stash$input, vec)
  }
})