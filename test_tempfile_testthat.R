library(testthat)
library(withr)
library(hdf5r)
library(neuroarchive)

# Test to see how local_tempfile behaves inside testthat
test_that("local_tempfile behavior in testthat", {
  cat("\n=== Testing local_tempfile in testthat context ===\n")
  
  # Test 1: Basic behavior
  tmp <- local_tempfile(fileext = ".h5")
  cat("Created temp file:", tmp, "\n")
  cat("File exists after creation:", file.exists(tmp), "\n")
  
  # Write something
  h5 <- H5File$new(tmp, mode = "w")
  h5$create_group("test")
  h5$close_all()
  
  cat("File exists after writing:", file.exists(tmp), "\n")
  
  # Test 2: Using write_lna and read_lna
  arr <- matrix(1:5, nrow = 1)
  tmp2 <- local_tempfile(fileext = ".h5")
  
  res <- write_lna(arr, file = tmp2, transforms = "delta",
                   transform_params = list(delta = list(axis = 1,
                                                      coding_method = "none")))
  
  cat("File exists after write_lna:", file.exists(tmp2), "\n")
  
  # Direct H5File access
  h5_direct <- H5File$new(tmp2, mode = "r")
  cat("Direct H5File open succeeded\n")
  h5_direct$close_all()
  
  # Try read_lna
  h <- read_lna(tmp2)
  cat("read_lna succeeded\n")
})

# Test the actual check_roundtrip function
check_roundtrip <- function(tmp, arr, axis = -1, coding = "none") {
  write_lna(arr, file = tmp, transforms = "delta",
            transform_params = list(delta = list(axis = axis,
                                               coding_method = coding)))
  h <- read_lna(tmp)
  h # Return the handle itself
}

test_that("check_roundtrip function", {
  cat("\n=== Testing check_roundtrip ===\n")
  
  vec <- 1:5
  tmp <- local_tempfile(fileext = ".h5")
  
  cat("Before check_roundtrip - file exists:", file.exists(tmp), "\n")
  
  # Call check_roundtrip
  h <- check_roundtrip(tmp, vec, axis = 1, coding = "none")
  
  cat("After check_roundtrip - file exists:", file.exists(tmp), "\n")
  cat("Result class:", class(h), "\n")
  
  expect_equal(h$stash$input, vec)
})