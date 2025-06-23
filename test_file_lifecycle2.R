library(testthat)
library(withr)
library(hdf5r)
library(neuroarchive)

# Test file lifecycle
test_that("file lifecycle tracking", {
  tmp <- local_tempfile(fileext = ".h5")
  cat("\n1. Created temp file:", tmp, "\n")
  cat("   File exists:", file.exists(tmp), "\n")
  
  # Track write_lna internals
  trace(neuroarchive:::open_output_h5, 
        exit = function() {
          ret <- returnValue()
          cat("2. open_output_h5 returned:\n")
          cat("   - file:", ret$file, "\n")
          cat("   - in_memory:", ret$in_memory, "\n")
          cat("   - File exists:", file.exists(ret$file), "\n")
        })
  
  trace(neuroarchive:::close_output_h5,
        tracer = function(info) {
          cat("3. close_output_h5 called with:\n")
          cat("   - file:", info$file, "\n")
          cat("   - in_memory:", info$in_memory, "\n")
          cat("   - File exists before close:", file.exists(info$file), "\n")
        })
  
  # Write data
  arr <- matrix(1:5, nrow = 1)
  res <- write_lna(arr, file = tmp, transforms = "delta",
                   transform_params = list(delta = list(axis = 1,
                                                      coding_method = "none")))
  
  untrace(neuroarchive:::open_output_h5)
  untrace(neuroarchive:::close_output_h5)
  
  cat("\n4. After write_lna:\n")
  cat("   - File exists:", file.exists(tmp), "\n")
  
  # Try to read
  tryCatch({
    h <- read_lna(tmp)
    cat("\n5. read_lna succeeded\n")
    cat("   - File exists:", file.exists(tmp), "\n")
  }, error = function(e) {
    cat("\n5. read_lna failed:", conditionMessage(e), "\n")
    cat("   - File exists:", file.exists(tmp), "\n")
  })
})

# Test check_roundtrip directly
check_roundtrip <- function(tmp, arr, axis = -1, coding = "none") {
  cat("\nIn check_roundtrip:\n")
  cat("  Before write_lna - file exists:", file.exists(tmp), "\n")
  
  write_lna(arr, file = tmp, transforms = "delta",
            transform_params = list(delta = list(axis = axis,
                                               coding_method = coding)))
  
  cat("  After write_lna - file exists:", file.exists(tmp), "\n")
  
  h <- read_lna(tmp)
  
  cat("  After read_lna - file exists:", file.exists(tmp), "\n")
  
  h # Return the handle itself
}

test_that("check_roundtrip with monitoring", {
  vec <- 1:5
  tmp <- local_tempfile(fileext = ".h5")
  
  h <- check_roundtrip(tmp, vec, axis = 1, coding = "none")
  expect_equal(h$stash$input, vec)
})