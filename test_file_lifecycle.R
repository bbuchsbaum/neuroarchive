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
        entry = function() {
          info <- get("info", parent.frame())
          cat("3. close_output_h5 called with:\n")
          cat("   - file:", info$file, "\n")
          cat("   - in_memory:", info$in_memory, "\n")
          cat("   - File exists before close:", file.exists(info$file), "\n")
        },
        exit = function() {
          info <- get("info", parent.frame(2))
          cat("4. close_output_h5 finished:\n")
          cat("   - File exists after close:", file.exists(info$file), "\n")
        })
  
  # Write data
  arr <- matrix(1:5, nrow = 1)
  res <- write_lna(arr, file = tmp, transforms = "delta",
                   transform_params = list(delta = list(axis = 1,
                                                      coding_method = "none")))
  
  untrace(neuroarchive:::open_output_h5)
  untrace(neuroarchive:::close_output_h5)
  
  cat("\n5. After write_lna:\n")
  cat("   - File exists:", file.exists(tmp), "\n")
  
  # Try to read
  tryCatch({
    h <- read_lna(tmp)
    cat("\n6. read_lna succeeded\n")
    cat("   - File exists:", file.exists(tmp), "\n")
  }, error = function(e) {
    cat("\n6. read_lna failed:", conditionMessage(e), "\n")
    cat("   - File exists:", file.exists(tmp), "\n")
  })
})

# Test with garbage collection
test_that("file lifecycle with gc", {
  tmp <- local_tempfile(fileext = ".h5")
  
  arr <- matrix(1:5, nrow = 1)
  res <- write_lna(arr, file = tmp, transforms = "delta",
                   transform_params = list(delta = list(axis = 1,
                                                      coding_method = "none")))
  
  cat("\nBefore gc - File exists:", file.exists(tmp), "\n")
  gc()
  cat("After gc - File exists:", file.exists(tmp), "\n")
  
  h <- read_lna(tmp)
  expect_equal(h$stash$input, arr)
})