#!/usr/bin/env Rscript
# Run testthat tests one file at a time to isolate failures
# Usage: Rscript run-tests-serial.R

library(testthat)
# Load development source code instead of installed package
# This ensures tests run against current code with all recent fixes
pkgload::load_all(".", quiet = TRUE)

# Get all test files
test_dir <- "tests/testthat"
test_files <- list.files(test_dir, pattern = "^test-.*\\.R$", full.names = TRUE)

cat(sprintf("\n=== Running %d test files serially ===\n\n", length(test_files)))

results <- list()
failed_files <- character()

for (i in seq_along(test_files)) {
  test_file <- test_files[i]
  test_name <- basename(test_file)

  cat(sprintf("[%d/%d] Testing %s...\n", i, length(test_files), test_name))

  # Force HDF5 cleanup BEFORE each test file
  try({
    hdf5r::h5garbage_collect()
    gc(full = TRUE, verbose = FALSE)
  }, silent = TRUE)
  Sys.sleep(0.1)  # Give OS time to release handles

  # Reload package before each test file to clear any mocks/pollution
  # This ensures test isolation between files
  pkgload::load_all(".", quiet = TRUE, reset = TRUE)

  # Run test file in isolation
  result <- tryCatch({
    test_results <- test_file(test_file, stop_on_failure = FALSE)

    # Force HDF5 cleanup after each file
    try({
      hdf5r::h5garbage_collect()
      gc(full = TRUE, verbose = FALSE)
    }, silent = TRUE)
    Sys.sleep(0.1)  # Give OS time to release handles

    # Count actual passes/failures from test results
    # Each element in test_results is a test with $results containing expectations
    n_failed <- 0
    n_passed <- 0
    n_skipped <- 0
    n_warnings <- 0

    for (test in test_results) {
      if (!is.null(test$results)) {
        for (expectation in test$results) {
          if (inherits(expectation, "expectation_failure") ||
              inherits(expectation, "expectation_error")) {
            n_failed <- n_failed + 1
          } else if (inherits(expectation, "expectation_warning")) {
            n_warnings <- n_warnings + 1
          } else if (inherits(expectation, "expectation_skip")) {
            n_skipped <- n_skipped + 1
          } else if (inherits(expectation, "expectation_success")) {
            n_passed <- n_passed + 1
          }
        }
      }
    }

    if (n_failed > 0) {
      failed_files <<- c(failed_files, test_name)
    }

    list(
      file = test_name,
      passed = n_passed,
      failed = n_failed,
      warnings = n_warnings,
      skipped = n_skipped,
      total = n_passed + n_failed,  # Only count pass/fail expectations, not warnings/skips
      success = n_failed == 0
    )
  }, error = function(e) {
    cat(sprintf("  ERROR in %s: %s\n", test_name, conditionMessage(e)))
    failed_files <<- c(failed_files, test_name)
    list(
      file = test_name,
      passed = 0,
      failed = 1,
      total = 1,
      success = FALSE,
      error = conditionMessage(e)
    )
  })

  results[[i]] <- result

  if (result$success) {
    cat(sprintf("  ✓ %d tests passed\n", result$passed))
  } else {
    cat(sprintf("  ✗ %d passed, %d FAILED\n", result$passed, result$failed))
  }
  cat("\n")
}

# Summary table
cat("\n=== RESULTS TABLE ===\n\n")
cat(sprintf("%-50s %8s %8s %8s %8s %10s\n", "Test File", "Passed", "Failed", "Warns", "Total", "Status"))
cat(paste(rep("=", 100), collapse = ""), "\n")

total_passed <- 0
total_failed <- 0
total_warnings <- 0
files_passed <- 0
files_failed <- 0

for (result in results) {
  status <- if (result$success) "✓ PASS" else "✗ FAIL"
  warns <- if (!is.null(result$warnings)) result$warnings else 0
  cat(sprintf("%-50s %8d %8d %8d %8d %10s\n",
              result$file, result$passed, result$failed, warns, result$total, status))

  total_passed <- total_passed + result$passed
  total_failed <- total_failed + result$failed
  total_warnings <- total_warnings + warns

  if (result$success) {
    files_passed <- files_passed + 1
  } else {
    files_failed <- files_failed + 1
  }
}

cat(paste(rep("=", 100), collapse = ""), "\n")
cat(sprintf("%-50s %8d %8d %8d %8d\n",
            "TOTALS", total_passed, total_failed, total_warnings, total_passed + total_failed))
cat("\n")

cat("=== SUMMARY ===\n")
cat(sprintf("Test files: %d passed, %d failed (out of %d total)\n",
            files_passed, files_failed, length(test_files)))
cat(sprintf("Total tests: %d passed, %d failed (out of %d total)\n",
            total_passed, total_failed, total_passed + total_failed))

if (files_failed > 0) {
  cat("\n❌ FAILED TEST FILES:\n")
  for (f in failed_files) {
    # Find result for this file to show counts
    file_result <- results[[which(sapply(results, function(r) r$file == f))]]
    cat(sprintf("  - %s (%d passed, %d FAILED)\n", f, file_result$passed, file_result$failed))
  }
  cat("\n")
  quit(status = 1)
} else {
  cat("\n✅ All tests passed!\n\n")
  quit(status = 0)
}
