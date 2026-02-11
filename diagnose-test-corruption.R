#!/usr/bin/env Rscript
# Diagnose which test corrupts HDF5 state
# Instruments tests #140-160 with detailed logging

library(testthat)
library(neuroarchive)

test_dir <- "tests/testthat"

# Helper to get HDF5 state info
get_h5_state <- function() {
  state <- list(
    timestamp = Sys.time(),
    version = tryCatch(hdf5r::h5version(), error = function(e) "unknown"),
    temp_h5_files = length(list.files(tempdir(), pattern = "\\.h5$"))
  )
  state
}

# Test that HDF5 is working
test_h5_health <- function(label) {
  tmp <- tempfile(fileext = ".h5")
  result <- tryCatch({
    # Try to create and immediately close a simple H5 file
    h5 <- hdf5r::H5File$new(tmp, mode = "w")
    h5$close_all()
    file.exists(tmp)
  }, error = function(e) {
    FALSE
  }, finally = {
    if (file.exists(tmp)) unlink(tmp)
  })

  cat(sprintf("[H5-HEALTH %s] Can create/close H5 files: %s\n",
              label, if (result) "YES" else "NO"))
  result
}

cat("\n=== HDF5 State Corruption Diagnostic ===\n\n")

# Initial health check
cat("Initial H5 health check:\n")
test_h5_health("INITIAL")
state_init <- get_h5_state()
cat(sprintf("  Version: %s\n", state_init$version))
cat(sprintf("  Temp .h5 files: %d\n\n", state_init$temp_h5_files))

# Run full suite up to test #140
cat("Running tests up to #140...\n")
result_140 <- test_dir(test_dir, stop_on_failure = FALSE)
summary_140 <- summary(result_140)
cat(sprintf("Tests run: %d, Passed: %d, Failed: %d\n",
            summary_140$n, summary_140$passed, summary_140$failed))

# Health check after #140
cat("\nH5 health check after first 140 tests:\n")
test_h5_health("AFTER-140")
hdf5r::h5garbage_collect()
Sys.sleep(0.1)

# Now run tests one-at-a-time from #141-160 with detailed logging
cat("\n=== Running tests #141-160 individually ===\n\n")

test_files <- list.files(test_dir, pattern = "^test-.*\\.R$", full.names = TRUE)

for (i in 141:min(160, length(test_files))) {
  test_file <- test_files[i]
  test_name <- basename(test_file)

  cat(sprintf("\n--- Test #%d: %s ---\n", i, test_name))

  # Pre-test health
  pre_health <- test_h5_health(sprintf("PRE-%d", i))
  pre_state <- get_h5_state()

  # Run test
  test_start <- Sys.time()
  test_result <- tryCatch({
    test_file(test_file, stop_on_failure = FALSE)
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", conditionMessage(e)))
    NULL
  })
  test_duration <- as.numeric(difftime(Sys.time(), test_start, units = "secs"))

  # Post-test health
  post_health <- test_h5_health(sprintf("POST-%d", i))
  post_state <- get_h5_state()

  # Force cleanup
  hdf5r::h5garbage_collect()
  Sys.sleep(0.05)

  # Check for health degradation
  if (pre_health && !post_health) {
    cat(sprintf("\n*** CULPRIT FOUND: Test #%d (%s) ***\n", i, test_name))
    cat("HDF5 state was healthy before this test but corrupted after!\n")
    cat(sprintf("Pre-test temp files: %d\n", pre_state$temp_h5_files))
    cat(sprintf("Post-test temp files: %d\n", post_state$temp_h5_files))
    cat(sprintf("Test duration: %.2f seconds\n", test_duration))

    # Try to recover
    cat("\nAttempting recovery...\n")
    gc()
    hdf5r::h5garbage_collect()
    Sys.sleep(0.5)
    recovery_health <- test_h5_health("RECOVERY")
    if (recovery_health) {
      cat("Recovery successful\n")
    } else {
      cat("Recovery failed - HDF5 state permanently corrupted\n")
      break
    }
  }

  cat(sprintf("Duration: %.2f sec, Temp files: %d -> %d\n",
              test_duration,
              pre_state$temp_h5_files,
              post_state$temp_h5_files))
}

cat("\n=== Diagnostic Complete ===\n")
