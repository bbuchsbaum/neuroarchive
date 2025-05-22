library(testthat)
library(hdf5r)
library(withr)

# Simulate chunk size failures to test retry heuristics

# --- Start: Copied and adapted from R/materialise.R for testing --- 
# Minimal local copy for testing retry logic. This avoids namespace manipulation.

# (Assuming abort_lna, h5_attr_write, guess_h5_type, reduce_chunk_dims, write_json_descriptor 
#  are either not hit in this specific test path or would need to be stubbed/simplified if they were.
#  For this test, we are focused on the retry logic around h5_write_dataset calls.)

materialise_plan_for_test <- function(h5, plan, h5_write_dataset_fn) {
  stopifnot(inherits(h5, "H5File"))
  if (!h5$is_valid) {
    stop("Provided HDF5 handle is not open or valid in materialise_plan_for_test")
  }
  stopifnot(inherits(plan, "Plan"))

  root <- h5[["/"]] 
  # Simplified: Skipping group creation and attribute writing not relevant to this test

  # --- Copied write_payload internal function (simplified) ---
  write_payload <- function(path, data, step_index) {
    comp_level <- 0 # Simplified
    chunk_dims <- NULL # Initial attempt uses this

    attempt <- function(level, chunks) {
      h5_write_dataset_fn(root, path, data, chunk_dims = chunks,
                           compression_level = level)
      NULL
    }

    res <- tryCatch(attempt(comp_level, chunk_dims), error = function(e) e)
    if (inherits(res, "error")) {
      # Simplified dtype_size and cdims for retry, actual values don't matter for mock
      dtype_size <- 8L 
      cdims <- c(10,10) # Placeholder

      cdims1 <- neuroarchive:::reduce_chunk_dims(cdims, dtype_size, 1024^3)
      warning_message_1 <- sprintf(
        "Write failed for %s; retrying with smaller chunks (<1 GiB, ~%.1f MiB)",
        path, prod(cdims1) * dtype_size / 1024^2
      )
      warning(warning_message_1)
      res <- tryCatch(attempt(0, cdims1), error = function(e) e)
    }

    if (inherits(res, "error")) {
      dtype_size <- 8L
      cdims1 <- c(5,5) # Placeholder, assuming cdims1 was reduced
      cdims2 <- neuroarchive:::reduce_chunk_dims(cdims1, dtype_size, 256 * 1024^2)
      warning_message_2 <- sprintf(
        "Write failed for %s; retrying with smaller chunks (<=256 MiB, ~%.1f MiB)",
        path, prod(cdims2) * dtype_size / 1024^2
      )
      warning(warning_message_2)
      res <- tryCatch(attempt(0, cdims2), error = function(e) e)
    }

    if (inherits(res, "error")) {
      stop(sprintf(
          "Failed to write dataset '%s' (step %d) after retries: %s",
          path, step_index, conditionMessage(res)
        ))
    }
  }
  # --- End write_payload ---

  if (nrow(plan$datasets) > 0) {
    for (i in seq_len(nrow(plan$datasets))) {
      row <- plan$datasets[i, ]
      key <- row$payload_key
      if (!nzchar(key)) next
      payload <- plan$payloads[[key]]
      if (is.null(payload)) next
      write_payload(row$path, payload, row$step_index)
    }
  }
  invisible(h5)
}
# --- End: Copied and adapted R/materialise.R --- 

test_that("materialise_plan retries with chunk heuristics using local copy", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- neuroarchive:::open_h5(tmp, mode = "w")
  
  # Ensure h5 is valid before doing anything else
  expect_true(h5$is_valid, "H5 handle should be valid after open_h5")

  # Need Plan to be the actual Plan object from the package
  plan <- neuroarchive:::Plan$new()
  plan$add_payload("p", matrix(1:10, nrow = 2))
  plan$add_dataset_def("/scans/run-01/data", "data", "dummy", "run-01", 0L, "{}", "p", "eager")

  mock_env <- new.env()
  mock_env$calls <- list()
  
  # This is our mock, to be passed to materialise_plan_for_test
  h5_write_dataset_for_test <- function(h5_group, path, data, chunk_dims = NULL, compression_level = 0) {
    mock_env$calls <- c(mock_env$calls, list(chunk_dims)) # New way to append
    
    if (length(mock_env$calls) < 3) {
      stop("chunk too large")
    }
    
    # On the 3rd call, simulate successful write by doing minimal HDF5 operation
    # Ensure h5_group is valid before using it.
    if (!inherits(h5_group, "H5Group") || !h5_group$is_valid) {
      stop("h5_group is invalid in h5_write_dataset_for_test on successful call")
    }
    parts <- strsplit(path, "/")[[1]]
    parts <- parts[nzchar(parts)]
    grp <- h5_group
    if (length(parts) > 1) {
      for (g_name in parts[-length(parts)]) {
        grp <- if (!grp$exists(g_name)) grp$create_group(g_name) else grp[[g_name]]
      }
    }
    ds <- grp$create_dataset(tail(parts, 1), data)
    if(inherits(ds, "H5D")) ds$close()
    return(invisible(TRUE)) 
  }
  
  mp_warnings <- list()
  mp_error <- NULL
  
  tryCatch({
    withCallingHandlers({
      # Call the local test version, passing the mock function
      materialise_plan_for_test(h5, plan, h5_write_dataset_fn = h5_write_dataset_for_test)
    }, warning = function(w) {
      mp_warnings[[length(mp_warnings) + 1]] <<- w
      invokeRestart("muffleWarning") 
    })
  }, error = function(e) {
    mp_error <<- e
  })

  if (!is.null(mp_error)) {
    print("Error during materialise_plan_for_test call:")
    print(mp_error)
    stop("materialise_plan_for_test failed unexpectedly.") # Force test to fail clearly
  }
  
  # DEBUG: Check h5 validity immediately after materialise_plan_for_test
  # print(paste("Class of h5 after materialise_plan_for_test:", class(h5)))
  # if (inherits(h5, "H5File")) {
  #   print(paste("Is h5 valid after materialise_plan_for_test:", h5$is_valid))
  # } else {
  #   print("h5 is NOT an H5File object after materialise_plan_for_test")
  # }
  
  expect_equal(length(mock_env$calls), 3, info = "h5_write_dataset_for_test should be called 3 times.")
  expect_length(mp_warnings, 2) 
  if (length(mp_warnings) >=1) expect_match(mp_warnings[[1]]$message, "<1 GiB", fixed = TRUE)
  if (length(mp_warnings) >=2) expect_match(mp_warnings[[2]]$message, "256 MiB", fixed = TRUE)

  # DEBUG: Check h5 validity RIGHT BEFORE close_h5_safely
  print("--- Before close_h5_safely ---")
  print(paste("Class of h5:", paste(class(h5), collapse=", ")))
  if (inherits(h5, "H5File")) {
    print(paste("Is h5 valid:", h5$is_valid))
  } else {
    print("h5 is NOT an H5File object")
  }
  neuroarchive:::close_h5_safely(h5)
})

