library(testthat)
library(hdf5r)
library(withr)
library(neuroarchive)

# Ensure core_read and helpers are loaded
# source("../R/core_read.R")

# Helper to create empty transforms group
create_empty_lna <- function(path) {
  h5 <- neuroarchive:::open_h5(path, mode = "w")
  h5$create_group("transforms")
  neuroarchive:::close_h5_safely(h5)
}

# Helper to create lna with one dummy descriptor
create_dummy_lna <- function(path) {
  h5 <- neuroarchive:::open_h5(path, mode = "w")
  tf <- h5$create_group("transforms")
  write_json_descriptor(tf, "00_dummy.json", list(type = "dummy"))
  neuroarchive:::close_h5_safely(h5)
}

# Helper to create lna with a dummy run (for tests needing run_id resolution)
create_lna_with_dummy_run <- function(path) {
  h5 <- neuroarchive:::open_h5(path, mode = "w")
  h5$create_group("transforms") # Minimal /transforms
  scans_group <- h5$create_group("scans") # Create /scans
  scans_group$create_group("run-01")      # Create a dummy run
  neuroarchive:::close_h5_safely(h5)
}

# Helper to create lna with a dummy run and a dummy descriptor
create_lna_with_dummy_run_and_descriptor <- function(path) {
  # First, create the file with a dummy run and /transforms group
  create_lna_with_dummy_run(path)
  
  # Now, open it and add the dummy descriptor
  h5 <- neuroarchive:::open_h5(path, mode = "r+")
  on.exit(neuroarchive:::close_h5_safely(h5), add = TRUE)
  tf_group <- h5[["transforms"]]
  write_json_descriptor(tf_group, "00_dummy.json", list(type = "dummy"))
}

test_that("core_read handles empty /transforms group", {
  tmp <- local_tempfile(fileext = ".h5")
  create_lna_with_dummy_run(tmp)

  handle <- core_read(tmp)
  expect_true(inherits(handle, "DataHandle"))
  expect_false(handle$h5$is_valid)
})

test_that("core_read closes file if invert_step errors", {
  tmp <- local_tempfile(fileext = ".h5")
  create_dummy_lna(tmp) # Creates /transforms/00_dummy.json

  # Also need to ensure a run exists for core_read to proceed
  h5_temp <- neuroarchive:::open_h5(tmp, mode = "r+")
  if (!h5_temp$exists("scans")) {
    h5_temp$create_group("scans")
  }
  if (!h5_temp[["scans"]]$exists("run-01")){
    h5_temp[["scans"]]$create_group("run-01")
  }
  neuroarchive:::close_h5_safely(h5_temp)

  captured_h5 <- NULL

  # Define and locally register the S3 method mock for invert_step.dummy
  mock_invert_step_dummy_closes_file <- function(type, desc, handle) {
    captured_h5 <<- handle$h5
    stop("mock error")
  }
  # Ensure invert_step generic exists for local registration
  if (!exists("invert_step", mode = "function", envir = .GlobalEnv)) {
    .GlobalEnv$invert_step <- function(type, ...) UseMethod("invert_step", type)
    withr::defer(rm(invert_step, envir = .GlobalEnv))
  }
  # Save current global invert_step.dummy if it exists, then assign mock, then defer restoration/removal
  original_invert_step_dummy <- if(exists("invert_step.dummy", envir = .GlobalEnv, inherits = FALSE)) .GlobalEnv$invert_step.dummy else NA
  .GlobalEnv$invert_step.dummy <- mock_invert_step_dummy_closes_file
  if (identical(original_invert_step_dummy, NA)) {
    withr::defer(rm(invert_step.dummy, envir = .GlobalEnv))
  } else {
    withr::defer(assign("invert_step.dummy", original_invert_step_dummy, envir = .GlobalEnv))
  }

  expect_error(core_read(tmp), "mock error")

  expect_true(inherits(captured_h5, "H5File"))
  expect_false(captured_h5$is_valid)
})

test_that("core_read lazy=TRUE keeps file open", {
  tmp <- local_tempfile(fileext = ".h5")
  create_lna_with_dummy_run(tmp)

  handle <- core_read(tmp, lazy = TRUE)
  expect_true(handle$h5$is_valid)
  neuroarchive:::close_h5_safely(handle$h5)
})

test_that("core_read output_dtype stored and float16 check", {
  tmp <- local_tempfile(fileext = ".h5")
  create_lna_with_dummy_run(tmp)

  h <- core_read(tmp, output_dtype = "float64")
  expect_equal(h$meta$output_dtype, "float64")
  err <- expect_error(
    core_read(tmp, output_dtype = "float16"),
    class = "lna_error_float16_unsupported"
  )
  expect_true(grepl("core_read", err$location))
})

test_that("core_read allows float16 when support present", {
  tmp <- local_tempfile(fileext = ".h5")
  create_lna_with_dummy_run(tmp)

  local_mocked_bindings(
    has_float16_support = function() TRUE,
    .env = asNamespace("neuroarchive")
  )
  h <- core_read(tmp, output_dtype = "float16")
  expect_equal(h$meta$output_dtype, "float16")
  expect_false(h$h5$is_valid)
  neuroarchive:::close_h5_safely(h$h5)
})

test_that("core_read works with progress handlers", {
  tmp <- local_tempfile(fileext = ".h5")
  create_lna_with_dummy_run_and_descriptor(tmp)
  
  old_handlers <- progressr::handlers()
  withr::defer(progressr::handlers(old_handlers))
  
  progressr::handlers(progressr::handler_void()) # Set the void handler for the test

  # Define and locally register the S3 method mock for invert_step.dummy
  mock_invert_step_dummy_progress <- function(type, desc, handle) handle # Simple mock
  if (!exists("invert_step", mode = "function", envir = .GlobalEnv)) {
    .GlobalEnv$invert_step <- function(type, ...) UseMethod("invert_step", type)
    withr::defer(rm(invert_step, envir = .GlobalEnv))
  }
  original_isd_progress <- if(exists("invert_step.dummy", envir = .GlobalEnv, inherits = FALSE)) .GlobalEnv$invert_step.dummy else NA
  .GlobalEnv$invert_step.dummy <- mock_invert_step_dummy_progress
  if (identical(original_isd_progress, NA)) {
    withr::defer(rm(invert_step.dummy, envir = .GlobalEnv))
  } else {
    withr::defer(assign("invert_step.dummy", original_isd_progress, envir = .GlobalEnv))
  }
  
  # Capture output to diagnose non-silent behavior
  output <- testthat::capture_output_lines(progressr::with_progress(core_read(tmp)))
  if (length(output) > 0) {
    cat("\nCaptured output in 'core_read works with progress handlers':\n")
    cat(paste("  -", output), sep = "\n")
    cat("\n")
  }
  expect_silent(progressr::with_progress(core_read(tmp)))
})

test_that("core_read validate=TRUE calls validate_lna", {
  tmp <- local_tempfile(fileext = ".h5")
  create_lna_with_dummy_run(tmp)
  called <- FALSE
  local_mocked_bindings(
    validate_lna = function(file) { called <<- TRUE },
    .env = asNamespace("neuroarchive")
  )
  core_read(tmp, validate = TRUE)
  expect_true(called)
})

test_that("core_read allow_plugins='none' errors on unknown transform", {
  tmp <- local_tempfile(fileext = ".h5")
  create_lna_with_dummy_run_and_descriptor(tmp)
  expect_error(core_read(tmp, allow_plugins = "none"),
               class = "lna_error_no_method")
})

test_that("core_read allow_plugins='prompt' falls back when non-interactive", {
  tmp <- local_tempfile(fileext = ".h5")
  create_lna_with_dummy_run_and_descriptor(tmp)
  
  # Mock non-interactive environment
  local_mocked_bindings(
    is_interactive = function() FALSE,
    .package = "rlang"
  )
  
  # Should warn about missing method but continue
  expect_warning(
    core_read(tmp, allow_plugins = "prompt"),
    "Missing invert_step"
  )
})

test_that("core_read allow_plugins='prompt' interactive respects user choice", {
  tmp <- local_tempfile(fileext = ".h5")
  create_lna_with_dummy_run_and_descriptor(tmp)
  
  # Test with user answering "no"
  local_mocked_bindings(
    is_interactive = function() TRUE,
    .package = "rlang"
  )
  
  local_mocked_bindings(
    readline = function(prompt) "n",
    .package = "base"
  )
  
  # Should error since user declined
  expect_error(
    core_read(tmp, allow_plugins = "prompt"),
    class = "lna_error_no_method"
  )
  
  # Test with user answering "yes"
  local_mocked_bindings(
    is_interactive = function() TRUE,
    .package = "rlang"
  )
  
  local_mocked_bindings(
    readline = function(prompt) "y",
    .package = "base"
  )
  
  # Should warn about missing method but continue
  expect_warning(
    core_read(tmp, allow_plugins = "prompt"),
    "Missing invert_step"
  )
})

test_that("core_read stores subset parameters", {
  tmp <- local_tempfile(fileext = ".h5")
  create_lna_with_dummy_run(tmp) # Use helper that creates a run
  roi <- array(TRUE, dim = c(1,1,1))
  h <- core_read(tmp, roi_mask = roi, time_idx = 1:2)
  expect_identical(h$subset$roi_mask, roi)
  expect_identical(h$subset$time_idx, 1:2)
  expect_false(h$h5$is_valid)
})

test_that("core_read run_id globbing returns handles", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- neuroarchive:::open_h5(tmp, mode = "w")
  plan <- Plan$new()
  plan$add_descriptor("00_dummy.json", list(type = "dummy"))
  plan$add_payload("p1", 1L)
  plan$add_dataset_def("/scans/run-01/data", "data", "dummy", "run-01", 0L, "{}", "p1", "eager", dtype = NA_character_)
  plan$add_payload("p2", 2L)
  plan$add_dataset_def("/scans/run-02/data", "data", "dummy", "run-02", 0L, "{}", "p2", "eager", dtype = NA_character_)
  materialise_plan(h5, plan)
  neuroarchive:::close_h5_safely(h5)

  mock_invert_step_dummy_glob <- function(type, desc, handle) {
    # Simple mock, can add more sophisticated behavior if needed for other tests
    return(handle)
  }

  # Ensure invert_step generic exists
  if (!exists("invert_step", mode = "function", envir = .GlobalEnv)) {
    .GlobalEnv$invert_step <- function(type, ...) UseMethod("invert_step", type)
    withr::defer(rm(invert_step, envir = .GlobalEnv))
  }

  original_isd_glob <- if(exists("invert_step.dummy", envir = .GlobalEnv, inherits = FALSE)) {
    .GlobalEnv$invert_step.dummy
  } else {
    NA # Sentinel
  }

  .GlobalEnv$invert_step.dummy <- mock_invert_step_dummy_glob

  if (identical(original_isd_glob, NA)) {
    withr::defer(rm(invert_step.dummy, envir = .GlobalEnv))
  } else {
    withr::defer(assign("invert_step.dummy", original_isd_glob, envir = .GlobalEnv))
  }

  handles <- core_read(tmp, run_id = "run-*")
  expect_length(handles, 2)
  expect_true(all(sapply(handles, inherits, "DataHandle")))
})

test_that("core_read run_id globbing lazy returns first", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- neuroarchive:::open_h5(tmp, mode = "w")
  plan <- Plan$new()
  plan$add_descriptor("00_dummy.json", list(type = "dummy"))
  plan$add_payload("p1", 1L)
  plan$add_dataset_def("/scans/run-01/data", "data", "dummy", "run-01", 0L, "{}", "p1", "eager", dtype = NA_character_)
  plan$add_payload("p2", 2L)
  plan$add_dataset_def("/scans/run-02/data", "data", "dummy", "run-02", 0L, "{}", "p2", "eager", dtype = NA_character_)
  materialise_plan(h5, plan)
  neuroarchive:::close_h5_safely(h5)

  mock_invert_step_dummy_glob_lazy <- function(type, desc, handle) {
    return(handle)
  }

  if (!exists("invert_step", mode = "function", envir = .GlobalEnv)) {
    .GlobalEnv$invert_step <- function(type, ...) UseMethod("invert_step", type)
    withr::defer(rm(invert_step, envir = .GlobalEnv))
  }

  original_isd_lazy <- if(exists("invert_step.dummy", envir = .GlobalEnv, inherits = FALSE)) {
    .GlobalEnv$invert_step.dummy
  } else {
    NA # Sentinel
  }

  .GlobalEnv$invert_step.dummy <- mock_invert_step_dummy_glob_lazy

  if (identical(original_isd_lazy, NA)) {
    withr::defer(rm(invert_step.dummy, envir = .GlobalEnv))
  } else {
    withr::defer(assign("invert_step.dummy", original_isd_lazy, envir = .GlobalEnv))
  }

  # For lazy globbing with multiple matches, core_read issues a warning and processes only the first.
  expect_warning(
    handle <- core_read(tmp, run_id = "run-*", lazy = TRUE),
    "Multiple runs matched; using first match in lazy mode"
  )
  expect_true(inherits(handle, "DataHandle"))
  expect_true(handle$h5$is_valid) # File should be open
  neuroarchive:::close_h5_safely(handle$h5)
})

test_that("core_read validate=TRUE checks dataset existence", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- neuroarchive:::open_h5(tmp, mode = "w")
  tf <- h5$create_group("transforms")
  desc <- list(type = "dummy",
               datasets = list(list(path = "/scans/run-01/missing", role = "data")))
  write_json_descriptor(tf, "00_dummy.json", desc)

  # Add root attributes to pass validate_lna()
  root <- h5[["/"]]
  neuroarchive:::h5_attr_write(root, "lna_spec", "LNA R v2.0")
  neuroarchive:::h5_attr_write(root, "creator", "neuroarchive test")
  neuroarchive:::h5_attr_write(root, "required_transforms", "dummy")

  neuroarchive:::close_h5_safely(h5)

  # Ensure the run group itself exists for core_read to find the run
  # before it checks for the dataset within the run.
  h5_temp <- neuroarchive:::open_h5(tmp, mode = "r+")
  if (!h5_temp$exists("scans")) {
    h5_temp$create_group("scans")
  }
  if (!h5_temp[["scans"]]$exists("run-01")){
    h5_temp[["scans"]]$create_group("run-01")
  }
  neuroarchive:::close_h5_safely(h5_temp)

  # Define and locally register the S3 method mock for invert_step.dummy
  mock_invert_step_dummy_validate <- function(type, desc, handle) handle

  # Ensure invert_step generic exists for local registration if not already present
  # This global manipulation is not ideal for tests but reflects current structure.
  generic_created_here <- FALSE
  if (!exists("invert_step", mode = "function", envir = .GlobalEnv)) {
    .GlobalEnv$invert_step <- function(type, ...) UseMethod("invert_step", type)
    withr::defer(rm(invert_step, envir = .GlobalEnv))
    generic_created_here <- TRUE
  }

  # Assign a placeholder for invert_step.dummy before mocking
  .GlobalEnv$invert_step.dummy <- function(...) stop("Placeholder, should be mocked")
  withr::defer(rm(invert_step.dummy, envir = .GlobalEnv))
  testthat::local_mocked_bindings(
    invert_step.dummy = mock_invert_step_dummy_validate,
    .env = .GlobalEnv # Mocking the S3 method in the global environment
  )
  
  # If the generic was created here, and no S3 method for dummy existed before,
  # local_mocked_bindings might not clean up the .GlobalEnv$invert_step.dummy if the generic itself is removed.
  # However, testthat's cleanup of bindings in .GlobalEnv should be robust.
  # If issues persist with the dynamic generic, a more involved setup might be needed
  # or the generic should be part of the package.

  expect_error(core_read(tmp, validate = TRUE),
               regexp = "HDF5 path '/scans/run-01/missing' not found")

  expect_no_error(core_read(tmp, validate = FALSE))
})

test_that("core_read validate=TRUE checks required params", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- neuroarchive:::open_h5(tmp, mode = "w")
  tf <- h5$create_group("transforms")
  desc <- list(
    type = "embed",
    datasets = list(
      list(path = "/basis/B", role = "basis_matrix"),
      list(path = "/scans/run-01/coeff", role = "coefficients")
    ),
    params = list(
      basis_path = "/basis/B",
      coeff_path = "/scans/run-01/coeff"
    )
  )
  write_json_descriptor(tf, "00_embed.json", desc)
  root <- h5[["/"]]
  h5_write_dataset(root, "/basis/B", matrix(1))
  h5_write_dataset(root, "/scans/run-01/coeff", matrix(1))
  
  # Also ensure the run group /scans/run-01 exists for core_read
  if (!h5$exists("scans/run-01")) { # Check relative to root, or ensure scans exists first
    if (!h5$exists("scans")) {
        h5$create_group("scans")
    }
    h5[["scans"]]$create_group("run-01")
  }
  
  # Add root attributes to pass validate_lna()
  # Ensure root is correctly referenced from the still open h5 handle
  neuroarchive:::h5_attr_write(h5[["/"]], "lna_spec", "LNA R v2.0")
  neuroarchive:::h5_attr_write(h5[["/"]], "creator", "neuroarchive test")
  neuroarchive:::h5_attr_write(h5[["/"]], "required_transforms", "embed")
  
  neuroarchive:::close_h5_safely(h5)

  # Define and locally register the S3 method mock for invert_step.embed
  mock_invert_step_embed_params <- function(type, desc, handle) handle # Simple mock
  # Ensure invert_step generic exists for local registration
  if (!exists("invert_step", mode = "function", envir = .GlobalEnv)) {
    .GlobalEnv$invert_step <- function(type, ...) UseMethod("invert_step", type)
    withr::defer(rm(invert_step, envir = .GlobalEnv))
  }
  # Save current global invert_step.embed if it exists, then assign mock, then defer restoration/removal
  original_invert_step_embed <- if(exists("invert_step.embed", envir = .GlobalEnv, inherits = FALSE)) .GlobalEnv$invert_step.embed else NA
  .GlobalEnv$invert_step.embed <- mock_invert_step_embed_params
  if (identical(original_invert_step_embed, NA)) {
    withr::defer(rm(invert_step.embed, envir = .GlobalEnv))
  } else {
    withr::defer(assign("invert_step.embed", original_invert_step_embed, envir = .GlobalEnv))
  }

  # Mock jsonvalidate::json_validate to prevent verbose output during this test
  testthat::local_mocked_bindings(
    json_validate = function(json, schema, verbose = FALSE, greedy = FALSE, engine = c("ajv", "is-my-json-valid"), ...) {
      # This mock assumes that for the purpose of this specific test (checking parameter validation silence),
      # schema validation itself is not under test and would pass silently.
      # It effectively makes the verbose = TRUE in validate_lna have no output effect.
      return(TRUE)
    },
    .package = "jsonvalidate"
  )
  expect_no_error(core_read(tmp, validate = TRUE))
})

test_that("core_read reports step provenance from failing method", {
  tmp <- local_tempfile(fileext = ".h5")
  
  # Setup a file that core_read can process up to the transform step
  # and includes a "fail" transform descriptor.
  h5_failing <- neuroarchive:::open_h5(tmp, mode = "w")
  root_failing <- h5_failing[["/"]]
  neuroarchive:::h5_attr_write(root_failing, "lna_spec", "LNA R v2.0")
  neuroarchive:::h5_attr_write(root_failing, "creator", "neuroarchive test provenance")
  neuroarchive:::h5_attr_write(root_failing, "required_transforms", "fail") # Or character(0) if "fail" is not globally required
  
  tf_group_failing <- h5_failing$create_group("transforms")
  # Assuming "00_fail.json" will result in step_idx 0 in the reversed loop
  neuroarchive:::write_json_descriptor(tf_group_failing, "00_fail.json", list(type = "fail"))
  
  scans_group_failing <- h5_failing$create_group("scans")
  scans_group_failing$create_group("run-01") # Needs at least one run
  if (!h5_failing$exists("basis")) { # Ensure basis group exists for validate_lna if called
      h5_failing$create_group("basis")
  }
  neuroarchive:::close_h5_safely(h5_failing)

  # Define an invert_step.fail S3 method that itself throws an error
  mock_invert_step_fail_throws_error <- function(type, desc, handle) {
    stop("Intentional error from inside invert_step.fail")
  }

  # Ensure invert_step generic exists for S3 dispatch
  # This global manipulation is not ideal but matches patterns elsewhere in the file.
  if (!exists("invert_step", mode = "function", envir = .GlobalEnv)) {
    .GlobalEnv$invert_step <- function(type, ...) UseMethod("invert_step", type)
    # Defer its removal within the testthat environment
    withr::defer(rm(invert_step, envir = .GlobalEnv))
  }

  # Assign a placeholder for invert_step.fail before mocking
  .GlobalEnv$invert_step.fail <- function(...) stop("Placeholder, should be mocked")
  withr::defer(rm(invert_step.fail, envir = .GlobalEnv))
  # Mock the .fail method using testthat::local_mocked_bindings
  testthat::local_mocked_bindings(
    invert_step.fail = mock_invert_step_fail_throws_error,
    .env = .GlobalEnv # Mocking S3 method in .GlobalEnv
  )

  # core_read should discover and attempt to run invert_step.fail, which then errors.
  # The default allow_plugins = "installed" is fine because the method *is* found.
  err_cr_prov <- expect_error(core_read(tmp))

  cat("\nDEBUG - Actual location for core_read reports step provenance from failing method:", err_cr_prov$location, "\n\n")
  # Check that the error object has the location field set by run_transform_step
  # The error message "Intentional error from inside invert_step.fail" will be part of the wrapped message.
  # The location should be "invert_step:fail[0]" (assuming index 0 for the "00_fail.json" descriptor).
  # We need to escape [ and ] for grepl if not using fixed = TRUE.
  # Using fixed = FALSE (default) for regexp behavior with brackets.
  expect_true(grepl("invert_step:fail[0]", trimws(err_cr_prov$location), fixed = TRUE))
  expect_s3_class(err_cr_prov, "lna_error_transform_step") # This is the class abort_lna in run_transform_step uses
})
