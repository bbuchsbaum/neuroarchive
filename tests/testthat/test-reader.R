library(testthat)
library(withr)
library(hdf5r)
library(neuroarchive)

# Helper to create simple LNA file with no transforms
create_empty_lna <- function(path) {
  h5 <- neuroarchive:::open_h5(path, mode = "w")
  h5$create_group("transforms")
  scans_group <- h5$create_group("scans")
  scans_group$create_group("run-01") # Add a dummy run
  if (inherits(h5, "H5File") && h5$is_valid) {
    tryCatch(h5$close_all(), error = function(e) {
      warning(paste("Error closing HDF5 handle (inlined create_empty_lna):", conditionMessage(e)))
    })
  }
}

# Helper to create lna with one dummy descriptor
create_dummy_lna <- function(path) {
  h5 <- neuroarchive:::open_h5(path, mode = "w")
  tf <- h5$create_group("transforms")
  scans_group <- h5$create_group("scans")
  scans_group$create_group("run-01") # Add a dummy run
  neuroarchive:::write_json_descriptor(tf, "00_dummy.json", list(type = "dummy"))
  if (inherits(h5, "H5File") && h5$is_valid) {
    tryCatch(h5$close_all(), error = function(e) {
      warning(paste("Error closing HDF5 handle (inlined create_dummy_lna):", conditionMessage(e)))
    })
  }
}


test_that("read_lna(lazy=TRUE) returns lna_reader", {
  tmp <- local_tempfile(fileext = ".h5")
  create_empty_lna(tmp)

  reader <- read_lna(tmp, lazy = TRUE)
  expect_s3_class(reader, "lna_reader")
  expect_true(reader$h5$is_valid)
  reader$close()
})

test_that("lna_reader initialize closes file on failure", {
  tmp <- local_tempfile(fileext = ".h5")
  create_empty_lna(tmp)
  
  # Original function from neuroarchive namespace
  original_open_h5 <- getFromNamespace("open_h5", "neuroarchive")
  captured_h5 <- NULL
  
  # Mock implementation
  mock_open_h5 <- function(file, mode = "r") {
    #message(paste("Mock open_h5 called for file:", file, "with mode:", mode))
    # Call the original function to actually open the file and capture the handle
    h5_obj <- original_open_h5(file, mode)
    captured_h5 <<- h5_obj # Capture the handle
    #message(paste("Captured h5 class:", class(captured_h5), "is_valid:", if(inherits(captured_h5, "H5File")) captured_h5$is_valid else "NA"))
    return(h5_obj) # Return the handle as normal
  }
  
  # Temporarily replace open_h5 in neuroarchive's namespace
  unlockBinding("open_h5", asNamespace("neuroarchive"))
  assignInNamespace("open_h5", mock_open_h5, ns = "neuroarchive")
  on.exit({
    #message("Restoring original open_h5")
    unlockBinding("open_h5", asNamespace("neuroarchive"))
    assignInNamespace("open_h5", original_open_h5, ns = "neuroarchive")
    lockBinding("open_h5", asNamespace("neuroarchive"))
  }, add = TRUE) # Add = TRUE to ensure it runs even if other on.exit calls exist

  expect_error({
    # This call to read_lna should use our mocked open_h5
    read_lna(tmp, run_id = "run-nonexistent", lazy = TRUE)
  },
  class = "lna_error_run_id")
  
  #message(paste("After expect_error, captured_h5 class:", class(captured_h5), "is_valid:", if(inherits(captured_h5, "H5File")) captured_h5$is_valid else "NA"))
  
  # Check conditions on the captured H5 handle
  expect_true(inherits(captured_h5, "H5File"), info = "captured_h5 should be an H5File object.")
  # The lna_reader constructor should close the file if it errors out after opening
  expect_false(captured_h5$is_valid, info = "captured_h5 should be closed (invalid) after LNAECCReader initialization error.")
})


test_that("lna_reader close is idempotent", {
  tmp <- local_tempfile(fileext = ".h5")
  create_empty_lna(tmp)

  reader <- read_lna(tmp, lazy = TRUE)
  expect_true(reader$h5$is_valid)
  reader$close()
  expect_null(reader$h5)
  expect_silent(reader$close())
})


test_that("lna_reader close clears caches", {
  tmp <- local_tempfile(fileext = ".h5")
  create_empty_lna(tmp)

  reader <- read_lna(tmp, lazy = TRUE)
  reader$data()
  expect_false(is.null(reader$data_cache))
  expect_false(is.null(reader$cache_params))
  reader$close()
  expect_null(reader$data_cache)
  expect_null(reader$cache_params)
})


test_that("lna_reader data caches result and respects subset", {
  tmp <- local_tempfile(fileext = ".h5")
  create_empty_lna(tmp)

  reader <- read_lna(tmp, lazy = TRUE)
  h1 <- reader$data()
  expect_s3_class(h1, "DataHandle")
  h2 <- reader$data()
  expect_identical(h1, h2)

  reader$subset(roi_mask = 1)
  h3 <- reader$data()
  expect_false(identical(h1, h3))
  expect_identical(h3$subset$roi_mask, 1)

  h4 <- reader$data()
  expect_identical(h3, h4)

  reader$close()
})

test_that("read_lna lazy passes subset params", {
  tmp <- local_tempfile(fileext = ".h5")
  create_empty_lna(tmp)
  msk <- array(TRUE, dim = c(1,1,1))
  reader <- read_lna(tmp, lazy = TRUE, roi_mask = msk, time_idx = 2)
  expect_identical(reader$subset_params$roi_mask, msk)
  expect_identical(reader$subset_params$time_idx, 2L)
  reader$close()
})

test_that("lna_reader$subset validates parameters", {
  tmp <- local_tempfile(fileext = ".h5")
  create_empty_lna(tmp)
  reader <- read_lna(tmp, lazy = TRUE)

  expect_error(reader$subset(bad = 1), class = "lna_error_validation")
  expect_error(reader$subset(1), class = "lna_error_validation")

  reader$close()
})

test_that("lna_reader$data allow_plugins='none' errors on unknown transform", {
  tmp <- local_tempfile(fileext = ".h5")
  create_dummy_lna(tmp)
  reader <- read_lna(tmp, lazy = TRUE, allow_plugins = "none")
  expect_error(reader$data(), class = "lna_error_no_method")
  reader$close()
})

test_that("lna_reader$data allow_plugins='prompt' falls back when non-interactive", {
  tmp <- local_tempfile(fileext = ".h5")
  create_dummy_lna(tmp)
  reader <- read_lna(tmp, lazy = TRUE, allow_plugins = "prompt")
  
  testthat::with_mocked_bindings(
    is_interactive = function() FALSE,
    .package = "rlang",
    code = {
      expect_warning(reader$data(), "Missing invert_step")
    }
  )
  reader$close()
})

test_that("lna_reader$data allow_plugins='prompt' interactive respects choice", {
  tmp <- local_tempfile(fileext = ".h5")
  create_dummy_lna(tmp)
  
  # Test 'n' case
  reader_n <- read_lna(tmp, lazy = TRUE, allow_plugins = "prompt")
  testthat::with_mocked_bindings(
    is_interactive = function() TRUE,
    .package = "rlang",
    code = {
      testthat::with_mocked_bindings(
        readline = function(prompt = "") "n",
        .package = "base",
        code = {
          expect_error(reader_n$data(), class = "lna_error_no_method")
        }
      )
    }
  )
  reader_n$close()

  # Test 'y' case
  reader_y <- read_lna(tmp, lazy = TRUE, allow_plugins = "prompt")
  testthat::with_mocked_bindings(
    is_interactive = function() TRUE,
    .package = "rlang",
    code = {
      testthat::with_mocked_bindings(
        readline = function(prompt = "") "y",
        .package = "base",
        code = {
          expect_warning(reader_y$data(), "Missing invert_step")
        }
      )
    }
  )
  reader_y$close()
})

test_that("lna_reader$data errors when called after close", {
  tmp <- local_tempfile(fileext = ".h5")
  create_empty_lna(tmp)
  reader <- read_lna(tmp, lazy = TRUE)
  reader$close()
  expect_error(reader$data(), class = "lna_error_closed_reader")
})
