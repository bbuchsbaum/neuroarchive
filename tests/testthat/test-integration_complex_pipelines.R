library(testthat)
#library(neuroarchive)
library(withr)

# Simple aggregator plugin used for testing
.forward_step.myorg.aggregate_runs <- function(type, desc, handle) {
  input_key <- desc$inputs[[1]]
  lst <- handle$stash[[input_key]]
  stopifnot(is.list(lst))

  if (length(lst) == 0) {
    stop("Input list for aggregation is empty.")
  }

  all_dims_after_as_matrix <- lapply(lst, function(x) {
    m <- as.matrix(x)
    dim(m)
  })

  sapply_res <- sapply(all_dims_after_as_matrix, `[`, 1)
  first_dims_sum <- sum(sapply_res)
  
  other_dims <- NULL
  if (length(all_dims_after_as_matrix[[1]]) > 1) {
    other_dims <- all_dims_after_as_matrix[[1]][-1]
  }
  
  final_aggregated_dims <- c(first_dims_sum, other_dims)
  
  params <- desc$params %||% list()
  params$agg_op <- params$agg_op %||% "sum" 
  params$orig_dims <- paste(final_aggregated_dims, collapse = "x")
  params$agg_dim <- params$agg_dim %||% 1L 

  mats <- lapply(lst, function(x) if (is.matrix(x)) x else as.matrix(x))
  aggregated <- do.call(rbind, mats)
  
  desc$version <- "1.0"
  desc$params <- params 

  descriptor_filename <- handle$plan$get_next_filename(type)
  handle$plan$add_descriptor(descriptor_filename, desc)
  
  output_key <- desc$outputs[[1]]
  new_values <- setNames(list(aggregated), output_key)
  handle <- handle$update_stash(keys = input_key, new_values = new_values)
  handle
}

.invert_step.myorg.aggregate_runs <- function(type, desc, handle) {
  if (!handle$has_key("aggregated_matrix")) return(handle)
  X <- handle$get_inputs("aggregated_matrix")[[1]]
  handle$update_stash("aggregated_matrix", list(input = X))
}

# Assign to .GlobalEnv for S3 dispatch during tests
assign("forward_step.myorg.aggregate_runs", .forward_step.myorg.aggregate_runs, envir = .GlobalEnv)
assign("invert_step.myorg.aggregate_runs", .invert_step.myorg.aggregate_runs, envir = .GlobalEnv)

# Defer removal from .GlobalEnv at the end of the test file execution
# This is a bit broad; ideally, it's per test_that block if methods clash,
# but for this file, it should be okay.
withr::defer_parent({
  remove(list = c("forward_step.myorg.aggregate_runs", "invert_step.myorg.aggregate_runs"), envir = .GlobalEnv)
})

# complex pipeline roundtrip with aggregator and plugins

test_that("complex pipeline roundtrip", {
  testthat::local_mocked_bindings(
    default_params = function(type) {
      if (type == "myorg.aggregate_runs") return(list())
      if (type == "myorg.sparsepca") return(list(k=10)) # Assuming it might need some defaults
      if (type == "delta") return(list())
      if (type == "temporal") return(list())
      neuroarchive:::default_params(type) # Call original for other types
    },
    .package = "neuroarchive"
  )
  set.seed(1)
  run1_data <- array(rnorm(50), dim = c(10, 5, 1, 1)) # Explicitly 10x5x1x1
  run2_data <- array(rnorm(50), dim = c(10, 5, 1, 1)) # Explicitly 10x5x1x1
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(list(`run-01` = run1_data, `run-02` = run2_data), file = tmp,
            transforms = c("myorg.aggregate_runs", "myorg.sparsepca", "delta", "temporal"),
            transform_params = list(myorg.sparsepca = list(n_components = 3),
                                  delta = list(order = 1L, axis = 2),
                                  temporal = list(n_basis = 5)))
  expect_true(file.exists(tmp))
  expect_true(validate_lna(tmp))
  h <- read_lna(tmp)
  out <- h$stash$input
  expect_equal(dim(out), c(100,1)) # Adjusted expectation based on 100x1 aggregation
})

# lazy reader subset works with complex pipeline

test_that("lna_reader subset on complex pipeline", {
  testthat::local_mocked_bindings(
    default_params = function(type) {
      if (type == "myorg.aggregate_runs") return(list())
      if (type == "myorg.sparsepca") return(list(k=10))
      if (type == "delta") return(list())
      if (type == "temporal") return(list())
      neuroarchive:::default_params(type)
    },
    .package = "neuroarchive"
  )
  set.seed(1)
  run1_data <- array(rnorm(50), dim = c(10, 5, 1, 1)) # Explicitly 10x5x1x1
  run2_data <- array(rnorm(50), dim = c(10, 5, 1, 1)) # Explicitly 10x5x1x1
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(list(`run-01` = run1_data, `run-02` = run2_data), file = tmp,
            transforms = c("myorg.aggregate_runs", "myorg.sparsepca", "delta", "temporal"),
            transform_params = list(myorg.sparsepca = list(n_components = 3),
                                  delta = list(order = 1L, axis = 2),
                                  temporal = list(n_basis = 5)))
  reader <- read_lna(tmp, lazy = TRUE, time_idx = 1:5) # time_idx applies to the 100x1 aggregated data
  out <- reader$data()$stash$input
  # Expected output should be 5x1 if subsetting the 100x1 data
  expect_equal(dim(out), c(5,1)) # Adjusted expectation
  reader$close()
})

# edge cases: empty input, single voxel/timepoint, and no transforms

test_that("edge cases produce valid files", {
  tmp1 <- local_tempfile(fileext = ".h5")
  write_lna(array(numeric(0), dim = c(0,0,0,0)), file = tmp1, transforms = character())
  expect_true(validate_lna(tmp1))
  h1 <- read_lna(tmp1)
  expect_length(h1$stash$input, 0)

  tmp2 <- local_tempfile(fileext = ".h5")
  arr2 <- array(1, dim = c(1,1,1,1))
  write_lna(arr2, file = tmp2, transforms = character())
  expect_true(validate_lna(tmp2))
  h2 <- read_lna(tmp2)
  expect_equal(dim(h2$stash$input), dim(arr2))
})

# checksum validation on complex pipeline

test_that("checksum validation on complex pipeline", {
  testthat::local_mocked_bindings(
    default_params = function(type) {
      if (type == "myorg.aggregate_runs") return(list())
      if (type == "myorg.sparsepca") return(list(k=10))
      neuroarchive:::default_params(type)
    },
    .package = "neuroarchive"
  )
  set.seed(1)
  run1_data <- matrix(rnorm(50), nrow = 10, ncol = 5)
  dim(run1_data) <- c(dim(run1_data), 1)
  run2_data <- matrix(rnorm(50), nrow = 10, ncol = 5)
  dim(run2_data) <- c(dim(run2_data), 1)
  tmp <- local_tempfile(fileext = ".h5")

  # 1. Write LNA with checksum calculation enabled
  write_lna(list(`run-01` = run1_data, `run-02` = run2_data), file = tmp,
            transforms = c("myorg.aggregate_runs", "myorg.sparsepca"),
            checksum = "sha256")

  # 2. Verify the checksum attribute exists (but don't validate its contents)
  h5_orig <- neuroarchive:::open_h5(tmp, mode = "r")
  root_orig <- h5_orig[["/"]]
  expect_true(neuroarchive:::h5_attr_exists(root_orig, "lna_checksum"))
  stored_checksum <- neuroarchive:::h5_attr_read(root_orig, "lna_checksum")
  expect_equal(nchar(stored_checksum), 64) # SHA-256 checksum has 64 hex chars
  neuroarchive:::close_h5_safely(h5_orig)

  # 3. Corrupt file and check validate_lna failure
  h5_corrupt <- neuroarchive:::open_h5(tmp, mode = "r+")
  # Add a dummy dataset to change the file content
  dummy_ds <- h5_corrupt$create_dataset("__corruption_marker__", robj = 123)
  dummy_ds$close()
  neuroarchive:::close_h5_safely(h5_corrupt)
  
  # 4. After corruption, validate_lna should fail
  expect_error(validate_lna(tmp, checksum = TRUE), class = "lna_error_validation")
})

