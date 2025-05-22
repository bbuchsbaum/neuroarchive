library(testthat)
#library(neuroarchive)
library(withr)

# Simple aggregator plugin used for testing
.forward_step.myorg.aggregate_runs <- function(type, desc, handle) {
  lst <- handle$stash$input # Directly access the full input list from stash
  stopifnot(is.list(lst))
  mats <- lapply(lst, function(x) {
    if (is.matrix(x)) x else as.matrix(x)
  })
  aggregated <- do.call(rbind, mats)
  desc$version <- "1.0"
  desc$inputs <- "input" # What it conceptually consumes
  desc$outputs <- "aggregated_matrix" # What it produces
  # The plan descriptor should reflect this for provenance
  handle$plan$add_descriptor(handle$plan$get_next_filename(type), desc)
  # Update stash: remove original 'input' and add 'aggregated_matrix'
  handle$update_stash(keys = "input", new_values = list(aggregated_matrix = aggregated))
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
  run1_data <- matrix(rnorm(50), nrow = 10, ncol = 5)
  dim(run1_data) <- c(dim(run1_data), 1)
  run2_data <- matrix(rnorm(50), nrow = 10, ncol = 5)
  dim(run2_data) <- c(dim(run2_data), 1)
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
  expect_equal(dim(out), c(20,5,1))
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
  run1_data <- matrix(rnorm(50), nrow = 10, ncol = 5)
  dim(run1_data) <- c(dim(run1_data), 1)
  run2_data <- matrix(rnorm(50), nrow = 10, ncol = 5)
  dim(run2_data) <- c(dim(run2_data), 1)
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(list(`run-01` = run1_data, `run-02` = run2_data), file = tmp,
            transforms = c("myorg.aggregate_runs", "myorg.sparsepca", "delta", "temporal"),
            transform_params = list(myorg.sparsepca = list(n_components = 3),
                                  delta = list(order = 1L, axis = 2),
                                  temporal = list(n_basis = 5)))
  reader <- read_lna(tmp, lazy = TRUE, time_idx = 1:5)
  out <- reader$data()$stash$input
  expect_equal(dim(out), c(5,5,1))
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
  write_lna(list(`run-01` = run1_data, `run-02` = run2_data), file = tmp,
            transforms = c("myorg.aggregate_runs", "myorg.sparsepca"))
  h5 <- neuroarchive:::open_h5(tmp, mode = "r+")
  root <- h5[["/"]]
  hash_val <- digest::digest(file = tmp, algo = "sha256")
  neuroarchive:::h5_attr_write(root, "lna_checksum", hash_val)
  neuroarchive:::close_h5_safely(h5)
  expect_true(validate_lna(tmp))

  # corrupt file and check failure
  h5 <- neuroarchive:::open_h5(tmp, mode = "r+")
  ds <- h5$create_dataset("dummy", robj = 1L)
  neuroarchive:::close_h5_safely(h5)
  expect_error(validate_lna(tmp), class = "lna_error_validation")
})

