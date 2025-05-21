library(testthat)
library(neuroarchive)
library(withr)

# Simple aggregator plugin used for testing
forward_step.myorg.aggregate_runs <- function(type, desc, handle) {
  lst <- handle$get_inputs("input")[[1]]
  stopifnot(is.list(lst))
  mats <- lapply(lst, function(x) {
    if (is.matrix(x)) x else as.matrix(x)
  })
  aggregated <- do.call(rbind, mats)
  desc$version <- "1.0"
  desc$inputs <- "input"
  desc$outputs <- "aggregated_matrix"
  handle$plan$add_descriptor(handle$plan$get_next_filename(type), desc)
  handle$update_stash("input", list(aggregated_matrix = aggregated))
}

invert_step.myorg.aggregate_runs <- function(type, desc, handle) {
  if (!handle$has_key("aggregated_matrix")) return(handle)
  X <- handle$get_inputs("aggregated_matrix")[[1]]
  handle$update_stash("aggregated_matrix", list(input = X))
}

# complex pipeline roundtrip with aggregator and plugins

test_that("complex pipeline roundtrip", {
  set.seed(1)
  run1 <- matrix(rnorm(50), nrow = 10, ncol = 5)
  run2 <- matrix(rnorm(50), nrow = 10, ncol = 5)
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(list(run1 = run1, run2 = run2), file = tmp,
            transforms = c("myorg.aggregate_runs", "myorg.sparsepca", "delta", "temporal"),
            transform_params = list(
              `myorg.sparsepca` = list(k = 3),
              delta = list(axis = 2),
              temporal = list(n_basis = 5)
            ))
  expect_true(file.exists(tmp))
  expect_true(validate_lna(tmp))
  h <- read_lna(tmp)
  out <- h$stash$input
  expect_equal(dim(out), c(20,5))
})

# lazy reader subset works with complex pipeline

test_that("lna_reader subset on complex pipeline", {
  set.seed(1)
  run1 <- matrix(rnorm(50), nrow = 10, ncol = 5)
  run2 <- matrix(rnorm(50), nrow = 10, ncol = 5)
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(list(run1 = run1, run2 = run2), file = tmp,
            transforms = c("myorg.aggregate_runs", "myorg.sparsepca", "delta", "temporal"),
            transform_params = list(
              `myorg.sparsepca` = list(k = 3),
              delta = list(axis = 2),
              temporal = list(n_basis = 5)
            ))
  reader <- read_lna(tmp, lazy = TRUE, time_idx = 1:5)
  out <- reader$data()$stash$input
  expect_equal(dim(out), c(5,5))
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
  set.seed(1)
  run1 <- matrix(rnorm(50), nrow = 10, ncol = 5)
  run2 <- matrix(rnorm(50), nrow = 10, ncol = 5)
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(list(run1 = run1, run2 = run2), file = tmp,
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

