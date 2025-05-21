library(testthat)
library(hdf5r)
library(tibble)
library(withr)

# Source functions if not running via devtools::test()
# source("../R/discover.R")

# Helper function to create dummy datasets within an OPEN HDF5 group
create_dummy_transforms_in_group <- function(h5_group, transform_names) {
  if (length(transform_names) > 0) {
    # Use create_dataset + assign pattern for robustness
    dtype <- hdf5r::h5types$H5T_NATIVE_INT
    space <- H5S$new("scalar")
    # Cleanup happens when the group/file is closed
    # on.exit({ if(!is.null(space)) space$close(); if(!is.null(dtype)) dtype$close() }, add = TRUE)

    for (name in transform_names) {
      dset <- NULL
      tryCatch({
        # Create scalar integer dataset, preventing chunking
        dset <- h5_group$create_dataset(name, dtype = dtype, space = space, chunk_dim = NULL)
        # Write dummy data using slice assignment
        dset[] <- 0L
      }, finally = {
        # Ensure dataset is closed immediately after write
        if(!is.null(dset) && inherits(dset, "H5D")) dset$close()
      })
    }
    # Close space and dtype after the loop
    if (!is.null(space) && inherits(space, "H5S")) space$close()
    if (!is.null(dtype) && inherits(dtype, "H5T")) dtype$close()
  }
  invisible(NULL)
}

test_that("discover_transforms handles empty group", {
  temp_h5_file <- local_tempfile(fileext = ".h5")
  h5_file <- H5File$new(temp_h5_file, mode = "w")
  transforms_group <- h5_file$create_group("transforms")
  # Helper operates on the open group
  create_dummy_transforms_in_group(transforms_group, character(0))

  # Run discovery on the still-open group
  result <- discover_transforms(transforms_group)

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0)
  expect_equal(names(result), c("name", "type", "index"))
  expect_equal(result$name, character())
  expect_equal(result$type, character())
  expect_equal(result$index, integer())

  h5_file$close_all()
})

test_that("discover_transforms finds and orders correct sequence", {
  temp_h5_file <- local_tempfile(fileext = ".h5")
  h5_file <- H5File$new(temp_h5_file, mode = "w")
  transforms_group <- h5_file$create_group("transforms")
  transform_names <- c("01_scale.json", "00_mask.json", "02_pca.json")
  create_dummy_transforms_in_group(transforms_group, transform_names)

  result <- discover_transforms(transforms_group)

  expected <- tibble::tibble(
    name = c("00_mask.json", "01_scale.json", "02_pca.json"),
    type = c("mask", "scale", "pca"),
    index = c(0L, 1L, 2L)
  )

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 3)
  expect_identical(result, expected)

  h5_file$close_all()
})

test_that("discover_transforms errors on non-contiguous sequence", {
  temp_h5_file <- local_tempfile(fileext = ".h5")
  h5_file <- H5File$new(temp_h5_file, mode = "w")
  transforms_group <- h5_file$create_group("transforms")
  transform_names <- c("00_mask.json", "02_pca.json")
  create_dummy_transforms_in_group(transforms_group, transform_names)

  expect_error(
    discover_transforms(transforms_group),
    "Transform descriptor indices are not contiguous starting from 0. Found indices: 0, 2"
    # TODO: Update expected error class to lna_error_sequence when implemented
  )

  h5_file$close_all()
})

test_that("discover_transforms errors if sequence doesn't start at 0", {
  temp_h5_file <- local_tempfile(fileext = ".h5")
  h5_file <- H5File$new(temp_h5_file, mode = "w")
  transforms_group <- h5_file$create_group("transforms")
  transform_names <- c("01_mask.json", "02_pca.json")
  create_dummy_transforms_in_group(transforms_group, transform_names)

  expect_error(
    discover_transforms(transforms_group),
    "Transform descriptor indices are not contiguous starting from 0. Found indices: 1, 2"
  )

  h5_file$close_all()
})

test_that("discover_transforms errors on invalid names (no match)", {
  temp_h5_file <- local_tempfile(fileext = ".h5")
  h5_file <- H5File$new(temp_h5_file, mode = "w")
  transforms_group <- h5_file$create_group("transforms")
  transform_names <- c("00_mask.json", "invalid_name.txt", "01_scale.json")
  create_dummy_transforms_in_group(transforms_group, transform_names)

  expect_error(
    discover_transforms(transforms_group),
    "Invalid object name found in /transforms: invalid_name.txt. Expected format NN_type.json."
  )

  h5_file$close_all()
})

test_that("discover_transforms errors on non-numeric prefix", {
  temp_h5_file <- local_tempfile(fileext = ".h5")
  h5_file <- H5File$new(temp_h5_file, mode = "w")
  transforms_group <- h5_file$create_group("transforms")
  transform_names <- c("aa_mask.json", "01_scale.json")
  create_dummy_transforms_in_group(transforms_group, transform_names)

  expect_error(
    discover_transforms(transforms_group),
    "Invalid object name found in /transforms: aa_mask.json. Expected format NN_type.json."
  )

  h5_file$close_all()
})

test_that("discover_transforms errors if only invalid names found", {
  temp_h5_file <- local_tempfile(fileext = ".h5")
  h5_file <- H5File$new(temp_h5_file, mode = "w")
  transforms_group <- h5_file$create_group("transforms")
  transform_names <- c("aa.json", "bb.txt")
  create_dummy_transforms_in_group(transforms_group, transform_names)

  # Currently errors on the first invalid name found
  expect_error(
    discover_transforms(transforms_group),
    "Invalid object name found in /transforms: aa.json. Expected format NN_type.json."
  )

  h5_file$close_all()
}) 