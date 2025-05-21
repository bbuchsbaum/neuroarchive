library(testthat)
library(hdf5r)
library(withr)

# Source functions if not running via devtools::test()
# source("../R/utils_hdf5.R")

test_that("HDF5 attribute helpers work on H5Group", {
  temp_h5_file <- local_tempfile(fileext = ".h5")
  h5_file <- H5File$new(temp_h5_file, mode = "w")
  root_group <- h5_file[["/"]]

  # --- Test Data ---
  attr_int    <- 123L
  attr_dbl    <- 456.789
  attr_logi   <- TRUE
  attr_char   <- "Test String"
  attr_int_v  <- c(1L, 2L, 3L)
  attr_dbl_v  <- c(1.1, 2.2, 3.3)
  attr_char_v <- c("a", "b", "c")

  # --- Initial State Checks ---
  expect_false(h5_attr_exists(root_group, "attr_int"))
  expect_error(h5_attr_read(root_group, "attr_int"), "Attribute .* not found")
  expect_null(h5_attr_delete(root_group, "attr_int")) # Deleting non-existent is no-op

  # --- Write Attributes ---
  expect_null(h5_attr_write(root_group, "attr_int",    attr_int))
  expect_null(h5_attr_write(root_group, "attr_dbl",    attr_dbl))
  expect_null(h5_attr_write(root_group, "attr_logi",   attr_logi))
  expect_null(h5_attr_write(root_group, "attr_char",   attr_char))
  expect_null(h5_attr_write(root_group, "attr_int_v",  attr_int_v))
  expect_null(h5_attr_write(root_group, "attr_dbl_v",  attr_dbl_v))
  expect_null(h5_attr_write(root_group, "attr_char_v", attr_char_v))

  # --- Existence Checks After Write ---
  expect_true(h5_attr_exists(root_group, "attr_int"))
  expect_true(h5_attr_exists(root_group, "attr_char_v"))

  # --- Read and Verify Attributes ---
  expect_identical(h5_attr_read(root_group, "attr_int"),    attr_int)
  expect_identical(h5_attr_read(root_group, "attr_dbl"),    attr_dbl)
  expect_identical(h5_attr_read(root_group, "attr_logi"),   attr_logi)
  expect_identical(h5_attr_read(root_group, "attr_char"),   attr_char)
  expect_identical(h5_attr_read(root_group, "attr_int_v"),  attr_int_v)
  expect_identical(h5_attr_read(root_group, "attr_dbl_v"),  attr_dbl_v)
  expect_identical(h5_attr_read(root_group, "attr_char_v"), attr_char_v)

  # --- Test Overwrite ---
  new_char <- "Overwritten"
  expect_null(h5_attr_write(root_group, "attr_char", new_char))
  expect_true(h5_attr_exists(root_group, "attr_char"))
  expect_identical(h5_attr_read(root_group, "attr_char"), new_char)

  # --- Test Delete ---
  expect_true(h5_attr_exists(root_group, "attr_int"))
  expect_null(h5_attr_delete(root_group, "attr_int"))
  expect_false(h5_attr_exists(root_group, "attr_int"))
  expect_error(h5_attr_read(root_group, "attr_int"), "Attribute .* not found")

  # Delete remaining attributes
  expect_null(h5_attr_delete(root_group, "attr_dbl"))
  expect_null(h5_attr_delete(root_group, "attr_logi"))
  expect_null(h5_attr_delete(root_group, "attr_char")) # Already overwritten & deleted above? No, re-wrote
  expect_null(h5_attr_delete(root_group, "attr_int_v"))
  expect_null(h5_attr_delete(root_group, "attr_dbl_v"))
  expect_null(h5_attr_delete(root_group, "attr_char_v"))

  # Final existence check
  expect_false(h5_attr_exists(root_group, "attr_dbl_v"))

  # --- Close File ---
  h5_file$close_all()
})

test_that("HDF5 attribute helpers work on H5D (Dataset)", {
  temp_h5_file <- local_tempfile(fileext = ".h5")
  h5_file <- H5File$new(temp_h5_file, mode = "w")
  test_data <- matrix(1:12, nrow = 3, ncol = 4)
  dset <- h5_file$create_dataset("test_dset", test_data)

  # --- Test Data ---
  attr_ds <- "Attribute on dataset"

  # --- Initial State Checks ---
  expect_false(h5_attr_exists(dset, "ds_attr"))

  # --- Write, Exist, Read ---
  expect_null(h5_attr_write(dset, "ds_attr", attr_ds))
  expect_true(h5_attr_exists(dset, "ds_attr"))
  expect_identical(h5_attr_read(dset, "ds_attr"), attr_ds)

  # --- Overwrite ---
  new_ds_attr <- "New DS Attribute"
  expect_null(h5_attr_write(dset, "ds_attr", new_ds_attr))
  expect_identical(h5_attr_read(dset, "ds_attr"), new_ds_attr)

  # --- Delete ---
  expect_true(h5_attr_exists(dset, "ds_attr"))
  expect_null(h5_attr_delete(dset, "ds_attr"))
  expect_false(h5_attr_exists(dset, "ds_attr"))

  # --- Close Dataset and File ---
  dset$close()
  h5_file$close_all()
})

test_that("HDF5 attribute helpers handle edge cases and errors", {
  temp_h5_file <- local_tempfile(fileext = ".h5")
  h5_file <- H5File$new(temp_h5_file, mode = "w")
  root_group <- h5_file[["/"]]

  # Invalid object type
  expect_error(h5_attr_write(h5_file, "bad", 1), "must be an H5Group or H5D object")
  expect_error(h5_attr_read(h5_file, "bad"), "must be an H5Group or H5D object")
  expect_error(h5_attr_exists(h5_file, "bad"), "must be an H5Group or H5D object")
  expect_error(h5_attr_delete(h5_file, "bad"), "must be an H5Group or H5D object")

  # Invalid name type
  expect_error(h5_attr_write(root_group, 123, 1), "is.character\\(name\\) is not TRUE")
  expect_error(h5_attr_write(root_group, c("a","b"), 1), "length\\(name\\) == 1 is not TRUE")

  # Read non-existent
  expect_error(h5_attr_read(root_group, "does_not_exist"), "Attribute .* not found")

  # Delete non-existent (should be silent)
  expect_null(h5_attr_delete(root_group, "does_not_exist"))

  h5_file$close_all()
}) 