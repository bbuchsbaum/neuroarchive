library(testthat)
library(hdf5r)
library(jsonlite)
library(withr)

# Load functions if not running via devtools::test()
# source("../R/utils_json.R")

test_that("write_json_descriptor and read_json_descriptor round-trip works", {
  # Create a temporary HDF5 file path
  temp_h5_file <- withr::local_tempfile(fileext = ".h5")

  # Define test list
  test_list <- list(
    name = "Test Descriptor",
    version = 1.0,
    params = list(alpha = 0.05, n_comp = 10L),
    nested = list(a = TRUE, b = list(), c = c(1,2,3))
  )

  # Write the descriptor
  h5_file <- H5File$new(temp_h5_file, mode = "w")
  root_group_w <- h5_file[["/"]]
  write_json_descriptor(root_group_w, "desc_1", test_list)
  h5_file$close_all()

  # Read the descriptor back
  h5_file_read <- H5File$new(temp_h5_file, mode = "r")
  root_group_r <- h5_file_read[["/"]]
  read_list <- read_json_descriptor(root_group_r, "desc_1")

  # Check HDF5 type and space before closing
  dset <- root_group_r[["desc_1"]]
  dtype <- dset$get_type()
  dspace <- dset$get_space()

  expect_equal(as.character(dtype$get_class()), "H5T_STRING")
  

  # Close HDF5 objects
  dtype$close()
  dspace$close()
  dset$close()
  h5_file_read$close_all()

  # Compare original and read list using expect_equal for tolerance
  expect_equal(read_list, test_list)
})

test_that("write_json_descriptor is idempotent", {
  temp_h5_file <- withr::local_tempfile(fileext = ".h5")
  list1 <- list(a = 1, b = "first")
  list2 <- list(c = 2, d = "second")

  # Write initial list
  h5_file <- H5File$new(temp_h5_file, mode = "w")
  write_json_descriptor(h5_file[["/"]], "desc", list1)
  h5_file$close_all()

  # Write second list to the same name
  h5_file <- H5File$new(temp_h5_file, mode = "a") # Reopen in append mode
  root_group_a <- h5_file[["/"]]
  write_json_descriptor(root_group_a, "desc", list2)

  # Check that only one object named "desc" exists
  expect_equal(length(root_group_a$ls()$name), 1)
  expect_equal(root_group_a$ls()$name[1], "desc")
  h5_file$close_all()

  # Read back and check content
  h5_file_read <- H5File$new(temp_h5_file, mode = "r")
  read_list <- read_json_descriptor(h5_file_read[["/"]], "desc")
  h5_file_read$close_all()

  # Use expect_equal for the idempotency check due to potential type differences
  expect_equal(read_list, list2) # Should contain the second list
})

test_that("read_json_descriptor error handling works", {
  temp_h5_file <- withr::local_tempfile(fileext = ".h5")

  # Setup: Create file with one valid desc and one invalid JSON string
  h5_file <- H5File$new(temp_h5_file, mode = "w")
  root_group_w <- h5_file[["/"]]
  write_json_descriptor(root_group_w, "good_desc", list(a=1))

  # Manually create a dataset with bad JSON
  bad_json_string <- "{ bad json : "
  root_group_w[["bad_desc"]] <- bad_json_string
  
  # Manually create a dataset that isn't a string
  root_group_w[["numeric_desc"]] <- 123L
  
  h5_file$close_all()

  # --- Start tests ---
  h5_file_read <- H5File$new(temp_h5_file, mode = "r")
  root_group_r <- h5_file_read[["/"]]

  # 1. Read non-existent descriptor
  expect_error(
    read_json_descriptor(root_group_r, "missing_desc"),
    class = "lna_error_missing_path"
  )

  # 2. Read descriptor with invalid JSON
  expect_error(
    read_json_descriptor(root_group_r, "bad_desc"),
    class = "lna_error_json_parse"
  )

  # 3. Read descriptor that is not a string (or not scalar char)
  expect_error(
      read_json_descriptor(root_group_r, "numeric_desc"),
      class = "lna_error_invalid_descriptor"
  )

  # Close file
  h5_file_read$close_all()
}) 
