library(testthat)
library(hdf5r)
library(withr)
library(digest)

# Test checksum writing

test_that("materialise_plan writes sha256 checksum attribute that matches pre-attribute state", {
  tmp_main_file <- local_tempfile(fileext = ".h5")
  tmp_for_reference <- local_tempfile(fileext = ".h5")

  # Define how to create the plan content
  create_test_plan <- function() {
    p <- Plan$new()
    p$add_descriptor("00_dummy.json", list(type = "dummy"))
    p$add_payload("payload", matrix(1:4, nrow = 2))
    p$add_dataset_def("/scans/run-01/data", "data", "dummy", "run-01", 0L, "{}", "payload", "eager")
    p
  }

  # --- Create a reference file in the exact same way materialise_plan does ---
  # First, create with placeholder attribute
  planA <- create_test_plan()
  h5A <- neuroarchive:::open_h5(tmp_for_reference, mode = "w")
  root <- h5A[["/"]]
  
  # Write data (simplified version of materialise_plan without checksumming)
  h5A$create_group("transforms")
  h5A$create_group("basis")
  h5A$create_group("scans")
  
  neuroarchive:::h5_attr_write(root, "lna_spec", "LNA R v2.0")
  neuroarchive:::h5_attr_write(root, "creator", "lna R package v0.0.1")
  neuroarchive:::h5_attr_write(root, "required_transforms", character(0))
  
  h5_group <- h5A[["transforms"]]
  for (nm in names(planA$descriptors)) {
    neuroarchive:::write_json_descriptor(h5_group, nm, planA$descriptors[[nm]])
  }
  
  # Add payload
  for (i in seq_len(nrow(planA$datasets))) {
    row <- planA$datasets[i, ]
    key <- row$payload_key
    if (!nzchar(key)) next
    payload <- planA$payloads[[key]]
    neuroarchive:::h5_write_dataset(root, row$path, payload)
  }

  # Write placeholder checksum - THIS IS THE KEY STEP
  placeholder_checksum <- paste(rep("0", 64), collapse = "")
  neuroarchive:::h5_attr_write(root, "lna_checksum", placeholder_checksum)
  
  # Close file and calculate reference hash
  neuroarchive:::close_h5_safely(h5A)
  expected_hash_value <- digest::digest(file = tmp_for_reference, algo = "sha256")

  # --- Run materialise_plan with checksumming enabled on the main temp file ---
  plan2 <- create_test_plan() # Use a fresh plan object
  h5_actual_write <- neuroarchive:::open_h5(tmp_main_file, mode = "w")
  neuroarchive:::materialise_plan(h5_actual_write, plan2, checksum = "sha256")

  # --- Verify ---
  h5_verify <- neuroarchive:::open_h5(tmp_main_file, mode = "r")
  root_verify <- h5_verify[["/"]]
  expect_true(neuroarchive:::h5_attr_exists(root_verify, "lna_checksum"))
  actual_hash_in_attr <- neuroarchive:::h5_attr_read(root_verify, "lna_checksum")
  neuroarchive:::close_h5_safely(h5_verify)

  expect_identical(actual_hash_in_attr, expected_hash_value)
})
