library(testthat)
library(hdf5r)
library(withr)
library(digest)

# Test checksum writing

test_that("materialise_plan writes sha256 checksum attribute that matches pre-attribute state", {
  tmp_main_file <- local_tempfile(fileext = ".h5")
  tmp_for_expected_hash1 <- local_tempfile(fileext = ".h5")
  tmp_for_expected_hash2 <- local_tempfile(fileext = ".h5") # For determinism check

  # Define how to create the plan content
  create_test_plan <- function() {
    p <- Plan$new()
    p$add_descriptor("00_dummy.json", list(type = "dummy"))
    p$add_payload("payload", matrix(1:4, nrow = 2))
    p$add_dataset_def("/scans/run-01/data", "data", "dummy", "run-01", 0L, "{}", "payload", "eager")
    p
  }

  # --- Sanity Check: Determinism of materialise_plan(checksum="none") ---
  planA <- create_test_plan()
  h5A <- neuroarchive:::open_h5(tmp_for_expected_hash1, mode = "w")
  neuroarchive:::materialise_plan(h5A, planA, checksum = "none")
  neuroarchive:::close_h5_safely(h5A)
  hashA <- digest::digest(file = tmp_for_expected_hash1, algo = "sha256")

  planB <- create_test_plan()
  h5B <- neuroarchive:::open_h5(tmp_for_expected_hash2, mode = "w")
  neuroarchive:::materialise_plan(h5B, planB, checksum = "none")
  neuroarchive:::close_h5_safely(h5B)
  hashB <- digest::digest(file = tmp_for_expected_hash2, algo = "sha256")

  expect_identical(hashA, hashB,
                   info = "Checksums from two identical materialise_plan('none') calls should match.")

  # --- Calculate Expected Checksum Matching materialise_plan("sha256") logic ---
  placeholder <- paste(rep("0", 64), collapse = "")
  h5_expected <- neuroarchive:::open_h5(tmp_for_expected_hash1, mode = "r+")
  root_expected <- h5_expected[["/"]]
  neuroarchive:::h5_attr_write(root_expected, "lna_checksum", placeholder)
  neuroarchive:::close_h5_safely(h5_expected)
  expected_hash_value <- digest::digest(file = tmp_for_expected_hash1, algo = "sha256")

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
