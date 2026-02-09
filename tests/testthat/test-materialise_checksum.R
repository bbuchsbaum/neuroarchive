library(testthat)
library(hdf5r)
library(withr)
library(digest)

# Test checksum writing

test_that("materialise_plan writes verifiable sha256 checksum", {
  tmp <- local_tempfile(fileext = ".h5")

  p <- Plan$new()
  p$add_descriptor("00_dummy.json", list(type = "dummy"))
  p$add_payload("payload", matrix(1:4, nrow = 2))
  p$add_dataset_def("/scans/run-01/data", "data", "dummy", "run-01", 0L, "{}", "payload", "eager", dtype = NA_character_)

  h5 <- neuroarchive:::open_h5(tmp, mode = "w")
  neuroarchive:::materialise_plan(h5, p, checksum = "sha256")

  # Read stored hash
  h5_verify <- neuroarchive:::open_h5(tmp, mode = "r")
  root_verify <- h5_verify[["/"]]
  expect_true(neuroarchive:::h5_attr_exists(root_verify, "lna_checksum"))
  stored_hash <- neuroarchive:::h5_attr_read(root_verify, "lna_checksum")
  neuroarchive:::close_h5_safely(h5_verify)

  expect_true(nchar(stored_hash) == 64)
  expect_true(grepl("^[0-9a-f]{64}$", stored_hash))

  # Verify: raw-replace the stored hash with placeholder, hash file, compare.
  # This mirrors exactly what write_sha256_checksum does in reverse.
  tmp_copy <- local_tempfile(fileext = ".h5")
  file.copy(tmp, tmp_copy)
  placeholder <- paste(rep("0", 64), collapse = "")
  raw_bytes <- readBin(tmp_copy, "raw", file.info(tmp_copy)$size)
  stored_raw <- charToRaw(stored_hash)
  placeholder_raw <- charToRaw(placeholder)
  pos <- grepRaw(stored_raw, raw_bytes, fixed = TRUE)
  expect_equal(length(pos), 1L)  # hash should appear exactly once
  raw_bytes[pos:(pos + 63L)] <- placeholder_raw
  writeBin(raw_bytes, tmp_copy)
  expected_hash <- digest::digest(file = tmp_copy, algo = "sha256")

  expect_identical(stored_hash, expected_hash)
})
