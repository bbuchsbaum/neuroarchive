library(testthat)
library(hdf5r)
library(withr)

test_that("validate_lna succeeds on valid file", {
  tmp <- local_tempfile(fileext = ".h5")
  create_valid_lna(tmp)
  expect_true(validate_lna(tmp))
})


test_that("validate_lna detects spec mismatch", {
  skip_if(requireNamespace("covr", quietly = TRUE) && covr::in_covr(), "validate_lna mutation tests unreliable under covr")
  tmp <- local_tempfile(fileext = ".h5")
  create_valid_lna(tmp)
  mutate_lna_file(tmp, function(h5) {
    root <- h5[["/"]]
    h5_attr_write(root, "lna_spec", "wrong")
  })
  expect_error(validate_lna(tmp), class = "lna_error_validation")
})


test_that("validate_lna detects checksum mismatch", {
  skip_if(requireNamespace("covr", quietly = TRUE) && covr::in_covr(), "validate_lna mutation tests unreliable under covr")
  tmp <- local_tempfile(fileext = ".h5")
  create_valid_lna(tmp)
  mutate_lna_file(tmp, function(h5) {
    root <- h5[["/"]]
    h5_attr_write(root, "lna_checksum", "bogus")
  })
  expect_error(validate_lna(tmp), class = "lna_error_validation")
  expect_true(validate_lna(tmp, checksum = FALSE))
})

test_that("validate_lna validates descriptor schemas", {
  tmp <- local_tempfile(fileext = ".h5")
  create_schema_lna(tmp)
  expect_true(validate_lna(tmp))
})

test_that("validate_lna fails on invalid descriptor", {
  skip_if(requireNamespace("covr", quietly = TRUE) && covr::in_covr(), "validate_lna mutation tests unreliable under covr")
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- neuroarchive:::open_h5(tmp, mode = "w")
  plan <- Plan$new()
  plan$add_descriptor("00_basis.json", list(type = "basis",
                                            params = list(method = "bogus")))
  materialise_plan(h5, plan, checksum = "none")
  expect_error(validate_lna(tmp), class = "lna_error_validation")
})

test_that("validate_lna strict=FALSE collects multiple issues", {
  skip_if(requireNamespace("covr", quietly = TRUE) && covr::in_covr(), "validate_lna mutation tests unreliable under covr")
  tmp <- local_tempfile(fileext = ".h5")
  create_valid_lna(tmp)
  mutate_lna_file(tmp, function(h5) {
    root <- h5[["/"]]
    h5_attr_write(root, "lna_spec", "wrong")
    h5_attr_write(root, "lna_checksum", "bogus")
  })

  res <- validate_lna(tmp, strict = FALSE)
  expect_s3_class(res, "lna_validation_report")
  expect_false(res$ok)
  expect_s3_class(res$issues, "data.frame")
  expect_equal(nrow(res$issues), 2)
  expect_named(res$issues, c("level", "code", "path", "layer", "message"))
  expect_equal(as.character(res), res$issues$message)
})

test_that("validate_lna strict=TRUE errors on first issue", {
  skip_if(requireNamespace("covr", quietly = TRUE) && covr::in_covr(), "validate_lna mutation tests unreliable under covr")
  tmp <- local_tempfile(fileext = ".h5")
  create_valid_lna(tmp)
  mutate_lna_file(tmp, function(h5) {
    root <- h5[["/"]]
    h5_attr_write(root, "lna_spec", "wrong")
    h5_attr_write(root, "lna_checksum", "bogus")
  })

  expect_error(validate_lna(tmp, strict = TRUE), class = "lna_error_validation")
})

test_that("validate_lna detects missing required groups", {
  skip_if(requireNamespace("covr", quietly = TRUE) && covr::in_covr(), "validate_lna mutation tests unreliable under covr")
  tmp <- local_tempfile(fileext = ".h5")
  create_valid_lna(tmp)
  mutate_lna_file(tmp, function(h5) h5$link_delete("basis"))
  expect_error(validate_lna(tmp), class = "lna_error_validation")
})

test_that("validate_lna detects missing dataset referenced by descriptor", {
  skip_if(requireNamespace("covr", quietly = TRUE) && covr::in_covr(), "validate_lna mutation tests unreliable under covr")
  tmp <- local_tempfile(fileext = ".h5")
  create_valid_lna(tmp)
  mutate_lna_file(tmp, function(h5) h5$link_delete("scans/run-01/data"))
  expect_error(validate_lna(tmp), class = "lna_error_validation")
})

test_that("validate_lna detects dimension mismatch hints", {
  # HDF5 file mutation tests are unreliable under covr instrumentation:
  # the checksum and validation paths behave differently when the package
  # is rebuilt with source instrumentation.
  skip_if(requireNamespace("covr", quietly = TRUE) && covr::in_covr(), "validate_lna mutation tests unreliable under covr")
  tmp <- local_tempfile(fileext = ".h5")
  create_valid_lna(tmp)
  mutate_lna_file(tmp, function(h5) {
    tf <- h5[["transforms"]]
    desc <- read_json_descriptor(tf, "00_dummy.json")
    desc$datasets[[1]]$dims <- c(1L, 1L)
    write_json_descriptor(tf, "00_dummy.json", desc)
  })
  expect_error(validate_lna(tmp), class = "lna_error_validation")
})

test_that("validate_lna errors when dataset cannot be read", {
  skip_if(requireNamespace("covr", quietly = TRUE) && covr::in_covr(), "validate_lna mutation tests unreliable under covr")
  tmp <- local_tempfile(fileext = ".h5")
  create_valid_lna(tmp)
  mutate_lna_file(tmp, function(h5) {
    h5$link_delete("scans/run-01/data")
    h5$create_group("scans/run-01/data")
  })

  expect_error(validate_lna(tmp), class = "lna_error_validation")
})
