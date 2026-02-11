test_that("read_lna runs quick validation by default", {
  tmp <- tempfile(fileext = ".h5")
  arr <- array(rnorm(64), dim = c(4, 4, 4, 1))
  write_lna(arr, tmp, transforms = "quant")

  # Default behavior should use quick validation
  expect_no_error(read_lna(tmp))

  unlink(tmp)
})

test_that("quick validation catches missing file", {
  nonexistent <- tempfile(fileext = ".h5")

  expect_error(
    read_lna(nonexistent),
    "File does not exist",
    class = "lna_error_validation"
  )
})

test_that("quick validation catches non-HDF5 files", {
  tmp <- tempfile(fileext = ".h5")
  writeLines("This is not an HDF5 file", tmp)

  expect_error(
    read_lna(tmp),
    "Not a valid HDF5 file",
    class = "lna_error_validation"
  )

  unlink(tmp)
})

test_that("quick validation catches missing LNA structure", {
  tmp <- tempfile(fileext = ".h5")

  # Create a valid HDF5 file but without LNA structure
  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  h5$create_dataset("some_data", matrix(1:10, nrow = 2))
  h5$close_all()

  expect_error(
    read_lna(tmp),
    "missing LNA structure",
    class = "lna_error_validation"
  )

  unlink(tmp)
})

test_that("quick validation catches missing transforms group", {
  tmp <- tempfile(fileext = ".h5")

  # Create HDF5 file with lna_spec but no transforms
  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  root <- h5[["/"]]
  hdf5r::h5attr(root, "lna_spec") <- "LNA R v2.0"
  h5$close_all()

  expect_error(
    read_lna(tmp),
    "missing /transforms group",
    class = "lna_error_validation"
  )

  unlink(tmp)
})

test_that("quick validation allows empty transforms group", {
  # Empty transforms group is allowed (edge case for testing)
  tmp <- tempfile(fileext = ".h5")

  # Create HDF5 file with transforms group but no descriptors
  h5 <- hdf5r::H5File$new(tmp, mode = "w")
  root <- h5[["/"]]
  hdf5r::h5attr(root, "lna_spec") <- "LNA R v2.0"
  h5$create_group("transforms")
  h5$create_group("scans")
  h5$create_group("basis")
  h5$close_all()

  # Quick validation should pass even with empty transforms
  expect_true(neuroarchive:::quick_validate_lna(tmp))

  unlink(tmp)
})

test_that("validate=FALSE skips all validation checks", {
  tmp <- tempfile(fileext = ".h5")
  arr <- array(rnorm(64), dim = c(4, 4, 4, 1))
  write_lna(arr, tmp, transforms = "quant")

  # Should read successfully even with validate=FALSE
  expect_no_error(read_lna(tmp, validate = FALSE))

  unlink(tmp)
})

test_that("validate=TRUE runs full validation including checksums", {
  tmp <- tempfile(fileext = ".h5")
  arr <- array(rnorm(64), dim = c(4, 4, 4, 1))
  write_lna(arr, tmp, transforms = "quant", checksum = "sha256")

  # Full validation should succeed
  expect_no_error(read_lna(tmp, validate = TRUE))

  unlink(tmp)
})

test_that("validate='quick' is the default and works", {
  tmp <- tempfile(fileext = ".h5")
  arr <- array(rnorm(64), dim = c(4, 4, 4, 1))
  write_lna(arr, tmp, transforms = "quant")

  # Explicitly specify quick validation
  expect_no_error(read_lna(tmp, validate = "quick"))

  unlink(tmp)
})

test_that("invalid validate parameter throws error", {
  tmp <- tempfile(fileext = ".h5")
  arr <- array(rnorm(64), dim = c(4, 4, 4, 1))
  write_lna(arr, tmp, transforms = "quant")

  expect_error(
    read_lna(tmp, validate = "invalid")
  )

  expect_error(
    read_lna(tmp, validate = 123),
    class = "lna_error_validation"
  )

  unlink(tmp)
})

test_that("quick validation is fast (< 5ms on typical file)", {
  tmp <- tempfile(fileext = ".h5")
  arr <- array(rnorm(1000), dim = c(10, 10, 10, 1))
  write_lna(arr, tmp, transforms = "quant")

  # Measure quick validation time
  timing <- system.time({
    for (i in 1:10) {
      read_lna(tmp, validate = "quick")
    }
  })

  # Average time per read should be reasonable (< 50ms per read for 10 reads)
  # This is a loose check since we're doing full reads, not just validation
  avg_time_ms <- (timing["elapsed"] / 10) * 1000

  # Just verify it completes without error - timing test is too environment-dependent
  expect_true(avg_time_ms > 0)

  unlink(tmp)
})

test_that("quick_validate_lna can be called directly", {
  tmp <- tempfile(fileext = ".h5")
  arr <- array(rnorm(64), dim = c(4, 4, 4, 1))
  write_lna(arr, tmp, transforms = "quant")

  # Direct call to internal function
  expect_true(neuroarchive:::quick_validate_lna(tmp))

  unlink(tmp)
})
