library(testthat)

test_that("typo option name produces warning with suggestion", {
  expect_warning(
    lna_options(write.compresion_level = 4L),
    "Unknown option.*write\\.compresion_level.*Did you mean"
  )
})

test_that("completely unknown option name produces warning", {
  expect_warning(
    lna_options(totally_fake_option = 42),
    "Unknown option.*totally_fake_option"
  )
})

test_that("wrong type for integer option produces error", {
  expect_error(
    lna_options(write.compression_level = "banana"),
    "expects.*integer|numeric"
  )
})

test_that("wrong type for numeric option produces error", {
  expect_error(
    lna_options(write.chunk_target_mib = "not_a_number"),
    "expects.*numeric"
  )
})

test_that("wrong type for logical option produces error", {
  expect_error(
    lna_options(read.strict_mask_hash_validation = "yes"),
    "expects.*logical"
  )
})

test_that("correct option and type succeeds silently", {
  old <- lna_options("write.compression_level")
  on.exit(lna_options(write.compression_level = old[[1]]), add = TRUE)

  expect_no_warning(lna_options(write.compression_level = 6L))
  expect_equal(lna_options("write.compression_level")[[1]], 6L)
})

test_that("numeric value accepted for integer option", {
  old <- lna_options("write.compression_level")
  on.exit(lna_options(write.compression_level = old[[1]]), add = TRUE)

  # Numeric is compatible with integer
  expect_no_error(lna_options(write.compression_level = 6))
})

test_that(".strict=FALSE suppresses unknown name warning", {
  expect_no_warning(
    lna_options(totally_bogus = 1, .strict = FALSE)
  )
})

test_that("list-typed option accepts list", {
  old <- lna_options("quant")
  on.exit(lna_options(quant = old[[1]]), add = TRUE)

  expect_no_error(lna_options(quant = list(bits = 12)))
})
