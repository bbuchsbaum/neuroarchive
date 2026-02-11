library(testthat)

test_that("write_lna returns lna_write_result object", {
  arr <- array(rnorm(256), dim = c(4, 4, 4, 4))
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  result <- write_lna(arr, tmp, transforms = "quant")

  expect_s3_class(result, "lna_write_result")
  expect_true(is.numeric(result$input_bytes))
  expect_true(is.numeric(result$output_bytes))
  expect_true(is.numeric(result$ratio))
  expect_equal(result$transforms_applied, "quant")
  expect_equal(result$file, tmp)
})

test_that("write result ratio > 1 for quant transform", {
  arr <- array(rnorm(4000), dim = c(10, 10, 10, 4))
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  result <- write_lna(arr, tmp, transforms = "quant")
  expect_gt(result$ratio, 1)
})

test_that("write result output_bytes matches file size", {
  arr <- array(rnorm(256), dim = c(4, 4, 4, 4))
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  result <- write_lna(arr, tmp, transforms = "quant")
  expect_equal(result$output_bytes, file.size(tmp))
})

test_that("print.lna_write_result produces output", {
  arr <- array(rnorm(256), dim = c(4, 4, 4, 4))
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  result <- write_lna(arr, tmp, transforms = "quant")
  output <- capture.output(print(result))
  expect_true(length(output) > 0)
  expect_true(grepl("LNA written:", output[1]))
  expect_true(grepl("x\\)", output[1]))
})

test_that("write_lna without assignment still works", {
  arr <- array(rnorm(256), dim = c(4, 4, 4, 4))
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  # Should not error when return value is ignored
  expect_no_error(write_lna(arr, tmp, transforms = "quant"))
  expect_true(file.exists(tmp))
})

test_that("write_lna with multiple transforms reports all in summary", {
  arr <- array(rnorm(4000), dim = c(10, 10, 10, 4))
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  result <- write_lna(arr, tmp, transforms = c("quant", "basis"))
  expect_equal(result$transforms_applied, c("quant", "basis"))

  output <- capture.output(print(result))
  expect_true(grepl("\\[quant\\+basis\\]", output[1]))
})

test_that("write_lna with no transforms reports [none]", {
  arr <- array(rnorm(256), dim = c(4, 4, 4, 4))
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  result <- write_lna(arr, tmp, transforms = character())

  output <- capture.output(print(result))
  expect_true(grepl("\\[none\\]", output[1]))
})

test_that("write_lna checksum field is preserved in result", {
  arr <- array(rnorm(256), dim = c(4, 4, 4, 4))
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  result <- write_lna(arr, tmp, transforms = "quant", checksum = "sha256")
  expect_equal(result$checksum, "sha256")

  output <- capture.output(print(result))
  expect_true(grepl("sha256:ok", output[1]))
})

test_that("write_lna with checksum='none' shows 'no checksum'", {
  arr <- array(rnorm(256), dim = c(4, 4, 4, 4))
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)

  result <- write_lna(arr, tmp, transforms = "quant", checksum = "none")
  expect_equal(result$checksum, "none")

  output <- capture.output(print(result))
  expect_true(grepl("no checksum", output[1]))
})
