library(testthat)

# Tests for .encode_rle and .decode_rle helpers

test_that("encode and decode rle round-trip", {
  vec_num <- c(1L, 1L, 2L, 2L, 2L, 3L)
  encoded <- neuroarchive:::.encode_rle(vec_num)
  expect_equal(encoded[, "lengths"], c(2L, 3L, 1L))
  expect_equal(encoded[, "values"], c(1L, 2L, 3L))
  decoded <- neuroarchive:::.decode_rle(encoded)
  expect_equal(decoded, vec_num)

  vec_empty <- numeric()
  encoded_empty <- neuroarchive:::.encode_rle(vec_empty)
  expect_equal(nrow(encoded_empty), 0)
  decoded_empty <- neuroarchive:::.decode_rle(encoded_empty)
  expect_equal(decoded_empty, vec_empty)
})


test_that("decode_rle validates input and expected length", {
  expect_error(
    neuroarchive:::.decode_rle(c(1, 1, 2)),
    class = "lna_error_runtime"
  )

  good_enc <- neuroarchive:::.encode_rle(c(5, 5, 6))
  expect_error(
    neuroarchive:::.decode_rle(good_enc, expected_length = 4,
                               location = "test"),
    class = "lna_error_runtime"
  )
})
