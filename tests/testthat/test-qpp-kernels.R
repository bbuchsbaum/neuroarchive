test_that("production QPP bytes match every independent golden fixture", {
  fixtures <- read_qpp_golden()
  for (fixture in fixtures) {
    q <- golden_q_matrix(fixture)
    bytes <- qpp_encode_tile(q)
    expect_identical(.qpp_ref_hex(bytes), fixture$hex, info = fixture$name)
    expect_identical(
      unname(qpp_decode_tile(.qpp_ref_unhex(fixture$hex))),
      unname(q),
      info = fixture$name
    )
  }
})

test_that("production kernels differentially match the Pure-R oracle", {
  set.seed(41023)
  for (iteration in seq_len(75L)) {
    n_time <- sample.int(41L, 1L)
    n_parcel <- sample.int(9L, 1L)
    values <- sample.int(2000001L, n_time * n_parcel, replace = TRUE) - 1000001L
    if (length(values) >= 2L && iteration %% 5L == 0L) {
      values[1:2] <- c(-2147483648, 2147483647)
    }
    q <- matrix(values, nrow = n_time, ncol = n_parcel)
    reference <- .qpp_ref_encode_tile(q)$bytes
    production <- qpp_encode_tile(q)
    expect_identical(production, reference, info = paste("iteration", iteration))
    expect_identical(unname(qpp_decode_tile(production)), unname(q * 1),
                     info = paste("iteration", iteration))
  }
})

test_that("independent tiles reset every predictor and anchor", {
  q <- cbind(cumsum(c(3, rep(1, 31))), rep(c(-8, 8), 16))
  first <- qpp_encode_tile(q[1:16, , drop = FALSE])
  second <- qpp_encode_tile(q[17:32, , drop = FALSE])
  expect_identical(qpp_decode_tile(first), q[1:16, , drop = FALSE] * 1)
  expect_identical(qpp_decode_tile(second), q[17:32, , drop = FALSE] * 1)
  expect_error(qpp_decode_tile(c(first, second)), "surplus")
})

test_that("QPP encoder rejects non-integers, nonfinite values, and limits", {
  expect_error(qpp_encode_tile(matrix(0.5, 1, 1)), "whole numbers")
  expect_error(qpp_encode_tile(matrix(NA_real_, 1, 1)), "whole numbers")
  expect_error(qpp_encode_tile(matrix(2147483648, 1, 1)), "signed-int32")
  expect_error(qpp_encode_tile(matrix(numeric(), 0, 1)), "dimensions")
  expect_error(qpp_encode_tile(matrix(1:4, 2, 2), max_values = 3),
               "max_values")
  expect_error(qpp_encode_tile(matrix(1, 1, 1), max_bytes = 8), "max_bytes")
})

test_that("QPP decoder fails on invalid modes widths buffers and padding", {
  bytes <- qpp_encode_tile(matrix(c(-1, -2, -4, -8, -16), ncol = 1))
  bad_mode <- bytes
  bad_mode[9L] <- as.raw(3L)
  expect_error(qpp_decode_tile(bad_mode), "header")

  bad_width <- bytes
  bad_width[10L] <- as.raw(35L)
  expect_error(qpp_decode_tile(bad_width), "header")

  bad_anchor <- qpp_encode_tile(matrix(1, 1, 1))
  bad_anchor[9L] <- as.raw(2L)
  bad_anchor[11L] <- as.raw(2L)
  expect_error(qpp_decode_tile(bad_anchor), "predictor lengths")

  for (n in 0:(length(bytes) - 1L)) {
    expect_error(qpp_decode_tile(bytes[seq_len(n)]), "truncated|magic")
  }
  expect_error(qpp_decode_tile(c(bytes, as.raw(0))), "surplus")

  padded <- qpp_encode_tile(matrix(1, 1, 1))
  padded[length(padded)] <- as.raw(bitwOr(as.integer(padded[length(padded)]), 128L))
  expect_error(qpp_decode_tile(padded), "padding")
})

test_that("QPP decoder checks lengths and resource ceilings before allocation", {
  bytes <- qpp_encode_tile(matrix(1:12, 4, 3))
  expect_error(qpp_decode_tile(bytes, max_time = 3L), "dimension")
  expect_error(qpp_decode_tile(bytes, max_parcels = 2L), "dimension")
  expect_error(qpp_decode_tile(bytes, max_values = 11), "max_values")
  expect_error(qpp_decode_tile(bytes, max_bytes = length(bytes) - 1), "max_bytes")

  malformed <- bytes
  # First record payload length is bytes 15..18 (little-endian uint32).
  malformed[15L] <- as.raw(as.integer(malformed[15L]) + 1L)
  expect_error(qpp_decode_tile(malformed), "payload length")
})

test_that("QPP decoder rejects inverse accumulation overflow", {
  overflow <- c(
    charToRaw("QPP1"), .qpp_ref_u16le(2), .qpp_ref_u16le(1),
    as.raw(c(1, 2, 1, 0)), .qpp_ref_u16le(1), .qpp_ref_u32le(1),
    .qpp_ref_i32le(2147483647), .qpp_ref_pack(2, 2)
  )
  expect_error(qpp_decode_tile(overflow), "overflows signed-int32")
})

test_that("QPP decoder supports valid 34-bit second differences", {
  residual <- 2^33 - 2
  frame <- c(
    charToRaw("QPP1"), .qpp_ref_u16le(3), .qpp_ref_u16le(1),
    as.raw(c(2, 34, 2, 0)), .qpp_ref_u16le(1), .qpp_ref_u32le(5),
    .qpp_ref_i32le(2147483647), .qpp_ref_i32le(-2147483648),
    .qpp_ref_pack(.qpp_ref_zigzag(residual), 34)
  )
  expect_identical(
    qpp_decode_tile(frame),
    matrix(c(2147483647, -2147483648, 2147483647), ncol = 1)
  )
})

test_that("QPP wrapper requires raw decoder input and valid limit scalars", {
  expect_error(qpp_decode_tile(1:3), "raw vector")
  bytes <- qpp_encode_tile(matrix(1, 1, 1))
  expect_error(qpp_decode_tile(bytes, max_time = 0L), "positive")
  expect_error(qpp_decode_tile(bytes, max_values = Inf), "exact integer")
})
