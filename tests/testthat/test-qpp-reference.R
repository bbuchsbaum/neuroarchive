test_that("QPP golden bytes and exact q roundtrips are platform-stable", {
  fixtures <- read_qpp_golden()
  expect_gt(length(fixtures), 0L)
  for (fixture in fixtures) {
    q <- golden_q_matrix(fixture)
    first <- .qpp_ref_encode_tile(q)
    second <- .qpp_ref_encode_tile(q)
    expect_identical(first$bytes, second$bytes, info = fixture$name)
    expect_identical(.qpp_ref_hex(first$bytes), fixture$hex, info = fixture$name)
    expect_identical(
      vapply(first$plan, `[[`, character(1), "predictor"),
      unlist(fixture$predictors, use.names = FALSE),
      info = fixture$name
    )
    expect_identical(
      vapply(first$plan, `[[`, integer(1), "width"),
      as.integer(unlist(fixture$widths, use.names = FALSE)),
      info = fixture$name
    )
    expect_identical(unname(.qpp_ref_decode_tile(.qpp_ref_unhex(fixture$hex))), unname(q),
                     info = fixture$name)
  }
})

test_that("P0 P1 and P2 residual oracles invert exactly", {
  q <- c(-12, -6, -1, 3, 8, 14, 21)
  for (predictor in c("P0", "P1", "P2")) {
    encoded <- .qpp_ref_residuals(q, predictor)
    decoded <- .qpp_ref_inverse_residuals(
      encoded$anchors,
      encoded$residuals,
      predictor
    )
    expect_identical(decoded, as.numeric(q))
  }
})

test_that("zig-zag widths and LSB packing cover exact boundaries", {
  signed <- c(-2^31, -65536, -128, -1, 0, 1, 127, 65535, 2^31 - 1)
  unsigned <- .qpp_ref_zigzag(signed)
  expect_identical(.qpp_ref_unzigzag(unsigned), signed)
  expect_identical(.qpp_ref_width(c(0)), 0L)
  for (width in c(1L, 2L, 7L, 8L, 15L, 16L, 31L, 32L, 33L)) {
    values <- c(0, 2^width - 1)
    packed <- .qpp_ref_pack(values, width)
    expect_identical(.qpp_ref_unpack(packed, 2L, width), values)
  }
  expect_identical(.qpp_ref_hex(.qpp_ref_pack(c(1, 2, 3), 2L)), "39")
})

test_that("predictor cost and ties are deterministic", {
  expect_identical(.qpp_ref_choose_predictor(rep(1, 16))$predictor, "P0")
  expect_identical(.qpp_ref_choose_predictor(rep(7, 12))$predictor, "P1")
  expect_identical(.qpp_ref_choose_predictor(-10:9)$predictor, "P2")
})

test_that("pointwise quantization uses ties-to-even and reports error bounds", {
  y <- matrix(c(-2.5, -1.5, -0.5, 0.5, 1.5, 2.5), ncol = 1)
  q <- .qpp_ref_quantize(y, 1)
  expect_identical(as.numeric(q), c(-2, -2, 0, 0, 2, 2))
  receipt <- .qpp_ref_reconstruct(q, 1, y)
  expect_lte(receipt$maximum_absolute_error, 0.5)

  synthetic <- .qpp_ref_synthetic_ar()
  delta <- c(0.01, 0.02, 0.04, 0.08)
  q2 <- .qpp_ref_quantize(synthetic, delta)
  receipt2 <- .qpp_ref_reconstruct(q2, delta, synthetic)
  expect_true(all(receipt2$parcelwise_maximum_absolute_error <=
                    receipt2$conservative_floating_bound))
})

test_that("synthetic calibration fixture is deterministic and autocorrelated", {
  first <- .qpp_ref_synthetic_ar()
  second <- .qpp_ref_synthetic_ar()
  expect_identical(first, second)
  lag_one <- vapply(seq_len(ncol(first)), function(k) {
    cor(first[-nrow(first), k], first[-1L, k])
  }, numeric(1))
  expect_true(all(lag_one > 0.25))
})

test_that("rate accounting separates payload padding and headers", {
  q <- cbind(rep(3, 20), -10:9, rep(c(-4, 4), 10))
  encoded <- .qpp_ref_encode_tile(q)
  metrics <- encoded$metrics
  expect_equal(metrics$total_bits, length(encoded$bytes) * 8)
  expect_equal(
    metrics$total_bits,
    metrics$header_and_anchor_bits + metrics$payload_storage_bits
  )
  expect_equal(metrics$effective_bits_per_parcel_tr,
               metrics$total_bits / length(q))
  expect_lte(metrics$payload_bits, metrics$payload_storage_bits)
})

test_that("chunk coverage accounts for bounded random-access reads", {
  index <- data.frame(
    time_start = c(1, 1, 65, 65),
    time_end = c(64, 64, 128, 128),
    parcel_start = c(1, 129, 1, 129),
    parcel_end = c(128, 256, 128, 256),
    byte_length = c(100, 110, 90, 95)
  )
  coverage <- .qpp_ref_chunk_coverage(index, 70:80, 120:140)
  expect_identical(coverage$tile_count, 2L)
  expect_identical(coverage$byte_count, 185)
  expect_identical(coverage$rows, c(3L, 4L))
})

test_that("reference decoder rejects malformed frames before reconstruction", {
  encoded <- .qpp_ref_encode_tile(matrix(-5:6, ncol = 1))$bytes
  expect_error(.qpp_ref_decode_tile(encoded[-length(encoded)]), "truncated")
  bad_magic <- encoded
  bad_magic[1L] <- as.raw(0)
  expect_error(.qpp_ref_decode_tile(bad_magic), "magic")
  expect_error(.qpp_ref_decode_tile(c(encoded, as.raw(0))), "trailing")
})
