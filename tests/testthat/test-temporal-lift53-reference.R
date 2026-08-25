lift53_reference_path <- test_inst_file(
  "prototypes", "temporal-lift53-reference.R"
)
source(lift53_reference_path, local = TRUE)

test_that("integer CDF 5/3 and Haar roundtrip odd short and multilevel units", {
  set.seed(20260825)
  for (wavelet in c("cdf53", "haar")) {
    for (length in 1:17) {
      series <- matrix(sample(-1000:1000, length * 3L, replace = TRUE),
                       nrow = length)
      for (levels in 0:5) {
        transformed <- lift_forward_unit(series, levels, wavelet)
        expect_identical(
          lift_inverse_unit(transformed), series * 1,
          info = paste(wavelet, "length", length, "levels", levels)
        )
      }
    }
  }
})

test_that("reference endpoint extension and rounding have golden coefficients", {
  cdf2 <- lift_forward_unit(matrix(c(4, 7), ncol = 1L), 1L, "cdf53")
  expect_identical(as.numeric(cdf2$coefficients), c(6, 3))
  cdf3 <- lift_forward_unit(matrix(c(4, 7, 10), ncol = 1L), 1L, "cdf53")
  expect_identical(as.numeric(cdf3$coefficients), c(4, 10, 0))
  haar3 <- lift_forward_unit(matrix(c(4, 7, 10), ncol = 1L), 1L, "haar")
  expect_identical(as.numeric(haar3$coefficients), c(5, 10, 3))
  expect_identical(
    cdf2$boundary_extension,
    "whole-sample symmetric endpoint duplication"
  )
  expect_match(cdf2$rounding, "floor")
})

test_that("coefficient ordering and band layout cover every sample exactly", {
  transformed <- lift_forward_unit(matrix(1:11, ncol = 1L), 3L, "cdf53")
  bands <- transformed$bands
  expect_identical(bands$kind, c("coarse", "detail", "detail", "detail"))
  expect_identical(bands$level, c(3L, 3L, 2L, 1L))
  expect_identical(bands$start[[1L]], 1L)
  expect_identical(tail(bands$end, 1L), 11L)
  expect_identical(sum(bands$length), 11L)
  expect_identical(
    unlist(Map(seq.int, bands$start, bands$end), use.names = FALSE),
    seq_len(11L)
  )
})

test_that("partitioned lifting never crosses declared units", {
  first <- matrix(c(1:9, 21:29), ncol = 2L)
  second <- matrix(c(101:107, 201:207), ncol = 2L)
  partitioned <- lift_forward_partitioned(
    list(run_01 = first, run_02 = second), 3L, "cdf53"
  )
  expect_false(partitioned$cross_unit_operations)
  expect_identical(partitioned$unit_ids, c("run_01", "run_02"))
  expect_identical(lift_inverse_partitioned(partitioned),
                   list(run_01 = first * 1, run_02 = second * 1))

  changed <- second
  changed[,] <- changed + 10000
  changed_partition <- lift_forward_partitioned(
    list(run_01 = first, run_02 = changed), 3L, "cdf53"
  )
  expect_identical(
    partitioned$units$run_01$coefficients,
    changed_partition$units$run_01$coefficients
  )
  global <- lift_forward_unit(rbind(first, second), 3L, "cdf53")
  expect_false(identical(
    global$coefficients[seq_len(nrow(first)), , drop = FALSE],
    partitioned$units$run_01$coefficients
  ))
})

test_that("lossy bands and progressive reconstruction are explicit", {
  series <- matrix(cumsum(rep(c(2, -1, 3, 0), 8L)), ncol = 2L)
  transformed <- lift_forward_unit(series, 3L, "cdf53")
  exact <- lift_quantize_bands(transformed, detail_steps = 1, dead_zone = 0)
  expect_identical(lift_inverse_quantized(exact), series * 1)

  lossy <- lift_quantize_bands(
    transformed, detail_steps = c(`1` = 4, `2` = 2, `3` = 1),
    dead_zone = 0.5, drop_finest = 1L
  )
  expect_identical(lossy$dropped_detail_levels, 1L)
  expect_identical(lossy$lossy_support, "complete-independence-unit")
  base <- lift_inverse_quantized(lossy, detail_levels = integer())
  medium <- lift_inverse_quantized(lossy, detail_levels = 3L)
  fuller <- lift_inverse_quantized(lossy, detail_levels = c(2L, 3L))
  expect_identical(dim(base), dim(series))
  expect_false(identical(base, medium))
  expect_false(identical(medium, fuller))
  expect_identical(
    lift_progressive_coefficients(lossy, c(1L, 2L, 3L)),
    lossy$reconstructed_coefficients
  )
})

test_that("integer reference rejects invalid domains and metadata", {
  expect_error(lift_forward_unit(matrix(0.5, 2L, 1L)), "exact integers")
  expect_error(
    lift_forward_unit(matrix(2147483648, 2L, 1L)), "signed-int32"
  )
  expect_error(lift_forward_unit(matrix(1:4, 2L), -1L), "non-negative")
  object <- lift_forward_unit(matrix(1:8, 4L), 2L)
  expect_error(lift_quantize_bands(object, detail_steps = 0), "positive")
  expect_error(lift_quantize_bands(object, drop_finest = 3L), "0..levels")
  quantized <- lift_quantize_bands(object)
  expect_error(lift_inverse_quantized(quantized, 3L), "subset")
})
