library(testthat)

# Tests for internal wavelet basis generation and temporal_basis.wavelet

test_that(".wavelet_basis aliases produce identical results", {
  expect_equal(
    neuroarchive:::.wavelet_basis(8, wavelet = "db4"),
    neuroarchive:::.wavelet_basis(8, wavelet = "d4")
  )
  expect_equal(
    neuroarchive:::.wavelet_basis(8, wavelet = "db1"),
    neuroarchive:::.wavelet_basis(8, wavelet = "haar")
  )
})

test_that(".wavelet_basis requires power-of-two length", {
  expect_error(
    neuroarchive:::.wavelet_basis(10),
    class = "lna_error_validation"
  )
})

test_that("temporal_basis.wavelet respects n_basis", {
  full_basis <- neuroarchive::temporal_basis("wavelet", n_time = 8, n_basis = 8)
  sub_basis <- neuroarchive::temporal_basis("wavelet", n_time = 8, n_basis = 4)
  expect_equal(dim(full_basis), c(8, 8))
  expect_equal(dim(sub_basis), c(8, 4))
  expect_equal(sub_basis, full_basis[, 1:4, drop = FALSE])
})

