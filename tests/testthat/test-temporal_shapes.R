library(testthat)

# Helper: write small data for a given temporal kind and inspect dims
inspect_temporal_dims <- function(kind, n_time, n_vox, n_basis = NULL, extra = list()) {
  set.seed(123)
  X <- matrix(rnorm(n_time * n_vox), nrow = n_time, ncol = n_vox)
  arr <- array(t(X), dim = c(n_vox, 1, 1, n_time))
  tmp <- withr::local_tempfile(fileext = ".h5")
  p <- c(list(kind = kind), if (!is.null(n_basis)) list(n_basis = n_basis), extra)
  invisible(write_lna(arr, file = tmp, transforms = "temporal",
                      transform_params = list(temporal = p)))

  h5 <- neuroarchive:::open_h5(tmp, mode = "r")
  on.exit(neuroarchive:::close_h5_safely(h5), add = TRUE)
  tfg <- h5[["/transforms"]]
  nm <- names(tfg)[1]
  desc <- neuroarchive:::read_json_descriptor(tfg, nm)
  root <- h5[["/"]]
  roles <- vapply(desc$datasets, function(d) d$role, character(1))
  B_path <- desc$datasets[[which(roles == "temporal_basis")]]$path
  C_path <- desc$datasets[[which(roles == "temporal_coefficients")]]$path
  B <- neuroarchive:::h5_read(root, B_path)
  C <- neuroarchive:::h5_read(root, C_path)
  list(B = B, C = C)
}

test_that("temporal forward uses T x K basis and K x V coeffs (dct, poly)", {
  res_dct <- inspect_temporal_dims("dct", n_time = 10, n_vox = 4, n_basis = 6)
  expect_equal(dim(res_dct$B), c(10, 6))
  expect_equal(dim(res_dct$C), c(6, 4))

  res_poly <- inspect_temporal_dims("polynomial", n_time = 12, n_vox = 4, n_basis = 5)
  expect_equal(dim(res_poly$B), c(12, 5))
  expect_equal(dim(res_poly$C), c(5, 4))
})

test_that("temporal forward uses T x K basis and K x V coeffs (bspline)", {
  # Request too-large df to trigger internal clamping; verify shapes
  res_bs <- inspect_temporal_dims("bspline", n_time = 10, n_vox = 3, n_basis = 20, extra = list(order = 3))
  expect_equal(nrow(res_bs$B), 10)
  # basis is T x K, coeff is K x V
  expect_equal(ncol(res_bs$B), nrow(res_bs$C))
  expect_equal(dim(res_bs$C)[2], 3)
  expect_lte(ncol(res_bs$B), 10)
})

test_that("temporal forward uses T x K basis and K x V coeffs (dpss)", {
  # Force NW/n_basis clamp parameters; verify shapes
  res_dp <- inspect_temporal_dims("dpss", n_time = 16, n_vox = 2, n_basis = 10,
                                  extra = list(time_bandwidth_product = 4, n_tapers = 6))
  expect_equal(nrow(res_dp$B), 16)
  expect_equal(ncol(res_dp$B), nrow(res_dp$C))
  expect_equal(dim(res_dp$C)[2], 2)
  expect_lte(ncol(res_dp$B), 2 * (16/2))
})

test_that("temporal forward uses T x K basis and K x V coeffs (wavelet)", {
  res_wv <- inspect_temporal_dims("wavelet", n_time = 16, n_vox = 3, n_basis = 8,
                                  extra = list(wavelet = "db4"))
  expect_equal(dim(res_wv$B), c(16, 8))
  expect_equal(dim(res_wv$C), c(8, 3))
})
