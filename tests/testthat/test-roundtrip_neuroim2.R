library(testthat)
library(neuroarchive)
library(neuroim2)

# Synthetic roundtrip tests using neuroim2::simulate_fmri

make_simulated_vec <- function(dims = c(6L, 6L, 4L), n_time = 10L, seed = 42L) {
  sp <- neuroim2::NeuroSpace(dims, spacing = c(2, 2, 2))
  mask_arr <- array(TRUE, dim = dims)
  mask_vol <- neuroim2::NeuroVol(as.integer(mask_arr), sp)
  vec <- neuroim2::simulate_fmri(mask_vol, n_time = n_time, seed = seed)
  list(vec = methods::as(vec, "DenseNeuroVec"), mask = mask_arr)
}

as_time_by_vox <- function(x) {
  as_dense_mat(as.array(x))
}

test_that("neuroim2 identity roundtrip is exact", {
  sim <- make_simulated_vec()
  tmp <- local_tempfile(fileext = ".h5")

  expect_no_error(write_lna(sim$vec, file = tmp, transforms = character(0)))
  h <- read_lna(tmp)
  out <- as_dense_mat(h$stash$input)

  expect_equal(out, as_time_by_vox(sim$vec), tolerance = 1e-7)
})

test_that("neuroim2 quantization roundtrip is similar (lossy)", {
  sim <- make_simulated_vec(seed = 7L)
  tmp <- local_tempfile(fileext = ".h5")

  expect_no_error(
    write_lna(
      sim$vec, file = tmp, transforms = c("quant"),
      transform_params = list(quant = list(bits = 8L, method = "range", center = TRUE,
                                           scale_scope = "global", allow_clip = TRUE))
    )
  )
  h <- read_lna(tmp)
  rec <- as_dense_mat(h$stash$input)
  orig <- as_time_by_vox(sim$vec)

  # Similarity metrics
  mse <- mean((rec - orig)^2)
  var_orig <- mean(orig^2)
  rel_mse <- if (var_orig > 0) mse/var_orig else 0
  cor_val <- suppressWarnings(cor(as.numeric(rec), as.numeric(orig)))

  expect_lt(rel_mse, 0.05)       # <5% relative MSE on small synthetic data
  expect_gt(cor_val, 0.98)
})

test_that("neuroim2 delta (lossless) roundtrip is exact", {
  sim <- make_simulated_vec(seed = 11L)
  tmp <- local_tempfile(fileext = ".h5")

  expect_no_error(
    write_lna(
      sim$vec, file = tmp, transforms = c("delta"),
      transform_params = list(delta = list(order = 1L))
    )
  )
  h <- read_lna(tmp)
  out <- as_dense_mat(h$stash$input)
  expect_equal(out, as_time_by_vox(sim$vec), tolerance = 1e-7)
})

test_that("neuroim2 temporal DCT (full basis) is exact", {
  sim <- make_simulated_vec(seed = 13L, n_time = 12L)
  tmp <- local_tempfile(fileext = ".h5")

  expect_no_error(
    write_lna(
      sim$vec, file = tmp, transforms = c("temporal"),
      transform_params = list(temporal = list(kind = "dct", n_basis = 12L))
    )
  )
  h <- read_lna(tmp)
  out <- h$stash$input
  expect_equal(out, as_time_by_vox(sim$vec), tolerance = 1e-7)
})

test_that("neuroim2 temporal DCT (compressed) is similar", {
  sim <- make_simulated_vec(seed = 17L, n_time = 12L)
  tmp <- local_tempfile(fileext = ".h5")

  expect_no_error(
    write_lna(
      sim$vec, file = tmp, transforms = c("temporal"),
      transform_params = list(temporal = list(kind = "dct", n_basis = 6L))
    )
  )
  h <- read_lna(tmp)
  rec <- h$stash$input
  orig <- as_time_by_vox(sim$vec)

  mse <- mean((rec - orig)^2)
  var_orig <- mean(orig^2)
  rel_mse <- if (var_orig > 0) mse/var_orig else 0
  cor_val <- suppressWarnings(cor(as.numeric(rec), as.numeric(orig)))

  expect_lt(rel_mse, 0.45)
  expect_gt(cor_val, 0.75)
})

test_that("neuroim2 spat.haar_octwave roundtrip (no threshold) is exact", {
  sim <- make_simulated_vec(seed = 23L)
  tmp <- local_tempfile(fileext = ".h5")

  # Coerce to a base array to avoid S4 NeuroVec dispatch in internal helpers
  arr <- array(as.numeric(as.array(sim$vec)), dim = dim(as.array(sim$vec)))
  expect_no_error(
    write_lna(
      arr, file = tmp, mask = sim$mask,
      transforms = c("spat.haar_octwave"),
      transform_params = list(spat.haar_octwave = list(levels = 3L))
    )
  )
  h <- read_lna(tmp)
  out <- as_dense_mat(h$stash$input)
  expect_equal(out, as_time_by_vox(sim$vec), tolerance = 1e-6)
})

test_that("neuroim2 spat.haar_octwave + quant roundtrip is similar", {
  sim <- make_simulated_vec(seed = 29L)
  tmp <- local_tempfile(fileext = ".h5")

  arr <- array(as.numeric(as.array(sim$vec)), dim = dim(as.array(sim$vec)))
  expect_no_error(
    write_lna(
      arr, file = tmp, mask = sim$mask,
      transforms = c("spat.haar_octwave", "quant"),
      transform_params = list(
        spat.haar_octwave = list(levels = 3L),
        quant = list(bits = 8L, method = "range", center = TRUE,
                     scale_scope = "global", allow_clip = TRUE)
      )
    )
  )
  h <- read_lna(tmp)
  rec <- as_dense_mat(h$stash$input)
  orig <- as_time_by_vox(sim$vec)

  mse <- mean((rec - orig)^2)
  var_orig <- mean(orig^2)
  rel_mse <- if (var_orig > 0) mse/var_orig else 0
  cor_val <- suppressWarnings(cor(as.numeric(rec), as.numeric(orig)))

  expect_lt(rel_mse, 0.08)
  expect_gt(cor_val, 0.98)
})
