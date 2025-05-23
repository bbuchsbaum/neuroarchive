library(testthat)
library(neuroarchive)

# Helper utilities for fake neuroim2 objects
FakeSpace <- function(dim, spacing_v) {
  structure(list(dim = dim, spacing = spacing_v, trans = diag(4), origin = c(0,0,0)),
            class = "FakeSpace")
}
space.LogicalNeuroVol <- function(x, ...) attr(x, "space")
spacing.FakeSpace <- function(x, ...) x$spacing
as.array.LogicalNeuroVol <- function(x, ...) x$arr


test_that("poisson_disk_sample_neuroim2 deterministic", {
  mask <- array(TRUE, dim = c(5,5,5))
  vol <- structure(list(arr = mask), class = "LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(5,5,5), c(1,1,1))

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.LogicalNeuroVol", space.LogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("as.array.LogicalNeuroVol", as.array.LogicalNeuroVol, envir=.GlobalEnv)
  withr::defer({
    rm(FakeSpace, space.LogicalNeuroVol, spacing.FakeSpace,
       as.array.LogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  c1 <- neuroarchive:::poisson_disk_sample_neuroim2(vol, radius_mm = 2, seed = 42)
  c2 <- neuroarchive:::poisson_disk_sample_neuroim2(vol, radius_mm = 2, seed = 42)
  expect_identical(c1, c2)
})


test_that("poisson_disk_sample_neuroim2 guard rail on tiny ROI", {
  mask <- array(TRUE, dim = c(2,2,2))
  vol <- structure(list(arr = mask), class = "LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(2,2,2), c(1,1,1))

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.LogicalNeuroVol", space.LogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("as.array.LogicalNeuroVol", as.array.LogicalNeuroVol, envir=.GlobalEnv)
  withr::defer({
    rm(FakeSpace, space.LogicalNeuroVol, spacing.FakeSpace,
       as.array.LogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  pts <- neuroarchive:::poisson_disk_sample_neuroim2(vol, radius_mm = 5, seed = 1)
  expect_equal(nrow(pts), 1)
})


test_that("poisson_disk_sample_neuroim2 handles disconnected components", {
  mask <- array(FALSE, dim = c(4,4,4))
  mask[1:2,1:2,1:2] <- TRUE
  mask[3:4,3:4,3:4] <- TRUE
  vol <- structure(list(arr = mask), class = "LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(4,4,4), c(1,1,1))

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.LogicalNeuroVol", space.LogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("as.array.LogicalNeuroVol", as.array.LogicalNeuroVol, envir=.GlobalEnv)
  withr::defer({
    rm(FakeSpace, space.LogicalNeuroVol, spacing.FakeSpace,
       as.array.LogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  pts <- neuroarchive:::poisson_disk_sample_neuroim2(vol, radius_mm = 1, seed = 99)
  expect_true(any(pts[,1] <= 2))
  expect_true(any(pts[,1] > 2))
})


test_that("generate_hrbf_atom gaussian normalisation", {
  coords <- matrix(rbind(c(0,0,0), c(1,0,0), c(0,1,0)), ncol = 3, byrow = TRUE)
  idx <- 1:3
  params <- list(kernel_type = "gaussian")
  res <- neuroarchive:::generate_hrbf_atom(coords, idx, c(0,0,0), sigma_mm = 1,
                                           current_level_j = 0, total_levels = 0,
                                           params = params)
  expect_equal(res$indices, idx)
  expect_equal(length(res$values), 3)
  expect_equal(sum(res$values^2), 1, tolerance = 1e-6)
})

test_that("generate_hrbf_atom wendland_c4 normalisation", {
  coords <- matrix(rbind(c(0,0,0), c(1,0,0), c(0,1,0)), ncol = 3, byrow = TRUE)
  idx <- 1:3
  params <- list(kernel_type = "gaussian",
                 kernel_type_fine_levels = "wendland_c4",
                 num_fine_levels_alt_kernel = 1L)
  res <- neuroarchive:::generate_hrbf_atom(coords, idx, c(0,0,0), sigma_mm = 2,
                                           current_level_j = 1, total_levels = 1,
                                           params = params)
  expect_equal(res$indices, idx)
  expect_equal(sum(res$values^2), 1, tolerance = 1e-6)
  expect_true(all(res$values >= 0))
})

test_that("compute_edge_map_neuroim2 self_mean", {
  mask <- array(TRUE, dim = c(3,3,3))
  vol <- structure(list(arr = mask), class = "LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(3,3,3), c(1,1,1))

  arr <- array(0, dim = c(3,3,3,2))
  arr[2,2,2,] <- 1

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.LogicalNeuroVol", space.LogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("as.array.LogicalNeuroVol", as.array.LogicalNeuroVol, envir=.GlobalEnv)
  withr::defer({
    rm(FakeSpace, space.LogicalNeuroVol, spacing.FakeSpace,
       as.array.LogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  h <- DataHandle$new(initial_stash = list(input_dense_mat = arr),
                      mask_info = list(mask = vol))
  edge <- neuroarchive:::compute_edge_map_neuroim2("self_mean", h,
                                                   list(edge_thresh_k = 0))
  expect_true(is.logical(edge))
  expect_equal(dim(edge), c(3,3,3))
  expect_true(any(edge))
})

test_that("compute_edge_map_neuroim2 structural_path", {
  mask <- array(TRUE, dim = c(2,2,2))
  vol <- structure(list(arr = mask), class = "LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(2,2,2), c(1,1,1))

  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)
  h5 <- H5File$new(tmp, mode = "w")
  h5[['grad']] <- array(2, dim = c(2,2,2))
  h5$close_all()
  h5 <- H5File$new(tmp, mode = "r")

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.LogicalNeuroVol", space.LogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("as.array.LogicalNeuroVol", as.array.LogicalNeuroVol, envir=.GlobalEnv)
  withr::defer({
    rm(FakeSpace, space.LogicalNeuroVol, spacing.FakeSpace,
       as.array.LogicalNeuroVol, envir=.GlobalEnv)
    h5$close_all()
  }, envir = parent.frame())

  h <- DataHandle$new(h5 = h5, mask_info = list(mask = vol))
  edge <- neuroarchive:::compute_edge_map_neuroim2(
    "structural_path", h,
    list(structural_path = "grad", edge_thresh_k = 0.5)
  )
  expect_true(all(edge))
})


test_that("poisson_disk_sample_neuroim2 edge adaptation favors edges", {
  mask <- array(TRUE, dim = c(5,5,1))
  vol <- structure(list(arr = mask), class = "LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(5,5,1), c(1,1,1))

  edge_map <- array(FALSE, dim = c(5,5,1))
  edge_map[1:2,,] <- TRUE

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.LogicalNeuroVol", space.LogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("as.array.LogicalNeuroVol", as.array.LogicalNeuroVol, envir=.GlobalEnv)
  withr::defer({
    rm(FakeSpace, space.LogicalNeuroVol, spacing.FakeSpace,
       as.array.LogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  pts_base <- neuroarchive:::poisson_disk_sample_neuroim2(vol, radius_mm = 2,
                                                           seed = 7)
  pts_edge <- neuroarchive:::poisson_disk_sample_neuroim2(vol, radius_mm = 2,
                                                           seed = 7,
                                                           edge_binary_map = edge_map,
                                                           density_factor = 2)
  frac_base <- mean(pts_base[,1] <= 2)
  frac_edge <- mean(pts_edge[,1] <= 2)
  expect_gt(frac_edge, frac_base)
})


test_that("combined upgrades improve reconstruction MSE", {
  mask <- array(TRUE, dim = c(3,3,3))
  vol <- structure(list(arr = mask), class = "LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(3,3,3), c(1,1,1))

  X <- matrix(rnorm(27 * 2), nrow = 2)

  params_base <- list(sigma0 = 6, levels = 0, radius_factor = 2.5,
                      kernel_type = "gaussian", seed = 1)
  params_adv <- list(sigma0 = 6, levels = 0, radius_factor = 2.5,
                     kernel_type = "gaussian",
                     kernel_type_fine_levels = "wendland_c4",
                     num_fine_levels_alt_kernel = 1L,
                     num_extra_fine_levels = 1L,
                     seed = 1)

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.LogicalNeuroVol", space.LogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("as.array.LogicalNeuroVol", as.array.LogicalNeuroVol, envir=.GlobalEnv)
  withr::defer({
    rm(FakeSpace, space.LogicalNeuroVol, spacing.FakeSpace,
       as.array.LogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  B_base <- hrbf_generate_basis(params_base, vol)
  coeff_base <- hrbf_project_matrix(X, vol, params_base)
  recon_base <- hrbf_reconstruct_matrix(coeff_base, vol, params_base)
  mse_base <- mean((X - recon_base)^2)

  B_adv <- hrbf_generate_basis(params_adv, vol)
  coeff_adv <- hrbf_project_matrix(X, vol, params_adv)
  recon_adv <- hrbf_reconstruct_matrix(coeff_adv, vol, params_adv)
  mse_adv <- mean((X - recon_adv)^2)

  expect_lt(mse_adv, 0.9 * mse_base)
})

