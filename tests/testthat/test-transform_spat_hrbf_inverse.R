library(neuroarchive)
library(withr)

FakeSpace <- function(dim, spacing_v, origin_v=c(0,0,0)) {
  structure(list(dim=dim, spacing=spacing_v, origin=origin_v), class="FakeSpace")
}
space.FakeLogicalNeuroVol <- function(x, ...) attr(x, "space")
spacing.FakeSpace <- function(x, ...) x$spacing
origin.FakeSpace <- function(x, ...) x$origin
as.array.FakeLogicalNeuroVol <- function(x, ...) x$arr


# Basic reconstruction ---------------------------------------------------------

test_that("invert_step.spat.hrbf reconstructs dense data", {
  mask <- array(TRUE, dim=c(2,2,2))
  vol <- structure(list(arr=mask), class="LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(2,2,2), c(1,1,1))

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.FakeLogicalNeuroVol", space.FakeLogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("origin.FakeSpace", origin.FakeSpace, envir=.GlobalEnv)
  assign("as.array.FakeLogicalNeuroVol", as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  defer({
    rm(FakeSpace, space.FakeLogicalNeuroVol, spacing.FakeSpace,
       origin.FakeSpace, as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  X <- matrix(1:8, nrow=2)
  plan <- Plan$new()
  h <- DataHandle$new(initial_stash=list(input_dense_mat=X), plan=plan,
                      mask_info=list(mask=vol, active_voxels=8))
  desc <- list(type="spat.hrbf",
               params=list(sigma0=6, levels=0, radius_factor=2.5,
                            kernel_type="gaussian", seed=1))
  h2 <- neuroarchive:::forward_step.spat.hrbf("spat.hrbf", desc, h)
  fdesc <- h2$plan$descriptors[[1]]
  coeff <- h2$stash$coefficients_hrbf

  tmp <- local_tempfile(fileext=".h5")
  h5 <- H5File$new(tmp, mode="w")
  on.exit(h5$close_all(), add=TRUE)
  h_inv <- DataHandle$new(initial_stash=list(coefficients_hrbf=coeff), h5=h5,
                          mask_info=list(mask=vol, active_voxels=8))
  out <- neuroarchive:::invert_step.spat.hrbf("spat.hrbf", fdesc, h_inv)
  expect_true(out$has_key("input_dense_mat"))
  expect_false(out$has_key("coefficients_hrbf"))
  expect_equal(out$stash$input_dense_mat, X)
})

# Mask hash mismatch ----------------------------------------------------------

test_that("invert_step.spat.hrbf mask hash mismatch warns/errors", {
  mask <- array(TRUE, dim=c(2,2,2))
  vol <- structure(list(arr=mask), class="LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(2,2,2), c(1,1,1))

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.FakeLogicalNeuroVol", space.FakeLogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("origin.FakeSpace", origin.FakeSpace, envir=.GlobalEnv)
  assign("as.array.FakeLogicalNeuroVol", as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  defer({
    rm(FakeSpace, space.FakeLogicalNeuroVol, spacing.FakeSpace,
       origin.FakeSpace, as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  X <- matrix(1:8, nrow=2)
  plan <- Plan$new()
  h <- DataHandle$new(initial_stash=list(input_dense_mat=X), plan=plan,
                      mask_info=list(mask=vol, active_voxels=8))
  desc <- list(type="spat.hrbf", params=list(sigma0=6, levels=0, radius_factor=2.5,
                                             kernel_type="gaussian", seed=1))
  h2 <- neuroarchive:::forward_step.spat.hrbf("spat.hrbf", desc, h)
  fdesc <- h2$plan$descriptors[[1]]
  coeff <- h2$stash$coefficients_hrbf

  mask2 <- mask
  mask2[1,1,1] <- FALSE
  vol2 <- structure(list(arr=mask2), class="LogicalNeuroVol")
  attr(vol2, "space") <- FakeSpace(c(2,2,2), c(1,1,1))

  tmp <- local_tempfile(fileext=".h5")
  h5 <- H5File$new(tmp, mode="w")
  on.exit(h5$close_all(), add=TRUE)

  h_inv <- DataHandle$new(initial_stash=list(coefficients_hrbf=coeff), h5=h5,
                          mask_info=list(mask=vol2, active_voxels=7))
  old <- lna_options("read.strict_mask_hash_validation")$read.strict_mask_hash_validation
  defer(lna_options(read.strict_mask_hash_validation = old), envir=parent.frame())
  lna_options(read.strict_mask_hash_validation = FALSE)
  expect_warning(neuroarchive:::invert_step.spat.hrbf("spat.hrbf", fdesc, h_inv))
  lna_options(read.strict_mask_hash_validation = TRUE)
  expect_error(neuroarchive:::invert_step.spat.hrbf("spat.hrbf", fdesc, h_inv),
               class="lna_error_validation")
})

# Determinism ------------------------------------------------------------

test_that("invert_step.spat.hrbf deterministic basis regeneration", {
  mask <- array(TRUE, dim = c(1,1,2))
  vol <- structure(list(arr = mask), class = "LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(1,1,2), c(1,1,1))

  assign("FakeSpace", FakeSpace, envir = .GlobalEnv)
  assign("space.FakeLogicalNeuroVol", space.FakeLogicalNeuroVol, envir = .GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir = .GlobalEnv)
  assign("origin.FakeSpace", origin.FakeSpace, envir = .GlobalEnv)
  assign("as.array.FakeLogicalNeuroVol", as.array.FakeLogicalNeuroVol, envir = .GlobalEnv)
  defer({
    rm(FakeSpace, space.FakeLogicalNeuroVol, spacing.FakeSpace,
       origin.FakeSpace, as.array.FakeLogicalNeuroVol, envir = .GlobalEnv)
  }, envir = parent.frame())

  X <- matrix(1:4, nrow = 2)
  plan <- Plan$new()
  h <- DataHandle$new(initial_stash = list(input_dense_mat = X), plan = plan,
                      mask_info = list(mask = vol, active_voxels = 2))
  desc <- list(type = "spat.hrbf",
               params = list(sigma0 = 6, levels = 0, radius_factor = 2.5,
                            kernel_type = "gaussian", seed = 1))
  h2 <- neuroarchive:::forward_step.spat.hrbf("spat.hrbf", desc, h)
  fdesc <- h2$plan$descriptors[[1]]
  coeff <- h2$stash$coefficients_hrbf

  tmp <- local_tempfile(fileext = ".h5")
  h5 <- H5File$new(tmp, mode = "w")
  on.exit(h5$close_all(), add = TRUE)

  h_inv1 <- DataHandle$new(initial_stash = list(coefficients_hrbf = coeff), h5 = h5,
                           mask_info = list(mask = vol, active_voxels = 2))
  out1 <- neuroarchive:::invert_step.spat.hrbf("spat.hrbf", fdesc, h_inv1)

  h_inv2 <- DataHandle$new(initial_stash = list(coefficients_hrbf = coeff), h5 = h5,
                           mask_info = list(mask = vol, active_voxels = 2))
  out2 <- neuroarchive:::invert_step.spat.hrbf("spat.hrbf", fdesc, h_inv2)

  expect_equal(out1$stash$input_dense_mat, out2$stash$input_dense_mat)
})

# Stored dense matrix ----------------------------------------------------

test_that("invert_step.spat.hrbf uses stored dense matrix when available", {
  mask <- array(TRUE, dim = c(1,1,2))
  vol <- structure(list(arr = mask), class = "LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(1,1,2), c(1,1,1))

  assign("FakeSpace", FakeSpace, envir = .GlobalEnv)
  assign("space.FakeLogicalNeuroVol", space.FakeLogicalNeuroVol, envir = .GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir = .GlobalEnv)
  assign("origin.FakeSpace", origin.FakeSpace, envir = .GlobalEnv)
  assign("as.array.FakeLogicalNeuroVol", as.array.FakeLogicalNeuroVol, envir = .GlobalEnv)
  defer({
    rm(FakeSpace, space.FakeLogicalNeuroVol, spacing.FakeSpace,
       origin.FakeSpace, as.array.FakeLogicalNeuroVol, envir = .GlobalEnv)
  }, envir = parent.frame())

  X <- matrix(1:4, nrow = 2)
  plan <- Plan$new()
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- H5File$new(tmp, mode = "w")

  h <- DataHandle$new(initial_stash = list(input_dense_mat = X), plan = plan,
                      h5 = h5,
                      mask_info = list(mask = vol, active_voxels = 2))
  desc <- list(type = "spat.hrbf",
               params = list(sigma0 = 6, levels = 0, radius_factor = 2.5,
                            kernel_type = "gaussian", seed = 1,
                            store_dense_matrix = TRUE))

  h2 <- neuroarchive:::forward_step.spat.hrbf("spat.hrbf", desc, h)

  mat_path <- "/basis/hrbf/analytic/matrix"
  neuroarchive:::h5_write_dataset(h5[["/"]], mat_path, h2$stash$hrbf_basis)

  fdesc <- h2$plan$descriptors[[1]]
  fdesc$params$seed <- 9999
  coeff <- h2$stash$coefficients_hrbf
  h_inv <- DataHandle$new(initial_stash = list(coefficients_hrbf = coeff), h5 = h5,
                          mask_info = list(mask = vol, active_voxels = 2))
  out <- neuroarchive:::invert_step.spat.hrbf("spat.hrbf", fdesc, h_inv)

  expect_equal(out$stash$input_dense_mat, X)
})

