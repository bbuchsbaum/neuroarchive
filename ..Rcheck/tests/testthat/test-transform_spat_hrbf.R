library(neuroarchive)

# Fake neuroim2 helpers for mask
FakeSpace <- function(dim, spacing_v, origin_v=c(0,0,0)) {
  structure(list(dim=dim, spacing=spacing_v, origin=origin_v), class="FakeSpace")
}
space.FakeLogicalNeuroVol <- function(x, ...) attr(x,"space")
spacing.FakeSpace <- function(x, ...) x$spacing
origin.FakeSpace <- function(x, ...) x$origin
as.array.FakeLogicalNeuroVol <- function(x, ...) x$arr


test_that("forward_step.spat.hrbf generates centres and hash", {
  mask <- array(TRUE, dim=c(2,2,2))
  vol <- structure(list(arr=mask), class="LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(2,2,2), c(1,1,1))

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.FakeLogicalNeuroVol", space.FakeLogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("origin.FakeSpace", origin.FakeSpace, envir=.GlobalEnv)
  assign("as.array.FakeLogicalNeuroVol", as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  withr::defer({
    rm(FakeSpace, space.FakeLogicalNeuroVol, spacing.FakeSpace,
       origin.FakeSpace, as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  plan <- Plan$new()
  h <- DataHandle$new(initial_stash=list(input=matrix(0, nrow=1, ncol=8)),
                      plan=plan, mask_info=list(mask=vol, active_voxels=8))
  desc <- list(type="spat.hrbf",
               params=list(sigma0=6, levels=1, radius_factor=2.5,
                            kernel_type="gaussian", seed=42))
  h2 <- neuroarchive:::forward_step.spat.hrbf("spat.hrbf", desc, h)

  expect_true(h2$has_key("hrbf_centres"))
  expect_true(h2$has_key("hrbf_sigmas"))
  dname <- names(h2$plan$descriptors)[1]
  stored_desc <- h2$plan$descriptors[[dname]]
  expect_true(startsWith(stored_desc$params$mask_hash, "sha256:"))
  expect_equal(stored_desc$params$k_actual, nrow(h2$stash$hrbf_centres))
})

test_that("forward_step.spat.hrbf stores basis matrix when requested", {
  mask <- array(TRUE, dim=c(2,2,2))
  vol <- structure(list(arr=mask), class="LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(2,2,2), c(1,1,1))

  plan <- Plan$new()
  h <- DataHandle$new(initial_stash=list(input=matrix(0, nrow=1, ncol=8)),
                      plan=plan, mask_info=list(mask=vol, active_voxels=8))
  desc <- list(type="spat.hrbf",
               params=list(sigma0=6, levels=0, radius_factor=2.5,
                            kernel_type="gaussian", seed=42,
                            store_dense_matrix=TRUE))
  h2 <- neuroarchive:::forward_step.spat.hrbf("spat.hrbf", desc, h)

  mat_path <- "/basis/hrbf/analytic/matrix"
  expect_true(mat_path %in% names(h2$plan$payloads))
  B <- h2$plan$payloads[[mat_path]]
  expect_true(inherits(B, "Matrix"))
  expect_equal(dim(B)[2], length(as.array(vol)))
})

test_that("forward_step.spat.hrbf computes coefficients", {
  mask <- array(TRUE, dim = c(1,1,2))
  vol <- structure(list(arr = mask), class = "LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(1,1,2), c(1,1,1))

  X <- matrix(1:4, nrow = 2)

  plan <- Plan$new()
  h <- DataHandle$new(initial_stash = list(input = X),
                      plan = plan, mask_info = list(mask = vol, active_voxels = 2))
  desc <- list(type = "spat.hrbf",
               params = list(sigma0 = 6, levels = 0, radius_factor = 2.5,
                            kernel_type = "gaussian", seed = 1))

  h2 <- neuroarchive:::forward_step.spat.hrbf("spat.hrbf", desc, h)

  coef_path <- "/scans/run-01/embedding/coefficients_hrbf"
  expect_true(coef_path %in% names(h2$plan$payloads))
  expect_true(h2$has_key("coefficients_hrbf"))
  C <- h2$stash$coefficients_hrbf
  expect_equal(nrow(C), nrow(X))
})

test_that("num_extra_fine_levels increases k_actual", {
  mask <- array(TRUE, dim = c(3,3,3))
  vol <- structure(list(arr = mask), class = "LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(3,3,3), c(1,1,1))

  plan1 <- Plan$new()
  h1 <- DataHandle$new(initial_stash = list(input = matrix(0, nrow = 1, ncol = 27)),
                       plan = plan1, mask_info = list(mask = vol, active_voxels = 27))
  desc1 <- list(type = "spat.hrbf",
                params = list(sigma0 = 6, levels = 0, radius_factor = 2.5,
                               kernel_type = "gaussian", seed = 1,
                               num_extra_fine_levels = 0))
  out1 <- neuroarchive:::forward_step.spat.hrbf("spat.hrbf", desc1, h1)
  k1 <- out1$plan$descriptors[[1]]$params$k_actual

  plan2 <- Plan$new()
  h2 <- DataHandle$new(initial_stash = list(input = matrix(0, nrow = 1, ncol = 27)),
                       plan = plan2, mask_info = list(mask = vol, active_voxels = 27))
  desc2 <- list(type = "spat.hrbf",
                params = list(sigma0 = 6, levels = 0, radius_factor = 2.5,
                               kernel_type = "gaussian", seed = 1,
                               num_extra_fine_levels = 1))
  out2 <- neuroarchive:::forward_step.spat.hrbf("spat.hrbf", desc2, h2)
  k2 <- out2$plan$descriptors[[1]]$params$k_actual

  expect_gt(k2, k1)
})

test_that("forward_step.spat.hrbf warns when anisotropic atoms requested", {
  mask <- array(TRUE, dim = c(1,1,1))
  vol <- structure(list(arr = mask), class = "LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(1,1,1), c(1,1,1))

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.FakeLogicalNeuroVol", space.FakeLogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("origin.FakeSpace", origin.FakeSpace, envir=.GlobalEnv)
  assign("as.array.FakeLogicalNeuroVol", as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  withr::defer({
    rm(FakeSpace, space.FakeLogicalNeuroVol, spacing.FakeSpace,
       origin.FakeSpace, as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  tmp <- local_tempfile(fileext = ".h5")
  h5 <- H5File$new(tmp, mode = "w")
  h5$create_group("tensor")

  plan <- Plan$new()
  h <- DataHandle$new(initial_stash = list(input = matrix(0, nrow = 1, ncol = 1)),
                      plan = plan, h5 = h5,
                      mask_info = list(mask = vol, active_voxels = 1))
  desc <- list(type = "spat.hrbf",
               params = list(sigma0 = 6, levels = 0, radius_factor = 2.5,
                              kernel_type = "gaussian", seed = 1,
                              use_anisotropic_atoms = TRUE,
                              anisotropy_source_path = "tensor"))

  expect_warning(
    neuroarchive:::forward_step.spat.hrbf("spat.hrbf", desc, h),
    regexp = "Anisotropic atoms"
  )
  h5$close_all()
})

test_that("forward_step.spat.hrbf warns for derivative Gaussian atoms", {
  mask <- array(TRUE, dim = c(1,1,1))
  vol <- structure(list(arr = mask), class = "LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(1,1,1), c(1,1,1))

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.FakeLogicalNeuroVol", space.FakeLogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("origin.FakeSpace", origin.FakeSpace, envir=.GlobalEnv)
  assign("as.array.FakeLogicalNeuroVol", as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  withr::defer({
    rm(FakeSpace, space.FakeLogicalNeuroVol, spacing.FakeSpace,
       origin.FakeSpace, as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  plan <- Plan$new()
  h <- DataHandle$new(initial_stash = list(input = matrix(0, nrow = 1, ncol = 1)),
                      plan = plan, mask_info = list(mask = vol, active_voxels = 1))
  desc <- list(type = "spat.hrbf",
               params = list(sigma0 = 6, levels = 0, radius_factor = 2.5,
                              kernel_type = "gaussian", seed = 1,
                              include_gaussian_derivatives = "first_order"))

  expect_warning(
    neuroarchive:::forward_step.spat.hrbf("spat.hrbf", desc, h),
    regexp = "Derivative-of-Gaussian"
  )
})

test_that("forward_step.spat.hrbf warns for centre steering map", {
  mask <- array(TRUE, dim = c(1,1,1))
  vol <- structure(list(arr = mask), class = "LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(1,1,1), c(1,1,1))

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.FakeLogicalNeuroVol", space.FakeLogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("origin.FakeSpace", origin.FakeSpace, envir=.GlobalEnv)
  assign("as.array.FakeLogicalNeuroVol", as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  withr::defer({
    rm(FakeSpace, space.FakeLogicalNeuroVol, spacing.FakeSpace,
       origin.FakeSpace, as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  tmp <- local_tempfile(fileext = ".h5")
  h5 <- H5File$new(tmp, mode = "w")
  h5$create_dataset("steer", array(0, dim = c(1,1,1)))

  plan <- Plan$new()
  h <- DataHandle$new(initial_stash = list(input = matrix(0, nrow = 1, ncol = 1)),
                      plan = plan, h5 = h5,
                      mask_info = list(mask = vol, active_voxels = 1))
  desc <- list(type = "spat.hrbf",
               params = list(sigma0 = 6, levels = 0, radius_factor = 2.5,
                              kernel_type = "gaussian", seed = 1,
                              centre_steering = list(map_path = "steer")))

  expect_warning(
    neuroarchive:::forward_step.spat.hrbf("spat.hrbf", desc, h),
    regexp = "Centre steering not fully implemented"
  )
  h5$close_all()
})

test_that("forward_step.spat.hrbf warns for differential encoding", {
  mask <- array(TRUE, dim = c(1,1,1))
  vol <- structure(list(arr = mask), class = "LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(1,1,1), c(1,1,1))

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.FakeLogicalNeuroVol", space.FakeLogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("origin.FakeSpace", origin.FakeSpace, envir=.GlobalEnv)
  assign("as.array.FakeLogicalNeuroVol", as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  withr::defer({
    rm(FakeSpace, space.FakeLogicalNeuroVol, spacing.FakeSpace,
       origin.FakeSpace, as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  plan <- Plan$new()
  h <- DataHandle$new(initial_stash = list(input = matrix(0, nrow = 1, ncol = 1)),
                      plan = plan, mask_info = list(mask = vol, active_voxels = 1))
  desc <- list(type = "spat.hrbf",
               params = list(sigma0 = 6, levels = 0, radius_factor = 2.5,
                              kernel_type = "gaussian", seed = 1,
                              use_differential_encoding = TRUE))

  expect_warning(
    neuroarchive:::forward_step.spat.hrbf("spat.hrbf", desc, h),
    regexp = "Differential encoding"
  )
})
