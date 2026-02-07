# This file tests forward_step.spat.hrbf_project and a basic roundtrip using a mock invert_step.

library(testthat)

# Helper neuroim2 replacements
FakeSpace <- function(dim, spacing_v, origin_v=c(0,0,0)) {
  structure(list(dim=dim, spacing=spacing_v, origin=origin_v), class="FakeSpace")
}
space.FakeLogicalNeuroVol <- function(x, ...) attr(x, "space")
spacing.FakeSpace <- function(x, ...) x$spacing
origin.FakeSpace <- function(x, ...) x$origin
as.array.LogicalNeuroVol <- function(x, ...) x$arr


# Test forward_step.spat.hrbf_project -----------------------------------------------------

test_that("forward_step.spat.hrbf_project outputs coefficients and descriptor", {
  mask <- array(TRUE, dim=c(2,2,2))
  vol <- structure(list(arr=mask), class="LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(2,2,2), c(1,1,1))

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.FakeLogicalNeuroVol", space.FakeLogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("origin.FakeSpace", origin.FakeSpace, envir=.GlobalEnv)
  assign("as.array.LogicalNeuroVol", as.array.LogicalNeuroVol, envir=.GlobalEnv)
  withr::defer({
    rm(FakeSpace, space.FakeLogicalNeuroVol, spacing.FakeSpace,
       origin.FakeSpace, as.array.LogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  X <- matrix(seq_len(16), nrow=2, ncol=8)
  plan <- Plan$new()
  h <- DataHandle$new(initial_stash=list(input=X),
                      plan=plan,
                      mask_info=list(mask=vol, active_voxels=8))
  desc <- list(type="spat.hrbf_project",
               params=list(sigma0=6, levels=0, radius_factor=2.5,
                            kernel_type="gaussian", seed=42))

  h2 <- neuroarchive:::forward_step.spat.hrbf_project("spat.hrbf_project", desc, h)

  expect_true(h2$has_key("hrbf_coefficients"))
  dname <- names(h2$plan$descriptors)[1]
  stored_desc <- h2$plan$descriptors[[dname]]
  expect_true(startsWith(stored_desc$params$mask_hash, "sha256:"))
  expect_equal(ncol(h2$stash$hrbf_coefficients), stored_desc$params$k_actual)
})


# Test roundtrip with mock invert_step -----------------------------------------------------

test_that("spat.hrbf descriptor-only roundtrip with quant", {
  set.seed(1)
  arr <- array(runif(8), dim=c(1,1,2,4))
  mask <- array(TRUE, dim=c(1,1,2))
  vol <- structure(list(arr=mask), class="LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(1,1,2), c(1,1,1))

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.FakeLogicalNeuroVol", space.FakeLogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("origin.FakeSpace", origin.FakeSpace, envir=.GlobalEnv)
  assign("as.array.LogicalNeuroVol", as.array.LogicalNeuroVol, envir=.GlobalEnv)
  withr::defer({
    rm(FakeSpace, space.FakeLogicalNeuroVol, spacing.FakeSpace,
       origin.FakeSpace, as.array.LogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  tmp <- local_tempfile(fileext = ".h5")
  write_lna(arr, file=tmp, mask=vol,
            transforms=c("spat.hrbf","quant"),
            transform_params=list(spat.hrbf=list(sigma0=6, levels=0, radius_factor=2.5,
                                                kernel_type="gaussian", seed=1)))

  h <- read_lna(tmp)
  out <- as.matrix(h$stash$input)
  # read_lna returns time x voxels (4 x 2), not original 4D shape
  expect_equal(nrow(out), dim(arr)[4])   # time dimension
  expect_equal(ncol(out), prod(dim(arr)[1:3]))  # spatial voxels
  # Reshape arr to 2D (time x voxels) for comparison
  arr_mat <- matrix(arr, nrow = prod(dim(arr)[1:3]), ncol = dim(arr)[4])
  arr_mat <- t(arr_mat)  # now time x voxels
  expect_lt(mean(abs(out - arr_mat)), 1)
})

