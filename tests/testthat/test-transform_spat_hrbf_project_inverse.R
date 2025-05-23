library(neuroarchive)
library(withr)

FakeSpace <- function(dim, spacing_v, origin_v=c(0,0,0)) {
  structure(list(dim=dim, spacing=spacing_v, origin=origin_v), class="FakeSpace")
}
space.FakeLogicalNeuroVol <- function(x, ...) attr(x, "space")
spacing.FakeSpace <- function(x, ...) x$spacing
origin.FakeSpace <- function(x, ...) x$origin
as.array.FakeLogicalNeuroVol <- function(x, ...) x$arr

# Basic reconstruction for spat.hrbf_project -----------------------------------

test_that("invert_step.spat.hrbf_project reconstructs dense data", {
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
  desc <- list(type="spat.hrbf_project",
               params=list(sigma0=6, levels=0, radius_factor=2.5,
                            kernel_type="gaussian", seed=1))
  h2 <- neuroarchive:::forward_step.spat.hrbf_project("spat.hrbf_project", desc, h)
  fdesc <- h2$plan$descriptors[[1]]
  coeff <- h2$stash$hrbf_coefficients

  h_inv <- DataHandle$new(initial_stash=list(hrbf_coefficients=coeff),
                          mask_info=list(mask=vol, active_voxels=8))
  out <- neuroarchive:::invert_step.spat.hrbf_project("spat.hrbf_project", fdesc, h_inv)
  expect_true(out$has_key("input_dense_mat"))
  expect_false(out$has_key("hrbf_coefficients"))
  expect_equal(out$stash$input_dense_mat, X)
})
