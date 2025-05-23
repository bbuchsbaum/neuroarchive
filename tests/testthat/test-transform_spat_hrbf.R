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
