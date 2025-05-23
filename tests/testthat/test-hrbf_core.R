library(neuroarchive)

FakeSpace <- function(dim, spacing_v, origin_v = c(0,0,0)) {
  structure(list(dim = dim, spacing = spacing_v, origin = origin_v), class = "FakeSpace")
}
space.FakeLogicalNeuroVol <- function(x, ...) attr(x, "space")
spacing.FakeSpace <- function(x, ...) x$spacing
origin.FakeSpace <- function(x, ...) x$origin
as.array.FakeLogicalNeuroVol <- function(x, ...) x$arr


test_that("hrbf helpers match transform outputs", {
  mask <- array(TRUE, dim = c(1,1,2))
  vol <- structure(list(arr = mask), class = "LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(1,1,2), c(1,1,1))

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.FakeLogicalNeuroVol", space.FakeLogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("origin.FakeSpace", origin.FakeSpace, envir=.GlobalEnv)
  assign("as.array.FakeLogicalNeuroVol", as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  withr::defer({
    rm(FakeSpace, space.FakeLogicalNeuroVol, spacing.FakeSpace,
       origin.FakeSpace, as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  X <- matrix(1:4, nrow = 2)
  params <- list(sigma0 = 6, levels = 0, radius_factor = 2.5,
                 kernel_type = "gaussian", seed = 1)

  B1 <- hrbf_generate_basis(params, vol)
  B2 <- neuroarchive:::hrbf_basis_from_params(params, vol)
  expect_equal(B1, B2)

  coeff_direct <- hrbf_project_matrix(X, vol, params)
  dense_direct <- hrbf_reconstruct_matrix(coeff_direct, vol, params)

  plan <- Plan$new()
  h <- DataHandle$new(initial_stash = list(input_dense_mat = X), plan = plan,
                      mask_info = list(mask = vol, active_voxels = 2))
  desc <- list(type = "spat.hrbf_project", params = params)
  out <- neuroarchive:::forward_step.spat.hrbf_project("spat.hrbf_project", desc, h)
  fdesc <- out$plan$descriptors[[1]]
  coeff_tr <- out$stash$hrbf_coefficients

  expect_equal(coeff_direct, coeff_tr)

  h_inv <- DataHandle$new(initial_stash = list(hrbf_coefficients = coeff_tr),
                          mask_info = list(mask = vol, active_voxels = 2))
  out2 <- neuroarchive:::invert_step.spat.hrbf_project("spat.hrbf_project", fdesc, h_inv)
  dense_tr <- out2$stash$input_dense_mat

  expect_equal(dense_direct, dense_tr)
})
