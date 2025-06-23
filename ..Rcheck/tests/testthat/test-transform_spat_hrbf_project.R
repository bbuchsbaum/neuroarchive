# This file tests forward_step.spat.hrbf_project and a basic roundtrip using a mock invert_step.

library(testthat)

# Helper neuroim2 replacements
FakeSpace <- function(dim, spacing_v, origin_v=c(0,0,0)) {
  structure(list(dim=dim, spacing=spacing_v, origin=origin_v), class="FakeSpace")
}
space.FakeLogicalNeuroVol <- function(x, ...) attr(x, "space")
spacing.FakeSpace <- function(x, ...) x$spacing
origin.FakeSpace <- function(x, ...) x$origin
as.array.FakeLogicalNeuroVol <- function(x, ...) x$arr


# Test forward_step.spat.hrbf_project -----------------------------------------------------

test_that("forward_step.spat.hrbf_project outputs coefficients and descriptor", {
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

  X <- matrix(1:8, nrow=2)
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
  assign("as.array.FakeLogicalNeuroVol", as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  withr::defer({
    rm(FakeSpace, space.FakeLogicalNeuroVol, spacing.FakeSpace,
       origin.FakeSpace, as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  invert_mock <- function(type, desc, handle) {
    p <- desc$params
    sigma0 <- p$sigma0
    levels <- p$levels
    radius_factor <- p$radius_factor
    seed <- p$seed
    mask_vol <- handle$mask_info$mask
    voxel_to_world <- function(vox_mat) {
      spc <- tryCatch(space(mask_vol), error=function(e) NULL)
      spacing_vec <- tryCatch(spacing(spc), error=function(e) c(1,1,1))
      origin_vec <- tryCatch(origin(spc), error=function(e) c(0,0,0))
      sweep(vox_mat - 1, 2, spacing_vec, `*`) +
        matrix(origin_vec, nrow(vox_mat), 3, byrow=TRUE)
    }
    centres_list <- list(); sigs <- numeric(); level_vec <- integer()
    for (j in seq_len(levels + 1L) - 1L) {
      sigma_j <- sigma0 / (2^j)
      r_j <- radius_factor * sigma_j
      vox_centres <- neuroarchive:::poisson_disk_sample_neuroim2(mask_vol, r_j, seed + j)
      if (nrow(vox_centres) > 0) {
        centres_list[[length(centres_list)+1L]] <- voxel_to_world(vox_centres)
        n_new <- nrow(vox_centres)
        sigs <- c(sigs, rep(sigma_j, n_new))
        level_vec <- c(level_vec, rep(j, n_new))
      }
    }
    C_total <- if (length(centres_list) > 0) do.call(rbind, centres_list) else matrix(numeric(0), ncol=3)
    sigma_vec <- sigs

    mask_arr <- as.array(mask_vol)
    mask_idx <- which(mask_arr)
    vox_coords <- which(mask_arr, arr.ind=TRUE)
    mask_coords_world <- voxel_to_world(vox_coords)
    k <- nrow(C_total)
    n_vox <- length(mask_arr)
    if (k > 0) {
      i_idx <- integer(); j_idx <- integer(); x_val <- numeric()
      for (kk in seq_len(k)) {
        atom <- neuroarchive:::generate_hrbf_atom(mask_coords_world, mask_idx,
                                                 C_total[kk,], sigma_vec[kk],
                                                 level_vec[kk], levels, p)
        i_idx <- c(i_idx, rep.int(kk, length(atom$indices)))
        j_idx <- c(j_idx, atom$indices)
        x_val <- c(x_val, atom$values)
      }
      B_final <- Matrix::sparseMatrix(i=i_idx, j=j_idx, x=x_val,
                                      dims=c(k, n_vox))
    } else {
      B_final <- Matrix::sparseMatrix(i=integer(), j=integer(), x=numeric(),
                                      dims=c(0, n_vox))
    }

    coeff_key <- desc$outputs[[1]] %||% "coefficients_hrbf"
    input_key  <- desc$inputs[[1]] %||% "input"
    C <- handle$get_inputs(coeff_key)[[1]]
    dense_mat <- tcrossprod(C, B_final)
    arr_out <- array(as.numeric(t(dense_mat)), dim=c(dim(mask_vol), nrow(C)))
    handle$update_stash(keys = coeff_key,
                        new_values = setNames(list(arr_out), input_key))
  }
  assign("invert_step.spat.hrbf", invert_mock, envir=.GlobalEnv)
  withr::defer(rm(invert_step.spat.hrbf, envir=.GlobalEnv), envir = parent.frame())

  tmp <- local_tempfile(fileext = ".h5")
  write_lna(arr, file=tmp, mask=vol,
            transforms=c("spat.hrbf","quant"),
            transform_params=list(spat.hrbf=list(sigma0=6, levels=0, radius_factor=2.5,
                                                kernel_type="gaussian", seed=1)))

  h <- read_lna(tmp)
  out <- h$stash$input
  expect_equal(dim(out), dim(arr))
  expect_lt(mean(abs(out - arr)), 1)
})

