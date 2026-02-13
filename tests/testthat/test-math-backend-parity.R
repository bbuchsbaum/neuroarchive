library(testthat)
library(neuroarchive)
library(neuroim2)

make_mask <- function(mask_arr, spacing_v = c(1, 1, 1)) {
  LogicalNeuroVol(mask_arr, NeuroSpace(dim(mask_arr), spacing = spacing_v))
}

compare_sparse <- function(a, b, tol = 1e-8) {
  expect_equal(dim(a), dim(b))
  expect_equal(Matrix::nnzero(a), Matrix::nnzero(b))
  expect_equal(as.matrix(a), as.matrix(b), tolerance = tol)
}

test_that("HRBF backend parity: seed path across mask types", {
  skip_if_not_installed("fmrilatent")
  old_opts <- options(
    lna.hrbf.use_rcpp_helpers = FALSE,
    fmrilatent.hrbf.use_rcpp = FALSE
  )
  on.exit(options(old_opts), add = TRUE)
  lna_options(engine.math = "internal")
  on.exit(lna_options(engine.math = "internal"), add = TRUE)

  set.seed(11)
  mask_full <- array(TRUE, dim = c(4, 4, 4))
  mask_sparse <- array(FALSE, dim = c(4, 4, 4))
  mask_sparse[sample(length(mask_sparse), 18)] <- TRUE
  mask_disconnected <- array(FALSE, dim = c(6, 6, 6))
  mask_disconnected[1:2, 1:2, 1:2] <- TRUE
  mask_disconnected[5:6, 5:6, 5:6] <- TRUE
  mask_tiny <- array(TRUE, dim = c(3, 3, 3))

  masks <- list(
    full = make_mask(mask_full),
    sparse = make_mask(mask_sparse),
    disconnected = make_mask(mask_disconnected),
    tiny = make_mask(mask_tiny),
    anisotropic = make_mask(mask_full, spacing_v = c(1, 2, 3))
  )

  params <- list(
    sigma0 = 3,
    levels = 2L,
    radius_factor = 2.5,
    kernel_type = "wendland_c4",
    seed = 23L,
    num_extra_fine_levels = 1L
  )

  for (nm in names(masks)) {
    mask_obj <- masks[[nm]]
    B_internal <- neuroarchive:::hrbf_basis_from_params_internal(params, mask_obj)
    lna_options(engine.math = "fmrilatent")
    B_fmrilatent <- neuroarchive:::hrbf_basis_from_params(params, mask_obj)
    lna_options(engine.math = "internal")
    compare_sparse(B_internal, B_fmrilatent, tol = 1e-8)
  }
})

test_that("HRBF centre generation parity through backend boundary", {
  skip_if_not_installed("fmrilatent")
  old_opts <- options(
    lna.hrbf.use_rcpp_helpers = FALSE,
    fmrilatent.hrbf.use_rcpp = FALSE
  )
  on.exit(options(old_opts), add = TRUE)
  lna_options(engine.math = "internal")
  on.exit(lna_options(engine.math = "internal"), add = TRUE)

  set.seed(41)
  mask_sparse <- array(FALSE, dim = c(5, 4, 3))
  mask_sparse[sample(length(mask_sparse), 19)] <- TRUE
  mask_obj <- make_mask(mask_sparse, spacing_v = c(1, 2, 1.5))

  params <- list(
    sigma0 = 3,
    levels = 2L,
    radius_factor = 2.5,
    seed = 17L,
    num_extra_fine_levels = 1L
  )

  lna_options(engine.math = "internal")
  cs_internal <- neuroarchive:::na_hrbf_centres_from_params(params, mask_obj)
  lna_options(engine.math = "fmrilatent")
  cs_fmrilatent <- neuroarchive:::na_hrbf_centres_from_params(params, mask_obj)

  expect_equal(cs_fmrilatent$centres_stored, cs_internal$centres_stored)
  expect_equal(cs_fmrilatent$centres, cs_internal$centres, tolerance = 1e-10)
  expect_equal(cs_fmrilatent$sigmas, cs_internal$sigmas, tolerance = 1e-10)
  expect_equal(cs_fmrilatent$level_vec, cs_internal$level_vec)
})

test_that("HRBF backend parity: stored centers path", {
  skip_if_not_installed("fmrilatent")
  old_opts <- options(
    lna.hrbf.use_rcpp_helpers = FALSE,
    fmrilatent.hrbf.use_rcpp = FALSE
  )
  on.exit(options(old_opts), add = TRUE)
  lna_options(engine.math = "internal")
  on.exit(lna_options(engine.math = "internal"), add = TRUE)

  mask_arr <- array(TRUE, dim = c(4, 4, 4))
  mask_obj <- make_mask(mask_arr, spacing_v = c(1, 2, 1.5))

  tf <- tempfile(fileext = ".h5")
  h5 <- hdf5r::H5File$new(tf, mode = "w")
  on.exit({
    try(h5$close_all(), silent = TRUE)
    if (file.exists(tf)) unlink(tf)
  }, add = TRUE)
  h5$create_group("aux")
  centres <- matrix(c(
    0, 0, 0,
    2, 2, 2,
    1.5, 3.0, 0.5
  ), ncol = 3, byrow = TRUE)
  sigmas <- c(3.0, 1.5, 0.75)
  h5[["/aux"]]$create_dataset("centres", robj = centres)
  h5[["/aux"]]$create_dataset("sigmas", robj = sigmas)

  params <- list(
    sigma0 = 3,
    levels = 2L,
    radius_factor = 2.5,
    kernel_type = "gaussian",
    centres_stored = TRUE,
    centres_path = "/aux/centres",
    sigma_vec_path = "/aux/sigmas"
  )

  root <- h5[["/"]]
  B_internal <- neuroarchive:::hrbf_basis_from_params_internal(params, mask_obj, h5_root = root)
  lna_options(engine.math = "fmrilatent")
  B_fmrilatent <- neuroarchive:::hrbf_basis_from_params(params, mask_obj, h5_root = root)
  compare_sparse(B_internal, B_fmrilatent, tol = 1e-8)
})

test_that("Haar backend parity: coefficients and roundtrip across masks", {
  skip_if_not_installed("fmrilatent")
  old_opts <- options(
    lna.hwt.use_rcpp = FALSE,
    lna.hwt.use_rcpp_batch = FALSE,
    fmrilatent.haar.use_rcpp = FALSE,
    fmrilatent.haar.use_rcpp_batch = FALSE
  )
  on.exit(options(old_opts), add = TRUE)
  lna_options(engine.math = "internal")
  on.exit(lna_options(engine.math = "internal"), add = TRUE)

  set.seed(101)
  mask_full <- array(TRUE, dim = c(4, 4, 4))
  mask_sparse <- array(FALSE, dim = c(4, 4, 4))
  mask_sparse[sample(length(mask_sparse), 22)] <- TRUE
  mask_disconnected <- array(FALSE, dim = c(6, 6, 6))
  mask_disconnected[1:2, 1:2, 1:2] <- TRUE
  mask_disconnected[5:6, 5:6, 5:6] <- TRUE
  mask_tiny <- array(TRUE, dim = c(3, 3, 3))

  masks <- list(full = mask_full, sparse = mask_sparse, disconnected = mask_disconnected, tiny = mask_tiny)

  for (nm in names(masks)) {
    mask_arr <- masks[[nm]]
    X <- matrix(rnorm(5 * sum(mask_arr)), nrow = 5)
    levels <- 2L

    lna_options(engine.math = "internal")
    fw_internal <- neuroarchive:::perform_haar_lift_analysis(X, mask_arr, levels, z_order_seed = 42L)
    reco_internal <- neuroarchive:::perform_haar_lift_synthesis(fw_internal, mask_arr, levels, z_order_seed = 42L)

    lna_options(engine.math = "fmrilatent")
    fw_fmrilatent <- neuroarchive:::perform_haar_lift_analysis(X, mask_arr, levels, z_order_seed = 42L)
    reco_fmrilatent <- neuroarchive:::perform_haar_lift_synthesis(fw_fmrilatent, mask_arr, levels, z_order_seed = 42L)

    expect_equal(fw_fmrilatent$root, fw_internal$root, tolerance = 1e-10)
    for (lv in seq_len(levels)) {
      expect_equal(fw_fmrilatent$detail[[lv]], fw_internal$detail[[lv]], tolerance = 1e-10)
    }
    expect_equal(reco_fmrilatent, reco_internal, tolerance = 1e-10)
  }
})
