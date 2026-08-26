library(testthat)
library(neuroarchive)

pr378_geometry <- function(dims, affine, n_volumes) {
  list(
    dims = as.integer(dims),
    transform = affine,
    spatial_units = "mm",
    voxel_ordering = "first-axis-fastest",
    n_volumes = as.integer(n_volumes)
  )
}

test_that("matrix voxel dequantization applies parameters by column", {
  x <- cbind(
    seq(-1, 1, length.out = 17),
    seq(-100, 100, length.out = 17),
    seq(250, 500, length.out = 17),
    seq(-0.01, 0.01, length.out = 17)
  )
  quantized <- quantize_data(
    x,
    bits = 8L,
    method = "range",
    scope = "voxel",
    return_stats = TRUE
  )

  expected <- sweep(quantized$q, 2, quantized$scale, FUN = "*")
  expected <- sweep(expected, 2, quantized$offset, FUN = "+")
  restored <- dequantize_data(quantized)

  expect_equal(restored, expected, tolerance = 1e-12)
  expect_equal(
    quantized$mse,
    mean((as.numeric(x) - as.numeric(expected))^2),
    tolerance = 1e-12
  )
  expect_lt(quantized$mse, 0.1)
})

test_that("hard-threshold sparse PCA scores match retained loadings", {
  set.seed(378)
  x <- matrix(rnorm(60 * 12), nrow = 60, ncol = 12)
  fit <- sparse_pca(
    x,
    n_components = 4L,
    sparsity_type = "hard_threshold",
    sparsity_param = 0.25,
    standardize = FALSE
  )

  expect_true(any(fit$loadings == 0))
  expect_equal(fit$scores, x %*% t(fit$loadings), tolerance = 1e-12)

  full_fit <- sparse_pca(
    x,
    n_components = 4L,
    sparsity_type = "hard_threshold",
    sparsity_param = 1,
    standardize = FALSE
  )
  expect_equal(full_fit$scores, x %*% t(full_fit$loadings), tolerance = 1e-12)
})

test_that("LNANeuroVec preserves the complete stored affine", {
  skip_if_not_installed("neuroim2")

  dims <- c(2L, 2L, 2L)
  n_time <- 4L
  affine <- rbind(
    c(0, -2, 0.25, 10),
    c(3, 0, 0, -5),
    c(0, 0, 4, 2),
    c(0, 0, 0, 1)
  )
  x <- array(rnorm(prod(dims) * n_time), dim = c(dims, n_time))
  file <- local_tempfile(fileext = ".lna.h5")
  write_lna(
    x,
    file,
    transforms = "quant",
    run_geometry = list(
      `run-01` = pr378_geometry(dims, affine, n_time)
    ),
    geometry_policy = "fmri",
    checksum = "none"
  )

  lazy <- LNANeuroVec(file)
  withr::defer(lna_close(lazy))

  expect_equal(neuroim2::trans(neuroim2::space(lazy)), affine)
  expect_equal(neuroim2::trans(neuroim2::space(lazy@mask)), affine)
})

test_that("as_neurovec intersects a stored mask with a strict ROI", {
  skip_if_not_installed("neuroim2")

  dims <- c(2L, 2L, 2L)
  n_time <- 4L
  affine <- diag(c(2, 2, 2, 1))
  x <- array(seq_len(prod(dims) * n_time), dim = c(dims, n_time))
  archive_mask <- array(FALSE, dim = dims)
  archive_mask[c(1, 2, 4, 5, 7)] <- TRUE
  roi <- array(FALSE, dim = dims)
  roi[c(2, 5, 7)] <- TRUE
  file <- local_tempfile(fileext = ".lna.h5")
  write_lna(
    x,
    file,
    transforms = "quant",
    mask = archive_mask,
    run_geometry = list(
      `run-01` = pr378_geometry(dims, affine, n_time)
    ),
    geometry_policy = "fmri",
    checksum = "none"
  )

  handle <- read_lna(file, roi_mask = roi, validate = FALSE)
  restored <- as_neurovec(handle, file = file)

  expect_s4_class(restored, "SparseNeuroVec")
  expect_equal(dim(restored), c(dims, n_time))
  expect_equal(sum(as.logical(restored@mask)), sum(archive_mask & roi))
  expect_equal(
    as.matrix(restored@data),
    t(neuroarchive:::as_dense_mat(handle$stash$input)),
    tolerance = 1e-12
  )
})

test_that("failed atomic publication retains the completed temporary archive", {
  target <- local_tempfile(fileext = ".lna.h5")
  captured <- NULL

  error <- testthat::with_mocked_bindings(
    code = {
      expect_error(
        write_lna(
          array(seq_len(16), dim = c(2, 2, 2, 2)),
          target,
          checksum = "none"
        ),
        class = "lna_error_write"
      )
    },
    file.rename = function(from, to) {
      captured <<- from
      FALSE
    },
    file.copy = function(from, to, overwrite = FALSE, ...) FALSE,
    .package = "base"
  )

  expect_identical(error$temporary_file, captured)
  expect_identical(error$target_file, target)
  expect_true(file.exists(error$temporary_file))
  expect_false(file.exists(target))
  withr::defer(unlink(error$temporary_file, force = TRUE))
})
