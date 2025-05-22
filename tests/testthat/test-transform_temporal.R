library(testthat)
#library(neuroarchive)
library(withr)


test_that("temporal transform forward and inverse roundtrip", {
  set.seed(1)
  X_matrix <- matrix(rnorm(40), nrow = 10, ncol = 4) # time x components
  # Reshape to components x 1 x time for core_write
  X <- array(X_matrix, dim = c(ncol(X_matrix), 1, nrow(X_matrix)))
  tmp <- local_tempfile(fileext = ".h5")

  res <- write_lna(X, file = tmp, transforms = "temporal",
                   transform_params = list(temporal = list(n_basis = 10)))
  expect_true(file.exists(tmp))

  h <- read_lna(tmp)
  out_matrix <- h$stash$input # Should be time x components (2D)
  
  expect_equal(dim(out_matrix), dim(X_matrix))
  expect_equal(out_matrix, X_matrix, tolerance = 1e-6)
})


test_that("invert_step.temporal applies time_idx subset", {
  X_matrix <- matrix(seq_len(40), nrow = 10, ncol = 4) # time x components
  # Reshape to components x 1 x time for core_write
  X <- array(X_matrix, dim = c(ncol(X_matrix), 1, nrow(X_matrix)))
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(X, file = tmp, transforms = "temporal",
            transform_params = list(temporal = list(n_basis = 10)))

  h <- read_lna(tmp, time_idx = c(1,5,10))
  out_matrix <- h$stash$input # Should be subsetted_time x components (2D)
  expect_equal(dim(out_matrix), c(3, ncol(X_matrix)))
  expect_equal(out_matrix, X_matrix[c(1,5,10), ])
})


test_that("default_params for temporal loads schema", {
  neuroarchive:::default_param_cache_clear()
  p <- neuroarchive:::default_params("temporal")
  expect_equal(p$kind, "dct")
  expect_equal(p$scope, "global")
  expect_true(is.numeric(p$n_basis))

})

test_that("temporal transform bspline roundtrip", {
  set.seed(1)
  X_matrix <- matrix(rnorm(60), nrow = 15, ncol = 4) # time x components
  # Reshape to components x 1 x time for core_write
  X <- array(X_matrix, dim = c(ncol(X_matrix), 1, nrow(X_matrix)))
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(X, file = tmp, transforms = "temporal",
            transform_params = list(temporal = list(kind = "bspline",
                                                    n_basis = 8,
                                                    order = 3)))
  h <- read_lna(tmp)
  out_matrix <- h$stash$input # Should be time x components (2D)
  expect_equal(dim(out_matrix), dim(X_matrix))
  expect_equal(out_matrix, X_matrix, tolerance = 1e-6)
})



test_that("temporal transform rejects unsupported kind", {
  X_matrix <- matrix(rnorm(10), nrow = 5, ncol = 2) # time x components
  # Reshape to components x 1 x time for core_write
  X <- array(X_matrix, dim = c(ncol(X_matrix), 1, nrow(X_matrix)))
  expect_error(
    core_write(X, transforms = "temporal",
               transform_params = list(temporal = list(kind = "unsupported_kind"))), # Corrected to unsupported_kind
    class = "lna_error_validation",
    regexp = "temporal kind"
  )
})
          
test_that("temporal transform dpss roundtrip", {
  set.seed(1)
  X_matrix <- matrix(rnorm(64), nrow = 16, ncol = 4) # time x components
  # Reshape to components x 1 x time for core_write
  X <- array(X_matrix, dim = c(ncol(X_matrix), 1, nrow(X_matrix)))
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(X, file = tmp, transforms = "temporal",
            transform_params = list(temporal = list(kind = "dpss",
                                                    n_basis = 4,
                                                    time_bandwidth_product = 3,
                                                    n_tapers = 4)))
  h <- read_lna(tmp)
  out_matrix <- h$stash$input # Should be time x components (2D)
  expect_equal(dim(out_matrix), dim(X_matrix))
  expect_equal(out_matrix, X_matrix, tolerance = 1e-6)
})

test_that("temporal transform polynomial roundtrip", {
  set.seed(1)
  X_matrix <- matrix(rnorm(48), nrow = 12, ncol = 4) # time x components
  # Reshape to components x 1 x time for core_write
  X <- array(X_matrix, dim = c(ncol(X_matrix), 1, nrow(X_matrix)))
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(X, file = tmp, transforms = "temporal",
            transform_params = list(temporal = list(kind = "polynomial",
                                                    n_basis = 5)))
  h <- read_lna(tmp)
  out_matrix <- h$stash$input # Should be time x components (2D)
  expect_equal(dim(out_matrix), dim(X_matrix))
  expect_equal(out_matrix, X_matrix, tolerance = 1e-6)
})

test_that("temporal transform wavelet roundtrip", {
  set.seed(1)
  X_matrix <- matrix(rnorm(64), nrow = 16, ncol = 4) # time x components
  # Reshape to components x 1 x time for core_write
  X <- array(X_matrix, dim = c(ncol(X_matrix), 1, nrow(X_matrix)))
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(X, file = tmp, transforms = "temporal",
            transform_params = list(temporal = list(kind = "wavelet",
                                                    wavelet = "d4")))
  h <- read_lna(tmp)
  out_matrix <- h$stash$input # Should be time x components (2D)
  expect_equal(dim(out_matrix), dim(X_matrix))
  expect_equal(out_matrix, X_matrix, tolerance = 1e-6)

})
