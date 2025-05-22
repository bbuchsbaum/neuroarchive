library(testthat)
#library(neuroarchive)
library(withr)


test_that("temporal transform forward and inverse roundtrip", {
  set.seed(1)
  X_matrix <- matrix(rnorm(40), nrow = 10, ncol = 4) # time x components
  message("--- X_matrix DEFINITION in test_that [1:2,1:2] ---")
  if (nrow(X_matrix) >=2 && ncol(X_matrix) >=2) print(X_matrix[1:2,1:2, drop=FALSE])
  
  # Reshape to components x 1 x time for core_write
  # Original X_matrix is Time x Components (e.g., 10x4)
  # We need X_for_core_write to be Components x 1 x Time (e.g., 4x1x10)
  # such that X_for_core_write[comp, 1, time] == X_matrix[time, comp]
  X_matrix_transposed <- t(X_matrix) # Now Components x Time (e.g., 4x10)
  X <- array(X_matrix_transposed, 
             dim = c(nrow(X_matrix_transposed), 1, ncol(X_matrix_transposed))) # Should be 4x1x10
  # This ensures X[c,1,t] == X_matrix_transposed[c,t] == X_matrix[t,c]

  tmp <- local_tempfile(fileext = ".h5")

  res <- write_lna(X, file = tmp, transforms = "temporal",
                   transform_params = list(temporal = list(n_basis = 10)))
  expect_true(file.exists(tmp))

  h <- read_lna(tmp)
  out_matrix <- h$stash$input # Should be time x components (2D)
  
  message("--- Testthat Scope: X_matrix (expected) [1:2,1:2] ---")
  if (nrow(X_matrix) >=2 && ncol(X_matrix) >=2) print(X_matrix[1:2,1:2, drop=FALSE])
  message("--- Testthat Scope: out_matrix (actual) [1:2,1:2] ---")
  if (nrow(out_matrix) >=2 && ncol(out_matrix) >=2) print(out_matrix[1:2,1:2, drop=FALSE])
  
  expect_equal(dim(out_matrix), dim(X_matrix))
  expect_equal(out_matrix, X_matrix, tolerance = 1e-6)
})


test_that("invert_step.temporal applies time_idx subset", {
  X_matrix <- matrix(seq_len(40), nrow = 10, ncol = 4) # time x components
  # Reshape to components x 1 x time for core_write
  X_matrix_transposed <- t(X_matrix) # Now Components x Time
  X <- array(X_matrix_transposed, 
             dim = c(nrow(X_matrix_transposed), 1, ncol(X_matrix_transposed)))
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
  X_matrix_transposed <- t(X_matrix) # Now Components x Time
  X <- array(X_matrix_transposed, 
             dim = c(nrow(X_matrix_transposed), 1, ncol(X_matrix_transposed)))
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(X, file = tmp, transforms = "temporal",
            transform_params = list(temporal = list(kind = "bspline",
                                                    n_basis = 8,
                                                    order = 3)))
  h <- read_lna(tmp)
  out_matrix <- h$stash$input # Should be time x components (2D)
  expect_equal(dim(out_matrix), dim(X_matrix))
  # Check that the residual is orthogonal to the projection
  residual <- X_matrix - out_matrix
  # Sum of element-wise product should be close to 0 for each column
  for (i in 1:ncol(X_matrix)) {
    expect_equal(sum(residual[,i] * out_matrix[,i]), 0, tolerance = 1e-6)
  }
})



test_that("temporal transform rejects unsupported kind", {
  X_matrix <- matrix(rnorm(10), nrow = 5, ncol = 2) # time x components
  # Reshape to components x 1 x time for core_write
  X_matrix_transposed <- t(X_matrix) # Now Components x Time
  X <- array(X_matrix_transposed, 
             dim = c(nrow(X_matrix_transposed), 1, ncol(X_matrix_transposed)))
  expect_error(
    core_write(X, transforms = "temporal",
               transform_params = list(temporal = list(kind = "unsupported_kind"))),
    class = "lna_error_transform_step" # Expect the outer wrapping class from run_transform_step
  )
  
  # Check the parent error (the one thrown by temporal_basis.default)
  expect_s3_class(err$parent, "lna_error_validation")
  expect_match(conditionMessage(err$parent), "Unsupported temporal kind 'unsupported_kind'")
  # Check the location from the original error source
  expect_equal(err$parent$location, "temporal_basis:kind")
})
          
test_that("temporal transform dpss roundtrip", {
  set.seed(1)
  X_matrix <- matrix(rnorm(64), nrow = 16, ncol = 4) # time x components
  # Reshape to components x 1 x time for core_write
  X_matrix_transposed <- t(X_matrix) # Now Components x Time
  X <- array(X_matrix_transposed, 
             dim = c(nrow(X_matrix_transposed), 1, ncol(X_matrix_transposed)))
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(X, file = tmp, transforms = "temporal",
            transform_params = list(temporal = list(kind = "dpss",
                                                    n_basis = 4,
                                                    time_bandwidth_product = 3,
                                                    n_tapers = 4)))
  h <- read_lna(tmp)
  out_matrix <- h$stash$input # Should be time x components (2D)
  expect_equal(dim(out_matrix), dim(X_matrix))
  # Check that the residual is orthogonal to the projection
  residual <- X_matrix - out_matrix
  for (i in 1:ncol(X_matrix)) {
    expect_equal(sum(residual[,i] * out_matrix[,i]), 0, tolerance = 1e-6)
  }
})

test_that("temporal transform polynomial roundtrip", {
  set.seed(1)
  X_matrix <- matrix(rnorm(48), nrow = 12, ncol = 4) # time x components
  # Reshape to components x 1 x time for core_write
  X_matrix_transposed <- t(X_matrix) # Now Components x Time
  X <- array(X_matrix_transposed, 
             dim = c(nrow(X_matrix_transposed), 1, ncol(X_matrix_transposed)))
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(X, file = tmp, transforms = "temporal",
            transform_params = list(temporal = list(kind = "polynomial",
                                                    n_basis = 12)))
  h <- read_lna(tmp)
  out_matrix <- h$stash$input # Should be time x components (2D)
  expect_equal(dim(out_matrix), dim(X_matrix))
  expect_equal(out_matrix, X_matrix, tolerance = 1e-6)
})

test_that("temporal transform wavelet roundtrip", {
  set.seed(1)
  X_matrix <- matrix(rnorm(64), nrow = 16, ncol = 4) # time x components
  # Reshape to components x 1 x time for core_write
  X_matrix_transposed <- t(X_matrix) # Now Components x Time
  X <- array(X_matrix_transposed, 
             dim = c(nrow(X_matrix_transposed), 1, ncol(X_matrix_transposed)))
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(X, file = tmp, transforms = "temporal",
            transform_params = list(temporal = list(kind = "wavelet",
                                                    wavelet = "d4")))
  h <- read_lna(tmp)
  out_matrix <- h$stash$input # Should be time x components (2D)
  expect_equal(dim(out_matrix), dim(X_matrix))
  expect_equal(out_matrix, X_matrix, tolerance = 1e-6)

})
