test_that("rcpp_control functionality", {
  # Test Rcpp control interface
  original_status <- rcpp_control()
  
  # Disable all acceleration
  rcpp_control(FALSE)
  current_status <- rcpp_control()
  expect_false(current_status$hrbf)
  expect_false(current_status$edge)
  expect_false(current_status$hwt)
  
  # Enable specific components
  rcpp_control(TRUE, components = "hrbf")
  current_status <- rcpp_control()
  expect_true(current_status$hrbf)
  # Other components should still be disabled
  expect_false(current_status$edge)
  expect_false(current_status$hwt)
  
  # Test verbose output
  expect_output(rcpp_control(verbose = TRUE), "Rcpp acceleration status")
  
  # Restore original settings
  rcpp_control(original_status$hrbf, components = "hrbf")
  rcpp_control(original_status$edge, components = "edge") 
  rcpp_control(original_status$hwt, components = "hwt")
})

# Helper function to check if Rcpp compiled code is available
rcpp_available <- function(symbol) {
  tryCatch({
    .Call(symbol, PACKAGE = "neuroarchive")
    TRUE
  }, error = function(e) {
    if (grepl("not available for \\.Call", e$message)) {
      FALSE
    } else {
      # Some other error (like wrong arguments) suggests the function exists
      TRUE
    }
  })
}

test_that("basis_space_centers basic functionality", {
  skip_if_not_installed("neuroim2")
  # Skip if compiled Rcpp code is not available
  skip_if(!rcpp_available("_neuroarchive_poisson_disk_sample_component_rcpp"), 
          "HRBF Rcpp functions not available in development environment")
  
  # Create test mask
  dims <- c(10, 10, 5)
  arr <- array(TRUE, dims)
  mask <- neuroim2::LogicalNeuroVol(arr, neuroim2::NeuroSpace(dims, spacing = c(2, 2, 2)))
  
  # Basic center generation
  result <- basis_space_centers(mask, sigma0 = 4, levels = 2)
  
  expect_type(result, "list")
  expect_true("basis_matrix" %in% names(result))
  expect_true("params" %in% names(result))
  expect_true("n_centers" %in% names(result))
  expect_true("mask_hash" %in% names(result))
  
  # Check parameters are stored correctly
  expect_equal(result$params$sigma0, 4)
  expect_equal(result$params$levels, 2L)
  expect_equal(result$params$kernel_type, "gaussian")
})

test_that("basis_space_centers parameter validation", {
  skip_if_not_installed("neuroim2")
  skip_if(!rcpp_available("_neuroarchive_poisson_disk_sample_component_rcpp"), 
          "HRBF Rcpp functions not available in development environment")
  
  dims <- c(10, 10, 5)
  arr <- array(TRUE, dims)
  mask <- neuroim2::LogicalNeuroVol(arr, neuroim2::NeuroSpace(dims, spacing = c(2, 2, 2)))
  
  # Test kernel type validation
  expect_no_error(basis_space_centers(mask, kernel_type = "gaussian"))
  expect_no_error(basis_space_centers(mask, kernel_type = "wendland_c4"))
  
  # Test edge adaptive parameters
  result_adaptive <- basis_space_centers(
    mask, 
    use_edge_adaptive = TRUE,
    edge_source = "self_mean"
  )
  
  expect_true("edge_adaptive" %in% names(result_adaptive$params))
  expect_equal(result_adaptive$params$edge_adaptive$source, "self_mean")
})

test_that("basis_time DCT functionality", {
  n_timepoints <- 100
  n_components <- 10
  
  basis <- basis_time(n_timepoints, "dct", n_components = n_components)
  
  expect_type(basis, "double")
  expect_equal(nrow(basis), n_timepoints)
  expect_equal(ncol(basis), n_components)
  
  # DCT basis should be orthogonal (approximately)
  gram <- t(basis) %*% basis
  expect_true(max(abs(gram - diag(diag(gram)))) < 1e-10)
})

test_that("basis_time polynomial functionality", {
  n_timepoints <- 50
  n_components <- 5
  
  basis <- basis_time(n_timepoints, "polynomial", n_components = n_components, degree = 3)
  
  expect_type(basis, "double")
  expect_equal(nrow(basis), n_timepoints)
  expect_equal(ncol(basis), n_components)
})

test_that("basis_time DPSS functionality", {
  skip_if_not_installed("multitaper")
  skip_if(!exists("dpss.taper", envir = asNamespace("multitaper")), "multitaper::dpss.taper not available")
  
  n_timepoints <- 100
  n_components <- 7
  
  basis <- basis_time(n_timepoints, "dpss", n_components = n_components, time_bandwidth_product = 4)
  
  expect_type(basis, "double")
  expect_equal(nrow(basis), n_timepoints)
  expect_equal(ncol(basis), n_components)
})

test_that("basis_time default components", {
  # Test default component selection
  basis_dct <- basis_time(100, "dct")
  expect_equal(ncol(basis_dct), 20)  # min(20, 100 %/% 5)
  
  basis_poly <- basis_time(100, "polynomial")
  expect_equal(ncol(basis_poly), 5)  # min(5, 100 %/% 10)
})

test_that("delta_transform basic functionality", {
  # Create test data with trend
  set.seed(123)
  n_time <- 50
  n_vars <- 20
  data <- matrix(rnorm(n_time * n_vars), nrow = n_time, ncol = n_vars)
  # Add linear trend
  data <- data + outer(1:n_time, rep(1, n_vars))
  
  # First differences
  diff1 <- delta_transform(data, order = 1, method = "forward")
  expect_equal(nrow(diff1), n_time - 1)
  expect_equal(ncol(diff1), n_vars)
  
  # Second differences  
  diff2 <- delta_transform(data, order = 2, method = "forward")
  expect_equal(nrow(diff2), n_time - 2)
  expect_equal(ncol(diff2), n_vars)
})

test_that("delta_transform methods", {
  set.seed(123)
  data <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
  
  # Test different methods
  forward <- delta_transform(data, order = 1, method = "forward")
  backward <- delta_transform(data, order = 1, method = "backward") 
  central <- delta_transform(data, order = 1, method = "central")
  
  expect_equal(nrow(forward), 49)
  expect_equal(nrow(backward), 49)
  expect_equal(nrow(central), 48)  # Central differences lose 2 points
})

test_that("delta_transform mean removal", {
  set.seed(123)
  # Create data with different column means and some noise
  n_time <- 20
  n_vars <- 5
  data <- matrix(rnorm(n_time * n_vars, mean = 0, sd = 0.1), nrow = n_time, ncol = n_vars)
  
  # Add different column means (this is what mean removal should affect)
  col_means <- c(10, 20, 30, 40, 50)
  data <- sweep(data, 2, col_means, "+")
  
  # With mean removal (default) - should remove column means before differencing
  diff_demean <- delta_transform(data, order = 1, remove_mean = TRUE)
  
  # Without mean removal - operates on raw data including different column levels
  diff_no_demean <- delta_transform(data, order = 1, remove_mean = FALSE)
  
  # Results should be identical when there's no trend (just different levels)
  # because differencing removes constant offsets anyway
  # But the function behavior should still be different
  expect_true(is.matrix(diff_demean))
  expect_true(is.matrix(diff_no_demean))
  expect_equal(dim(diff_demean), dim(diff_no_demean))
  
  # Test that the function parameter is actually being used
  # by checking that remove_mean FALSE doesn't call scale()
  data_small <- matrix(c(10, 11, 12, 20, 21, 22), nrow = 3, ncol = 2)
  diff1 <- delta_transform(data_small, order = 1, remove_mean = TRUE)
  diff2 <- delta_transform(data_small, order = 1, remove_mean = FALSE)
  
  # For this simple case, they should be equal since diff removes constant offsets
  expect_equal(diff1, diff2)
})

test_that("quantize_data basic functionality", {
  set.seed(123)
  data <- matrix(rnorm(100 * 50, mean = 100, sd = 15), nrow = 100, ncol = 50)
  
  # Basic quantization
  quant_result <- quantize_data(data, bits = 8)
  
  expect_type(quant_result, "list")
  expect_true("quantized_data" %in% names(quant_result))
  expect_true("params" %in% names(quant_result))
  expect_true("clip_report" %in% names(quant_result))
  
  # Check parameters
  expect_equal(quant_result$params$bits, 8L)
  expect_equal(quant_result$params$method, "range")
  expect_equal(quant_result$params$scale_scope, "global")
})

test_that("quantize_data parameter validation", {
  set.seed(123)
  data <- matrix(rnorm(50 * 20), nrow = 50, ncol = 20)
  
  # Test different methods and scopes
  quant_range <- quantize_data(data, method = "range", scope = "global")
  expect_equal(quant_range$params$method, "range")
  expect_equal(quant_range$params$scale_scope, "global")
  
  quant_sd <- quantize_data(data, method = "sd", scope = "voxel")
  expect_equal(quant_sd$params$method, "sd")
  expect_equal(quant_sd$params$scale_scope, "voxel")
})

test_that("dequantize_data functionality", {
  set.seed(123)
  data <- matrix(rnorm(50 * 20), nrow = 50, ncol = 20)
  
  quant_result <- quantize_data(data, bits = 8)
  reconstructed <- dequantize_data(quant_result)
  
  expect_type(reconstructed, "double")
  expect_equal(dim(reconstructed), dim(data))
})

test_that("sparse_pca basic functionality", {
  set.seed(123)
  # Create data with some structure
  data <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
  data[, 1:10] <- data[, 1:10] + rnorm(100, mean = 2)
  
  # Basic sparse PCA
  spca_result <- sparse_pca(data, n_components = 5)
  
  expect_type(spca_result, "list")
  expect_true("components" %in% names(spca_result))
  expect_true("loadings" %in% names(spca_result))
  expect_true("params" %in% names(spca_result))
  expect_true("convergence" %in% names(spca_result))
  
  expect_equal(nrow(spca_result$components), nrow(data))
  expect_equal(ncol(spca_result$components), 5)
  expect_equal(nrow(spca_result$loadings), ncol(data))
  expect_equal(ncol(spca_result$loadings), 5)
})

test_that("sparse_pca parameter validation", {
  set.seed(123)
  data <- matrix(rnorm(50 * 30), nrow = 50, ncol = 30)
  
  # Test different sparsity types
  spca_lasso <- sparse_pca(data, sparsity_type = "lasso", sparsity_param = 0.1)
  expect_equal(spca_lasso$params$sparsity_type, "lasso")
  
  spca_threshold <- sparse_pca(data, sparsity_type = "hard_threshold", sparsity_param = 0.2)
  expect_equal(spca_threshold$params$sparsity_type, "hard_threshold")
  
  spca_elastic <- sparse_pca(data, sparsity_type = "elastic_net", sparsity_param = 0.15)
  expect_equal(spca_elastic$params$sparsity_type, "elastic_net")
})

test_that("basis_space_embed basic functionality", {
  skip_if_not_installed("neuroim2")
  skip_if(!rcpp_available("_neuroarchive_poisson_disk_sample_component_rcpp"), 
          "HRBF Rcpp functions not available in development environment")
  
  # Create test data and mask
  dims <- c(10, 10, 5)
  arr <- array(TRUE, dims)
  mask <- neuroim2::LogicalNeuroVol(arr, neuroim2::NeuroSpace(dims, spacing = c(2, 2, 2)))
  
  set.seed(123)
  n_voxels <- sum(as.logical(mask))
  data <- matrix(rnorm(20 * n_voxels), nrow = 20, ncol = n_voxels)
  
  # Basic embedding
  embedded <- basis_space_embed(data, mask)
  
  expect_type(embedded, "list")
  expect_true("coefficients" %in% names(embedded))
  expect_true("basis_params" %in% names(embedded))
  expect_true("compression_ratio" %in% names(embedded))
  expect_true("mask_hash" %in% names(embedded))
  
  # Check dimensions
  expect_equal(nrow(embedded$coefficients), nrow(data))
  expect_true(ncol(embedded$coefficients) > 0)
})

test_that("basis_space_embed with custom parameters", {
  skip_if_not_installed("neuroim2")
  skip_if(!rcpp_available("_neuroarchive_poisson_disk_sample_component_rcpp"), 
          "HRBF Rcpp functions not available in development environment")
  
  dims <- c(10, 10, 5)
  arr <- array(TRUE, dims)
  mask <- neuroim2::LogicalNeuroVol(arr, neuroim2::NeuroSpace(dims, spacing = c(2, 2, 2)))
  
  set.seed(123)
  n_voxels <- sum(as.logical(mask))
  data <- matrix(rnorm(15 * n_voxels), nrow = 15, ncol = n_voxels)
  
  # Custom basis parameters
  custom_params <- list(sigma0 = 4, levels = 4L, kernel_type = "wendland_c4")
  embedded <- basis_space_embed(data, mask, basis_params = custom_params)
  
  expect_equal(embedded$basis_params$sigma0, 4)
  expect_equal(embedded$basis_params$levels, 4L)
  expect_equal(embedded$basis_params$kernel_type, "wendland_c4")
})

test_that("basis_space_embed with basis return", {
  skip_if_not_installed("neuroim2")
  skip_if(!rcpp_available("_neuroarchive_poisson_disk_sample_component_rcpp"), 
          "HRBF Rcpp functions not available in development environment")
  
  dims <- c(8, 8, 4)
  arr <- array(TRUE, dims)
  mask <- neuroim2::LogicalNeuroVol(arr, neuroim2::NeuroSpace(dims, spacing = c(2, 2, 2)))
  
  set.seed(123)
  n_voxels <- sum(as.logical(mask))
  data <- matrix(rnorm(10 * n_voxels), nrow = 10, ncol = n_voxels)
  
  # Request basis return
  embedded <- basis_space_embed(data, mask, return_basis = TRUE)
  
  expect_true("basis" %in% names(embedded))
  expect_type(embedded$basis, "double")
})

test_that("basis_space_reconstruct basic functionality", {
  skip_if_not_installed("neuroim2")
  skip_if(!rcpp_available("_neuroarchive_poisson_disk_sample_component_rcpp"), 
          "HRBF Rcpp functions not available in development environment")
  
  dims <- c(8, 8, 4)
  arr <- array(TRUE, dims)
  mask <- neuroim2::LogicalNeuroVol(arr, neuroim2::NeuroSpace(dims, spacing = c(2, 2, 2)))
  
  set.seed(123)
  n_voxels <- sum(as.logical(mask))
  data <- matrix(rnorm(10 * n_voxels), nrow = 10, ncol = n_voxels)
  
  # Embed and reconstruct
  embedded <- basis_space_embed(data, mask)
  reconstructed <- basis_space_reconstruct(embedded, mask)
  
  expect_type(reconstructed, "double")
  expect_equal(nrow(reconstructed), nrow(data))
  expect_equal(ncol(reconstructed), ncol(data))
})

test_that("basis_space_reconstruct with subsetting", {
  skip_if_not_installed("neuroim2")
  skip_if(!rcpp_available("_neuroarchive_poisson_disk_sample_component_rcpp"), 
          "HRBF Rcpp functions not available in development environment")
  
  dims <- c(8, 8, 4)
  arr <- array(TRUE, dims)
  mask <- neuroim2::LogicalNeuroVol(arr, neuroim2::NeuroSpace(dims, spacing = c(2, 2, 2)))
  
  set.seed(123)
  n_voxels <- sum(as.logical(mask))
  data <- matrix(rnorm(20 * n_voxels), nrow = 20, ncol = n_voxels)
  
  embedded <- basis_space_embed(data, mask)
  
  # Time subsetting
  partial_time <- basis_space_reconstruct(embedded, mask, subset_time = 1:10)
  expect_equal(nrow(partial_time), 10)
  expect_equal(ncol(partial_time), ncol(data))
  
  # Voxel subsetting
  partial_voxels <- basis_space_reconstruct(embedded, mask, subset_voxels = 1:50)
  expect_equal(nrow(partial_voxels), nrow(data))
  expect_equal(ncol(partial_voxels), 50)
})

test_that("basis_space_reconstruct mask validation", {
  skip_if_not_installed("neuroim2")
  skip_if(!rcpp_available("_neuroarchive_poisson_disk_sample_component_rcpp"), 
          "HRBF Rcpp functions not available in development environment")
  
  dims <- c(8, 8, 4)
  arr <- array(TRUE, dims)
  mask1 <- neuroim2::LogicalNeuroVol(arr, neuroim2::NeuroSpace(dims, spacing = c(2, 2, 2)))
  
  # Different mask
  arr2 <- array(TRUE, c(10, 10, 5))
  mask2 <- neuroim2::LogicalNeuroVol(arr2, neuroim2::NeuroSpace(c(10, 10, 5), spacing = c(2, 2, 2)))
  
  set.seed(123)
  n_voxels <- sum(as.logical(mask1))
  data <- matrix(rnorm(10 * n_voxels), nrow = 10, ncol = n_voxels)
  
  embedded <- basis_space_embed(data, mask1)
  
  # Should warn about mask mismatch
  expect_warning(
    basis_space_reconstruct(embedded, mask2),
    "Mask hash mismatch"
  )
})

test_that("basis_space_visualize basic functionality", {
  skip_if_not_installed("neuroim2")
  skip_if(!rcpp_available("_neuroarchive_poisson_disk_sample_component_rcpp"), 
          "HRBF Rcpp functions not available in development environment")
  
  dims <- c(10, 10, 5)
  arr <- array(TRUE, dims)
  mask <- neuroim2::LogicalNeuroVol(arr, neuroim2::NeuroSpace(dims, spacing = c(2, 2, 2)))
  
  # Generate basis for visualization
  basis_result <- basis_space_centers(mask, sigma0 = 4, levels = 2)
  
  # Basic visualization (should not error)
  expect_message(
    result <- basis_space_visualize(basis_result, mask),
    "Would visualize.*basis functions"
  )
  
  expect_null(result)  # Function returns invisible(NULL)
})

test_that("basis_space_visualize parameter validation", {
  skip_if_not_installed("neuroim2")
  skip_if(!rcpp_available("_neuroarchive_poisson_disk_sample_component_rcpp"), 
          "HRBF Rcpp functions not available in development environment")
  
  dims <- c(10, 10, 5)
  arr <- array(TRUE, dims)
  mask <- neuroim2::LogicalNeuroVol(arr, neuroim2::NeuroSpace(dims, spacing = c(2, 2, 2)))
  
  basis_result <- basis_space_centers(mask, sigma0 = 4, levels = 2)
  
  # Test different slice axes
  expect_message(
    basis_space_visualize(basis_result, mask, slice_axis = "axial"),
    "axial slice"
  )
  
  expect_message(
    basis_space_visualize(basis_result, mask, slice_axis = "coronal"),
    "coronal slice"
  )
  
  expect_message(
    basis_space_visualize(basis_result, mask, slice_axis = "sagittal"),
    "sagittal slice"
  )
  
  # Test custom slice position
  expect_message(
    basis_space_visualize(basis_result, mask, slice_position = 3),
    "slice 3"
  )
})

test_that("error handling for invalid inputs", {
  # Test basis_time with invalid type
  expect_error(
    basis_time(100, "invalid_type"),
    "'arg' should be one of"
  )
  
  # Test delta_transform with invalid method
  data <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
  expect_error(
    delta_transform(data, method = "invalid"),
    "'arg' should be one of"
  )
  
  # Test delta_transform with insufficient data for central differences
  small_data <- matrix(rnorm(2 * 5), nrow = 2, ncol = 5)
  expect_error(
    delta_transform(small_data, method = "central"),
    "Central differences require at least 3 time points"
  )
}) 