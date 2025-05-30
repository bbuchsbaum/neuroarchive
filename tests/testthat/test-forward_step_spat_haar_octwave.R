library(testthat)

# Unit tests for forward_step.spat.haar_octwave (HWT-S1-10)

test_that("forward_step.spat.haar_octwave basic IO", {
  mask <- array(TRUE, dim = c(2,2,2))
  X <- array(1:16, dim = c(2,2,2,2))
  plan <- Plan$new()
  h <- DataHandle$new(initial_stash = list(input = X),
                      plan = plan,
                      mask_info = list(mask = mask, active_voxels = 8),
                      run_ids = "run-01", current_run_id = "run-01")
  desc <- list(type = "spat.haar_octwave", params = list(levels = 1))

  called <- FALSE
  testthat::local_mocked_bindings(
    perform_haar_lift_analysis = function(data_matrix_T_x_Nmask, mask_3d_array, levels, z_order_seed = 42L) {
      called <<- TRUE
      list(
        root = matrix(1, nrow = 2, ncol = 1),
        detail = list(matrix(2, nrow = 2, ncol = 8))
      )
    },
    .env = asNamespace("neuroarchive")
  )

  out <- neuroarchive:::forward_step.spat.haar_octwave("spat.haar_octwave", desc, h)

  expect_true(called)

  dname <- names(out$plan$descriptors)[1]
  stored_desc <- out$plan$descriptors[[dname]]
  p <- stored_desc$params

  expect_equal(p$num_voxels_in_mask, 8)
  expect_equal(p$octree_bounding_box_mask_space, c(0,1,0,1,0,1))
  expect_match(p$morton_hash_mask_indices, "^sha1:[a-f0-9]{40}$")
  expect_equal(p$num_coeffs_per_level$lowpass, c(1L, 1L))
  expect_equal(p$num_coeffs_per_level$detail, c(8L))

  expect_true("/wavelet/level_ROOT/coefficients" %in% out$plan$datasets$path)
  expect_true("/wavelet/level_0/detail_coefficients" %in% out$plan$datasets$path)

  root_mat <- out$plan$payloads[["/wavelet/level_ROOT/coefficients"]]
  detail_mat <- out$plan$payloads[["/wavelet/level_0/detail_coefficients"]]
  expect_equal(dim(root_mat), c(2,1))
  expect_equal(dim(detail_mat), c(2,8))

  expect_true(out$has_key("wavelet_coefficients"))
  stash_mat <- out$stash$wavelet_coefficients
  expect_equal(dim(stash_mat), c(2,9))
})

test_that("forward_step.spat.haar_octwave sparsifies details", {
  mask <- array(TRUE, dim = c(2,2,2))
  X <- array(1:16, dim = c(2,2,2,2))
  plan <- Plan$new()
  h <- DataHandle$new(initial_stash = list(input = X),
                      plan = plan,
                      mask_info = list(mask = mask, active_voxels = 8),
                      run_ids = "run-01", current_run_id = "run-01")
  desc <- list(type = "spat.haar_octwave",
               params = list(levels = 1,
                              detail_threshold_type = "absolute",
                              detail_threshold_value = 0.5))

  detail_in <- matrix(c(0.1, 0.6, -0.2, 1.0, 0.4, -0.7, 2.0, 0.3,
                         0.4, -0.8, 0.1, 1.2, -0.1, 0.2, 2.5, -0.4),
                       nrow = 2, byrow = TRUE)

  testthat::local_mocked_bindings(
    perform_haar_lift_analysis = function(data_matrix_T_x_Nmask, mask_3d_array,
                                          levels, z_order_seed = 42L) {
      list(
        root = matrix(1, nrow = 2, ncol = 1),
        detail = list(detail_in)
      )
    },
    .env = asNamespace("neuroarchive")
  )

  out <- neuroarchive:::forward_step.spat.haar_octwave("spat.haar_octwave", desc, h)

  dmat <- out$plan$payloads[["/wavelet/level_0/detail_coefficients"]]
  expect_true(all(abs(dmat[abs(detail_in) < 0.5]) == 0))
  expect_true(all(dmat[abs(detail_in) >= 0.5] == detail_in[abs(detail_in) >= 0.5]))
})

