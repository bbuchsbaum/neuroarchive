library(testthat)
#library(neuroarchive)


test_that("default_params for basis loads schema", {
  neuroarchive:::default_param_cache_clear()
  p <- neuroarchive:::default_params("basis")
  expect_equal(p$method, "pca")
  expect_true(is.numeric(p$k))
  expect_true(p$center)
  expect_false(p$scale)
})


test_that("forward_step.basis validates storage_order", {
  plan <- Plan$new()
  h <- DataHandle$new(initial_stash = list(input = matrix(1:4, nrow = 2)),
                      plan = plan)
  desc <- list(type = "basis",
               params = list(storage_order = "invalid"),
               inputs = c("input"))

  expect_error(
    neuroarchive:::forward_step.basis("basis", desc, h),
    class = "lna_error_validation",
    regexp = "Invalid storage_order"
  )
})


test_that("forward_step.basis truncates k when PCA returns fewer components", {
  plan <- Plan$new()
  X_small <- matrix(rnorm(10), nrow = 2, ncol = 5) # 2x5 matrix
  # Store the original handle to access its modified plan later
  original_handle <- DataHandle$new(initial_stash = list(input = X_small), plan = plan)
  desc <- list(type = "basis", params = list(k = 5), inputs = c("input")) # Request k=5

  # Determine expected k and warning based on whether irlba is available
  expected_final_k <- NA_integer_
  expected_warning_msg_regex <- ""
  
  min_dim_X_small <- min(dim(X_small))

  if (requireNamespace("irlba", quietly = TRUE)) {
    # irlba: k must be < min_dim. Max k is min_dim - 1.
    expected_final_k <- max(1, min_dim_X_small - 1)
    expected_warning_msg_regex <- sprintf(
      "Requested k=5 but irlba::prcomp_irlba can only compute %d components .* truncating k to %d",
      expected_final_k, expected_final_k
    )
  } else {
    # stats::prcomp: k can be <= min_dim. Max k is min_dim.
    expected_final_k <- max(1, min_dim_X_small)
    expected_warning_msg_regex <- sprintf(
      "Requested k=5 but stats::prcomp can only compute %d components .* truncating k to %d",
      expected_final_k, expected_final_k
    )
  }

  h_after_forward_step <- NULL
  expect_warning(
    h_after_forward_step <- neuroarchive:::forward_step.basis("basis", desc, original_handle),
    regexp = expected_warning_msg_regex
  )

  # After expect_warning, original_handle's internal plan object will have been modified
  # by forward_step.basis. h_after_forward_step should be original_handle.
  # Let's verify h_after_forward_step is indeed the handle that was modified.
  # If neuroarchive:::forward_step.basis returns the modified handle, this is fine.

  testthat::expect_true(!is.null(h_after_forward_step$plan), "Plan object in returned handle should not be null")
  current_plan_datasets <- h_after_forward_step$plan$datasets
  testthat::expect_true(tibble::is_tibble(current_plan_datasets), "Plan datasets should be a tibble")
  testthat::expect_gt(nrow(current_plan_datasets), 0, "Plan datasets should not be empty")

  basis_matrix_def_row <- current_plan_datasets[current_plan_datasets$role == "basis_matrix", ]

  testthat::expect_equal(nrow(basis_matrix_def_row), 1,
                         info = "Should be exactly one 'basis_matrix' role definition in the plan.")

  params_json_str <- basis_matrix_def_row$params_json[1]
  params <- jsonlite::fromJSON(params_json_str)

  expect_equal(params$k, expected_final_k)

  payload_path <- basis_matrix_def_row$path[1]
  payload <- h_after_forward_step$plan$payloads[[payload_path]]

  # Default storage_order is "component_x_voxel", so basis matrix is k x n_voxels
  expect_equal(nrow(payload), expected_final_k) # Number of components
  expect_equal(ncol(payload), ncol(X_small))   # Number of voxels
})
