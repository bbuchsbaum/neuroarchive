library(testthat)
#library(neuroarchive)


test_that("default_params for embed loads schema", {
  neuroarchive:::default_param_cache_clear()
  p <- neuroarchive:::default_params("embed")
  expect_equal(p$basis_path, "")
  expect_null(p$center_data_with)
  expect_null(p$scale_data_with)
})

test_that("embed transform errors when basis_path missing", {
  X <- matrix(rnorm(10), nrow = 5)
  expect_error(
    core_write(X, transforms = "embed"),
    class = "lna_error_validation",
    regexp = "basis_path"
  )
})


test_that("embed transform forward computes coefficients", {
  set.seed(1)
  X <- matrix(rnorm(20), nrow = 5, ncol = 4)

  pr <- prcomp(X, rank. = 3, center = TRUE, scale. = FALSE)
  basis_mat <- t(pr$rotation[, seq_len(3), drop = FALSE])
  center_vec <- pr$center

  plan <- Plan$new()
  plan$add_payload("/basis/mat", basis_mat)
  plan$add_payload("/basis/center", center_vec)

  desc <- list(
    type = "embed",
    params = list(
      basis_path = "/basis/mat",
      center_data_with = "/basis/center"
    )
  )
  handle <- DataHandle$new(initial_stash = list(input = X), plan = plan)

  res_handle <- forward_step.embed("embed", desc, handle)
  plan <- res_handle$plan

  
  # Create a plan and handle directly
  plan <- Plan$new()
  
  # Create a simple basis matrix and add it to the plan
  basis_mat <- matrix(rnorm(8), nrow = 4, ncol = 2)  # 4×2 basis matrix
  basis_path <- "/basis/test/matrix"
  plan$add_payload(basis_path, basis_mat)
  
  # Create a handle with our matrix
  handle <- DataHandle$new(initial_stash = list(input = X), plan = plan)
  
  # Create a descriptor for embed transform
  desc <- list(
    type = "embed",
    params = list(basis_path = basis_path),
    inputs = c("input")
  )
  
  # Call forward_step.embed directly
  handle <- forward_step.embed("embed", desc, handle)
  
  # Now check the plan for coefficients
  coeff_idx <- which(plan$datasets$role == "coefficients")
  expect_length(coeff_idx, 1)
  coeff_path <- plan$datasets$path[[coeff_idx]]
  expect_true(coeff_path %in% names(plan$payloads))
  coeff <- plan$payloads[[coeff_path]]
  expect_equal(nrow(coeff), nrow(X))
})

test_that("embed transform requires numeric input", {
  plan <- Plan$new()

  basis_mat <- diag(2)
  plan$add_payload("/basis/mat", basis_mat)
  desc <- list(type = "embed", params = list(basis_path = "/basis/mat"))
  handle <- DataHandle$new(initial_stash = list(input = matrix("a", nrow = 2, ncol = 2)),
                           plan = plan)

  expect_error(
    forward_step.embed("embed", desc, handle),
    class = "lna_error_validation",
    regexp = "numeric input"
    )

  h <- DataHandle$new(initial_stash = list(input = matrix("a", nrow = 2)),
                      plan = plan)
  desc <- list(type = "embed", params = list(basis_path = "/basis/mat"),
               inputs = c("input"))
  expect_error(
    neuroarchive:::forward_step.embed("embed", desc, h),
    class = "lna_error_validation",
    regexp = "numeric"

  )
})
