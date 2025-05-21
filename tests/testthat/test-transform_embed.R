library(testthat)
library(neuroarchive)


test_that("default_params for embed loads schema", {
  cache_env <- get(".default_param_cache", envir = asNamespace("neuroarchive"))
  rm(list = ls(envir = cache_env), envir = cache_env)
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
  res <- core_write(X, transforms = c("basis", "embed"))
  plan <- res$plan
  coeff_idx <- which(plan$datasets$role == "coefficients")
  expect_length(coeff_idx, 1)
  coeff_path <- plan$datasets$path[[coeff_idx]]
  expect_true(coeff_path %in% names(plan$payloads))
  coeff <- plan$payloads[[coeff_path]]
  expect_equal(nrow(coeff), nrow(X))
})
