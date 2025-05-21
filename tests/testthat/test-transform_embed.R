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
