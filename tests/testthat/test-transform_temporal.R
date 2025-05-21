library(testthat)
library(neuroarchive)


test_that("default_params for temporal loads schema", {
  cache_env <- get(".default_param_cache", envir = asNamespace("neuroarchive"))
  rm(list = ls(envir = cache_env), envir = cache_env)
  p <- neuroarchive:::default_params("temporal")
  expect_equal(p$kind, "dct")
  expect_equal(p$scope, "global")
  expect_true(is.numeric(p$n_basis))
})
