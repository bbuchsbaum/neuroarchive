library(testthat)
library(neuroarchive)


test_that("default_params for basis loads schema", {
  cache_env <- get(".default_param_cache", envir = asNamespace("neuroarchive"))
  rm(list = ls(envir = cache_env), envir = cache_env)
  p <- neuroarchive:::default_params("basis")
  expect_equal(p$method, "pca")
  expect_true(is.numeric(p$k))
  expect_true(p$center)
  expect_false(p$scale)
})
