library(testthat)
library(neuroarchive)


test_that("default_params for delta loads schema", {
  cache_env <- get(".default_param_cache", envir = asNamespace("neuroarchive"))
  rm(list = ls(envir = cache_env), envir = cache_env)
  p <- neuroarchive:::default_params("delta")
  expect_equal(p$order, 1)
  expect_equal(p$axis, 4)
  expect_equal(p$reference_value_storage, "first_value_verbatim")
  expect_null(p$delta_quantization_bits)
  expect_equal(p$coding_method, "none")
})
