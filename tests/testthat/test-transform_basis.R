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


test_that("forward_step.basis truncates k when PCA returns fewer components", {
  plan <- Plan$new()
  X <- matrix(rnorm(10), nrow = 2)
  handle <- DataHandle$new(initial_stash = list(input = X), plan = plan)
  desc <- list(type = "basis", params = list(k = 5), inputs = c("input"))
  expect_warning(
    h <- forward_step.basis("basis", desc, handle),
    "truncating"
  )
  defs <- h$plan$datasets
  params <- jsonlite::fromJSON(defs$params_json)
  expect_equal(params$k, 2)
  payload <- h$plan$payloads[[defs$payload_key]]
  expect_equal(nrow(payload), params$k)
})
