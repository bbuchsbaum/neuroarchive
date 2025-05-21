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
