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
}


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
