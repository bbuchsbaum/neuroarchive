library(testthat)


test_that("schema_cache_clear empties internal cache", {
  schema_cache_clear()
  cache_env <- get(".schema_cache", envir = asNamespace("lna"))
  cache_env$foo <- 1
  cache_env$bar <- list(a = 2)
  expect_gt(length(ls(envir = cache_env)), 0)
  schema_cache_clear()
  expect_equal(length(ls(envir = cache_env)), 0)
})

