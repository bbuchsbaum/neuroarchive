library(testthat)
library(neuroarchive)

# Helper to access internal env
opts_env <- get(".lna_opts", envir = neuroarchive:::lna_options_env)

# Ensure a clean state
teardown({
  rm(list = ls(envir = opts_env), envir = opts_env)
  cache_env <- get(".default_param_cache", envir = asNamespace("neuroarchive"))
  rm(list = ls(envir = cache_env), envir = cache_env)
})

# Test lna_options set/get

test_that("lna_options set and get work", {
  lna_options(write.compression = 3)
  expect_equal(lna_options("write.compression")[[1]], 3)

  lna_options(foo = "bar", baz = 1)
  res <- lna_options("foo", "baz")
  expect_identical(res$foo, "bar")
  expect_identical(res$baz, 1)
})

# Test default_params caching behavior

test_that("default_params returns empty list and caches", {
  cache_env <- get(".default_param_cache", envir = asNamespace("neuroarchive"))
  rm(list = ls(envir = cache_env), envir = cache_env)

  p1 <- neuroarchive:::default_params("foo")
  expect_equal(p1, list())
  expect_true("foo" %in% ls(envir = cache_env))

  p2 <- neuroarchive:::default_params("foo")
  expect_identical(p1, p2)
})
