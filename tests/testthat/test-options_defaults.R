library(testthat)
library(neuroarchive)

# Helper to access internal env
opts_env <- get(".lna_opts", envir = neuroarchive:::lna_options_env)

# Ensure a clean state
teardown({
  rm(list = ls(envir = opts_env), envir = opts_env)
  neuroarchive:::default_param_cache_clear()
  schema_env <- get(".schema_cache", envir = asNamespace("neuroarchive"))
  rm(list = ls(envir = schema_env), envir = schema_env)
})

# Test lna_options set/get

test_that("lna_options set and get work", {
  lna_options(write.compression_level = 3)
  expect_equal(lna_options("write.compression_level")[[1]], 3)

  lna_options(write.chunk_target_mib = 2)
  expect_equal(lna_options("write.chunk_target_mib")[[1]], 2)

  lna_options(foo = "bar", baz = 1)
  res <- lna_options("foo", "baz")
  expect_identical(res$foo, "bar")
  expect_identical(res$baz, 1)
})

# Test default_params caching behavior

test_that("default_params warns and caches empty list when schema missing", {
  neuroarchive:::default_param_cache_clear()
  cache_env <- get(".default_param_cache", envir = asNamespace("neuroarchive"))
  schema_env <- get(".schema_cache", envir = asNamespace("neuroarchive"))
  rm(list = ls(envir = schema_env), envir = schema_env)

  expect_warning(p1 <- neuroarchive:::default_params("foo"), "not found")
  expect_equal(p1, list())
  expect_true("foo" %in% ls(envir = cache_env))
  expect_false("foo" %in% ls(envir = schema_env))

  p2 <- neuroarchive:::default_params("foo")
  expect_identical(p1, p2)
})

test_that("default_params loads defaults from schema and caches", {
  neuroarchive:::default_param_cache_clear()
  cache_env <- get(".default_param_cache", envir = asNamespace("neuroarchive"))
  schema_env <- get(".schema_cache", envir = asNamespace("neuroarchive"))
  rm(list = ls(envir = schema_env), envir = schema_env)

  expect_false("test" %in% ls(envir = cache_env))
  d1 <- neuroarchive:::default_params("test")
  expect_equal(d1, list(a = 1L, b = "x", nested = list(c = 0.5)))
  expect_true("test" %in% ls(envir = cache_env))
  expect_true("test" %in% ls(envir = schema_env))

  d2 <- neuroarchive:::default_params("test")
  expect_identical(d1, d2)
})

test_that("defaults are extracted from array items", {
  cache_env <- get(".default_param_cache", envir = asNamespace("neuroarchive"))
  schema_env <- get(".schema_cache", envir = asNamespace("neuroarchive"))
  rm(list = ls(envir = cache_env), envir = cache_env)
  rm(list = ls(envir = schema_env), envir = schema_env)

  d <- neuroarchive:::default_params("test_array")
  expect_equal(d$numArray$items, 2)
  expect_equal(d$objArray$items, list(flag = TRUE))
})
