library(testthat)
#library(neuroarchive)

# Helper to access internal options environment
opts_env <- get(".lna_opts", envir = neuroarchive:::lna_options_env)

# Ensure a clean state after each test
teardown({
  rm(list = ls(envir = opts_env), envir = opts_env)
  cache_env <- get(".default_param_cache", envir = asNamespace("neuroarchive"))
  rm(list = ls(envir = cache_env), envir = cache_env)
  schema_env <- get(".schema_cache", envir = asNamespace("neuroarchive"))
  rm(list = ls(envir = schema_env), envir = schema_env)
})

# Verify merge order defaults -> options -> user

test_that("resolve_transform_params merges in correct precedence", {
  cache_env <- get(".default_param_cache", envir = asNamespace("neuroarchive"))
  rm(list = ls(envir = cache_env), envir = cache_env)
  schema_env <- get(".schema_cache", envir = asNamespace("neuroarchive"))
  rm(list = ls(envir = schema_env), envir = schema_env)

  lna_options(test = list(a = 10L, nested = list(d = 5)))
  user <- list(test = list(a = 15L, nested = list(e = 9)))

  res <- neuroarchive:::resolve_transform_params("test", user)$test

  expect_equal(res$a, 15L)
  expect_equal(res$b, "x")
  expect_equal(res$nested$c, 0.5)
  expect_equal(res$nested$d, 5)
  expect_equal(res$nested$e, 9)
})

# Explicit NULL values are preserved with keep.null = TRUE

test_that("resolve_transform_params keeps explicit NULL values", {
  cache_env <- get(".default_param_cache", envir = asNamespace("neuroarchive"))
  rm(list = ls(envir = cache_env), envir = cache_env)
  schema_env <- get(".schema_cache", envir = asNamespace("neuroarchive"))
  rm(list = ls(envir = schema_env), envir = schema_env)

  lna_options(test = list(nested = list(c = NULL)))
  user <- list(test = list(b = NULL))

  res <- neuroarchive:::resolve_transform_params("test", user)$test

  expect_true("c" %in% names(res$nested))
  expect_null(res$nested$c)
  expect_true("b" %in% names(res))
  expect_null(res$b)
  expect_equal(res$a, 1L)
})

