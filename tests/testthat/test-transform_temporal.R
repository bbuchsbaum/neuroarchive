library(testthat)
library(neuroarchive)


test_that("default_params for temporal loads schema", {
  cache_env <- get(".default_param_cache", envir = asNamespace("neuroarchive"))
  rm(list = ls(envir = cache_env), envir = cache_env)
  p <- neuroarchive:::default_params("temporal")
  expect_equal(p$kind, "dct")
  expect_true(is.numeric(p$n_basis))
})


test_that("temporal bspline forward/inverse roundtrip", {
  X <- matrix(sin(seq(0, pi, length.out = 10)), nrow = 10, ncol = 1)
  tmp <- local_tempfile(fileext = ".h5")

  write_lna(X, file = tmp, transforms = "temporal",
            transform_params = list(temporal = list(kind = "bspline",
                                                    n_basis = 6, order = 3)))
  h <- read_lna(tmp)
  out <- h$stash$input
  expect_equal(dim(out), dim(X))
  expect_true(mean(abs(out - X)) < 1)
})
