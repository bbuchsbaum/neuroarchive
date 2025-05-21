library(testthat)
library(neuroarchive)
library(withr)


test_that("temporal transform forward and inverse roundtrip", {
  set.seed(1)
  X <- matrix(rnorm(40), nrow = 10, ncol = 4)
  tmp <- local_tempfile(fileext = ".h5")

  res <- write_lna(X, file = tmp, transforms = "temporal",
                   transform_params = list(temporal = list(n_basis = 10)))
  expect_true(file.exists(tmp))

  h <- read_lna(tmp)
  out <- h$stash$input
  expect_equal(dim(out), dim(X))
  expect_equal(out, X, tolerance = 1e-6)
})


test_that("invert_step.temporal applies time_idx subset", {
  X <- matrix(seq_len(40), nrow = 10, ncol = 4)
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(X, file = tmp, transforms = "temporal",
            transform_params = list(temporal = list(n_basis = 10)))

  h <- read_lna(tmp, time_idx = c(1,5,10))
  out <- h$stash$input
  expect_equal(dim(out), c(3, ncol(X)))
  expect_equal(out, X[c(1,5,10), ])
})
