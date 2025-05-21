library(testthat)
library(neuroarchive)
library(withr)


test_that("default_params for myorg.sparsepca loads schema", {
  neuroarchive:::default_param_cache_clear()
  p <- neuroarchive:::default_params("myorg.sparsepca")
  expect_equal(p$k, 50)
  expect_equal(p$alpha, 0.001)
  expect_identical(p$whiten, FALSE)
  expect_equal(p$storage_order, "component_x_voxel")
})


test_that("sparsepca forward and inverse roundtrip", {
  set.seed(1)
  X <- matrix(rnorm(40), nrow = 10, ncol = 4)
  tmp <- local_tempfile(fileext = ".h5")

  write_lna(X, file = tmp, transforms = "myorg.sparsepca",
            transform_params = list(`myorg.sparsepca` = list(k = 4)))
  h <- read_lna(tmp)
  out <- h$stash$input
  expect_equal(dim(out), dim(X))
  expect_equal(out, X, tolerance = 1e-6)
})
