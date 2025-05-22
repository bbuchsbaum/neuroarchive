library(testthat)
#library(neuroarchive)
library(withr)


test_that("default_params for myorg.sparsepca loads schema", {
  neuroarchive:::default_param_cache_clear()
  p <- neuroarchive:::default_params("myorg.sparsepca")
  expect_equal(p$k, 50)
  expect_equal(p$alpha, 0.001)
  expect_identical(p$whiten, FALSE)
  expect_equal(p$storage_order, "component_x_voxel")
  expect_equal(p$seed, 42)
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

test_that("whitening centers and scales the matrix", {
  set.seed(1)
  X <- matrix(rnorm(60), nrow = 15)
  plan <- Plan$new()
  h <- DataHandle$new(initial_stash = list(input = X), plan = plan)
  desc <- list(type = "myorg.sparsepca", params = list(k = 3, whiten = TRUE))
  h2 <- neuroarchive:::forward_step.myorg.sparsepca("myorg.sparsepca", desc, h)
  B <- h2$stash$sparsepca_basis
  E <- h2$stash$sparsepca_embedding
  if (identical(desc$params$storage_order %||% "component_x_voxel", "component_x_voxel")) {
    B <- t(B)
  }
  Xw <- E %*% B
  expect_true(all(abs(colMeans(Xw)) < 1e-6))
  expect_true(all(abs(apply(Xw, 2, sd) - 1) < 1e-6))
})

test_that("seed parameter yields deterministic results", {
  set.seed(2)
  X <- matrix(rnorm(40), nrow = 10)
  plan1 <- Plan$new(); plan2 <- Plan$new()
  desc <- list(type = "myorg.sparsepca", params = list(k = 2, seed = 123))
  h1 <- DataHandle$new(initial_stash = list(input = X), plan = plan1)
  h2 <- DataHandle$new(initial_stash = list(input = X), plan = plan2)
  r1 <- neuroarchive:::forward_step.myorg.sparsepca("myorg.sparsepca", desc, h1)
  r2 <- neuroarchive:::forward_step.myorg.sparsepca("myorg.sparsepca", desc, h2)
  expect_equal(r1$stash$sparsepca_basis, r2$stash$sparsepca_basis)
  expect_equal(r1$stash$sparsepca_embedding, r2$stash$sparsepca_embedding)
})

test_that("singular values dataset is written", {
  set.seed(3)
  X <- matrix(rnorm(40), nrow = 10)
  tmp <- local_tempfile(fileext = ".h5")
  res <- write_lna(X, file = tmp, transforms = "myorg.sparsepca",
                   transform_params = list(`myorg.sparsepca` = list(k = 2)))
  sv_idx <- which(res$plan$datasets$role == "singular_values")
  expect_length(sv_idx, 1)
  h5 <- hdf5r::H5File$new(tmp, mode = "r")
  vals <- h5[[res$plan$datasets$path[[sv_idx]]]]$read()
  expect_length(vals, 2)
  h5$close_all()
})
