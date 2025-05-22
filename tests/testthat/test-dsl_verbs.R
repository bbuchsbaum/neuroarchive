library(testthat)

# Access internal options environment for cleanup
opts_env <- get(".lna_opts", envir = neuroarchive:::lna_options_env)

test_that("quant() merges defaults, options, and user args", {
  rm(list = ls(envir = opts_env), envir = opts_env)
  neuroarchive:::default_param_cache_clear()

  lna_options(quant = list(bits = 10L, method = "sd"))
  pipe <- quant(array(1:4), center = FALSE)
  step <- pipe$steps[[1]]

  expect_equal(step$type, "quant")
  expect_equal(step$params$bits, 10L)
  expect_equal(step$params$method, "sd")
  expect_false(step$params$center)

  pipe2 <- quant(array(1:4), bits = 6)
  step2 <- pipe2$steps[[1]]
  expect_equal(step2$params$bits, 6)
})

test_that("pca() merges defaults, options, and user args", {
  rm(list = ls(envir = opts_env), envir = opts_env)
  neuroarchive:::default_param_cache_clear()

  lna_options(basis = list(k = 30L))
  pipe <- pca(matrix(rnorm(20), nrow = 4), center = FALSE)
  step <- pipe$steps[[1]]

  expect_equal(step$type, "basis")
  expect_equal(step$params$k, 30L)
  expect_equal(step$params$method, "pca")
  expect_false(step$params$center)

  pipe2 <- pca(matrix(rnorm(20), nrow = 4), k = 12)
  step2 <- pipe2$steps[[1]]
  expect_equal(step2$params$k, 12)
})

test_that("quant() pipeline executes via lna_write", {
  rm(list = ls(envir = opts_env), envir = opts_env)
  arr <- array(1, dim = c(1,1,1))

  captured <- list()
  local_mocked_bindings(
    write_lna = function(x, file, transforms, transform_params, run_id) {
      captured$transforms <<- transforms
      captured$transform_params <<- transform_params
      list(ok = TRUE)
    },
    .env = asNamespace("neuroarchive")
  )

  pipe <- quant(arr)
  lna_write(pipe, file = "foo.h5")

  expect_equal(captured$transforms, "quant")
  expect_equal(captured$transform_params$quant$bits, 8)
})

test_that("pca -> embed -> quant pipeline executes", {
  rm(list = ls(envir = opts_env), envir = opts_env)
  X <- matrix(rnorm(20), nrow = 5)

  captured <- list()
  local_mocked_bindings(
    write_lna = function(x, file, transforms, transform_params, run_id) {
      captured$transforms <<- transforms
      captured$transform_params <<- transform_params
      list(ok = TRUE)
    },
    .env = asNamespace("neuroarchive")
  )

  pipe <- pca(X, k = 2)
  pipe$add_step(list(type = "embed", params = list(basis_path = "/basis/00_basis/matrix")))
  pipe <- quant(pipe, bits = 6)
  lna_write(pipe, file = "bar.h5")

  expect_equal(captured$transforms, c("basis", "embed", "quant"))
  expect_equal(captured$transform_params$basis$k, 2)
  expect_equal(captured$transform_params$quant$bits, 6)
  expect_equal(captured$transform_params$embed$basis_path, "/basis/00_basis/matrix")
})
