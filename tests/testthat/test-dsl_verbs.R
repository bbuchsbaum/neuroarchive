library(testthat)

# Access internal options environment for cleanup
opts_env <- get(".lna_opts", envir = neuroarchive:::lna_options_env)

test_that("quant() merges defaults, options, and user args", {
  rm(list = ls(envir = opts_env), envir = opts_env)
  neuroarchive:::default_param_cache_clear()

  lna_options(quant = list(bits = 10L, method = "sd"))
  pipe <- quant(array(1:4), center = FALSE)
  step <- pipe$steps()[[1]]

  expect_equal(step$type, "quant")
  expect_equal(step$params$bits, 10L)
  expect_equal(step$params$method, "sd")
  expect_false(step$params$center)

  pipe2 <- quant(array(1:4), bits = 6)
  step2 <- pipe2$steps()[[1]]
  expect_equal(step2$params$bits, 6)
})

test_that("pca() merges defaults, options, and user args", {
  rm(list = ls(envir = opts_env), envir = opts_env)
  neuroarchive:::default_param_cache_clear()

  lna_options(basis = list(k = 30L))
  pipe <- pca(matrix(rnorm(20), nrow = 4), center = FALSE)
  step <- pipe$steps()[[1]]

  expect_equal(step$type, "basis")
  expect_equal(step$params$k, 30L)
  expect_equal(step$params$method, "pca")
  expect_false(step$params$center)

  pipe2 <- pca(matrix(rnorm(20), nrow = 4), k = 12)
  step2 <- pipe2$steps()[[1]]
  expect_equal(step2$params$k, 12)
})

test_that("quant() pipeline executes via lna_write", {
  rm(list = ls(envir = opts_env), envir = opts_env)
  arr <- array(1, dim = c(1,1,1))

  captured <- list()
  local_mock_write_lna(function(x, file, transforms, transform_params, run_id, checksum = "none") {
    captured$transforms <<- transforms
    captured$transform_params <<- transform_params
    list(ok = TRUE)
  })

  pipe <- quant(arr)
  lna_write(pipe, file = "foo.h5")

  expect_equal(captured$transforms, "quant")
  expect_equal(captured$transform_params$step_01$bits, 8)
})

test_that("pca -> embed -> quant pipeline executes", {
  rm(list = ls(envir = opts_env), envir = opts_env)
  X <- matrix(rnorm(20), nrow = 5)

  captured <- list()
  local_mock_write_lna(function(x, file, transforms, transform_params, run_id, checksum = "none") {
    captured$transforms <<- transforms
    captured$transform_params <<- transform_params
    list(ok = TRUE)
  })

  pipe <- pca(X, k = 2)
  pipe <- embed(pipe)
  pipe <- quant(pipe, bits = 6)
  lna_write(pipe, file = "bar.h5")

  expect_equal(captured$transforms, c("basis", "embed.pca", "quant"))
  expect_equal(captured$transform_params$step_01$k, 2)
  expect_equal(captured$transform_params$step_03$bits, 6)
  expect_equal(captured$transform_params$step_02$basis_path, "/basis/00_basis/matrix")
  expect_equal(captured$transform_params$step_02$center_data_with, "/basis/00_basis/center")
})

test_that("embed() without prior basis step errors", {
  rm(list = ls(envir = opts_env), envir = opts_env)
  expect_error(embed(as_pipeline(matrix(1, 2, 2))), "basis-producing")
})

test_that("pca -> quant errors about missing embed", {
  rm(list = ls(envir = opts_env), envir = opts_env)
  arr <- array(rnorm(20), dim = c(5, 4))

  local_mock_write_lna(function(x, file, transforms, transform_params, run_id, checksum = "none") {
    list(ok = TRUE)
  })

  pipe <- pca(arr, k = 2)
  pipe <- quant(pipe, bits = 6)

  expect_error(
    lna_write(pipe, file = "warn.h5"),
    regexp = "basis.*without a following 'embed'|allow this dense pass-through pipeline explicitly",
    class = "lna_error_validation"
  )
})

test_that("pca -> quant allows explicit dense pass-through opt-in", {
  rm(list = ls(envir = opts_env), envir = opts_env)
  arr <- array(rnorm(20), dim = c(5, 4))

  captured <- list()
  local_mock_write_lna(function(x, file, transforms, transform_params, run_id, checksum = "none") {
    captured$transform_params <<- transform_params
    list(ok = TRUE)
  })

  pipe <- pca(arr, k = 2, allow_dense_passthrough = TRUE)
  pipe <- quant(pipe, bits = 6)
  lna_write(pipe, file = "ok.h5")

  expect_null(captured$transform_params$step_01$allow_dense_passthrough)
})
