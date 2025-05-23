library(testthat)

opts_env <- get(".lna_opts", envir = neuroarchive:::lna_options_env)

clear_opts <- function() {
  rm(list = ls(envir = opts_env), envir = opts_env)
  neuroarchive:::default_param_cache_clear()
}

test_that("delta() merges parameters and executes", {
  clear_opts()
  lna_options(delta = list(order = 2L))
  arr <- array(1, dim = c(1,1,1))
  pipe <- delta(arr, order = 3L)
  step <- pipe$steps[[1]]
  expect_equal(step$type, "delta")
  expect_equal(step$params$order, 3L)

  captured <- list()
  local_mocked_bindings(
    write_lna = function(x, file, transforms, transform_params, run_id) {
      captured$transforms <<- transforms
      captured$transform_params <<- transform_params
      list(ok = TRUE)
    },
    .env = asNamespace("neuroarchive")
  )
  lna_write(pipe, file = "tmp.h5")
  expect_equal(captured$transforms, "delta")
  expect_equal(captured$transform_params$delta$order, 3L)
})

test_that("temporal() verb adds step and executes", {
  clear_opts()
  lna_options(temporal = list(kind = "dct"))
  arr <- array(1, dim = c(1,1,1))
  pipe <- temporal(arr, kind = "bspline")
  step <- pipe$steps[[1]]
  expect_equal(step$type, "temporal")
  expect_equal(step$params$kind, "bspline")

  captured <- list()
  local_mocked_bindings(
    write_lna = function(x, file, transforms, transform_params, run_id) {
      captured$transforms <<- transforms
      captured$transform_params <<- transform_params
      list(ok = TRUE)
    },
    .env = asNamespace("neuroarchive")
  )
  lna_write(pipe, file = "tmp.h5")
  expect_equal(captured$transforms, "temporal")
  expect_equal(captured$transform_params$temporal$kind, "bspline")
})

test_that("hrbf() verb works", {
  clear_opts()
  arr <- array(1, dim = c(1,1,1))
  pipe <- hrbf(arr, levels = 3)
  step <- pipe$steps[[1]]
  expect_equal(step$type, "spat.hrbf")
  expect_equal(step$params$levels, 3)
})

test_that("embed() infers path after hrbf", {
  clear_opts()
  arr <- array(1, dim = c(1,1,1))
  pipe <- hrbf(arr)
  pipe <- embed(pipe)
  step <- pipe$steps[[2]]
  expect_equal(step$type, "embed.hrbf_analytic")
  expect_true(grepl("/basis/00_spat.hrbf/matrix", step$params$basis_path))
})

test_that("register_lna_verb slugging and collision", {
  reg <- get(".verb_registry", envir = neuroarchive:::lna_verb_registry_env)
  rm(list = ls(envir = reg), envir = reg)

  res <- register_lna_verb(lna_transform_type = "my.org.filt")
  expect_equal(res$name, "my_org_filt")
  expect_true(exists("my_org_filt", envir = reg))

  expect_warning(register_lna_verb("my_org_filt", "other"))
})

test_that("template registration and application with overrides", {
  reg <- get(".template_registry", envir = neuroarchive:::lna_template_registry_env)
  rm(list = ls(envir = reg), envir = reg)

  simple_template <- function(pipe, bits = 4) {
    quant(pipe, bits = bits)
  }
  register_lna_template("simple", simple_template)

  pipe <- as_pipeline(array(1))
  pipe <- apply_template(pipe, "simple", bits = 6)
  step <- pipe$steps[[1]]
  expect_equal(step$type, "quant")
  expect_equal(step$params$bits, 6)

  pipe2 <- as_pipeline(array(1))
  pipe2 <- apply_template(pipe2, "simple", quant.bits = 7)
  step2 <- pipe2$steps[[1]]
  expect_equal(step2$params$bits, 7)
})
