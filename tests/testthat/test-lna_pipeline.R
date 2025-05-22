library(testthat)

# Tests for lna_pipeline$set_input

test_that("set_input handles single array", {
  pipe <- neuroarchive:::lna_pipeline$new()
  arr <- array(1:24, dim = c(2,3,4))
  pipe$set_input(arr)
  expect_identical(pipe$input, arr)
  expect_equal(pipe$runs, "run-01")
  expect_match(pipe$input_summary, "1 run")
  expect_match(pipe$input_summary, "4 TR")
  expect_null(pipe$engine_opts$chunk_mb_suggestion)
})

test_that("set_input handles named list", {
  pipe <- neuroarchive:::lna_pipeline$new()
  arr1 <- array(1:8, dim = c(2,2,2))
  arr2 <- array(1:8, dim = c(2,2,2))
  lst <- list(foo = arr1, bar = arr2)
  pipe$set_input(lst)
  expect_equal(pipe$runs, c("foo", "bar"))
  expect_identical(pipe$input, lst)
})

test_that("set_input handles unnamed list with run_ids", {
  pipe <- neuroarchive:::lna_pipeline$new()
  arr1 <- array(1:8, dim = c(2,2,2))
  arr2 <- array(1:8, dim = c(2,2,2))
  lst <- list(arr1, arr2)
  pipe$set_input(lst, run_ids = c("r1", "r2"))
  expect_equal(pipe$runs, c("r1", "r2"))
})

# Test lna_pipeline$new defaults

test_that("lna_pipeline fields initialise correctly", {
  pipe <- neuroarchive:::lna_pipeline$new()
  expect_null(pipe$input)
  expect_identical(pipe$input_summary, "")
  expect_equal(pipe$runs, character())
  expect_equal(pipe$steps, list())
  expect_equal(pipe$engine_opts, list())
})

# Test add_step behaviour

test_that("add_step appends step specification", {
  pipe <- neuroarchive:::lna_pipeline$new()
  step <- list(type = "quant", params = list(bits = 8))
  pipe$add_step(step)
  expect_equal(length(pipe$steps), 1L)
  expect_identical(pipe$steps[[1]], step)
})

# Tests for as_pipeline

test_that("as_pipeline creates pipeline from array", {
  arr <- array(1:8, dim = c(2,2,2))
  pipe <- as_pipeline(arr)
  expect_s3_class(pipe, "lna_pipeline")
  expect_identical(pipe$input, arr)
  expect_equal(pipe$runs, "run-01")
})

test_that("as_pipeline handles list input", {
  lst <- list(a = array(1, dim = c(1,1,1)), b = array(2, dim = c(1,1,1)))
  pipe <- as_pipeline(lst)
  expect_s3_class(pipe, "lna_pipeline")
  expect_equal(pipe$runs, c("a", "b"))
})

# Tests for lna_write argument forwarding

test_that("lna_write forwards arguments to write_lna", {
  arr <- array(1, dim = c(1,1,1))
  pipe <- as_pipeline(arr)
  pipe$add_step(list(type = "quant", params = list(bits = 8)))

  captured <- list()
  fake_result <- list(fake = TRUE)

  local_mocked_bindings(
    write_lna = function(x, file, transforms, transform_params, run_id, header) {
      captured$x <<- x
      captured$file <<- file
      captured$transforms <<- transforms
      captured$transform_params <<- transform_params
      captured$run_id <<- run_id
      captured$header <<- header
      fake_result
    },
    .env = asNamespace("neuroarchive")
  )

  res <- lna_write(pipe, file = "foo.h5", header = list(a = 1))

  expect_identical(res, fake_result)
  expect_identical(captured$x, arr)
  expect_equal(captured$file, "foo.h5")
  expect_equal(captured$transforms, "quant")
  expect_equal(captured$transform_params, list(quant = list(bits = 8)))
  expect_equal(captured$run_id, "run-01")
  expect_equal(captured$header, list(a = 1))
})

# Empty pipeline still calls write_lna with no transforms

test_that("lna_write works with empty pipeline", {
  arr <- array(1, dim = c(1,1,1))
  pipe <- as_pipeline(arr)

  captured <- list()
  local_mocked_bindings(
    write_lna = function(x, file, transforms, transform_params, run_id) {
      captured$transforms <<- transforms
      captured$transform_params <<- transform_params
      list(ok = TRUE)
    },
    .env = asNamespace("neuroarchive")
  )

  lna_write(pipe, file = "bar.h5")
  expect_length(captured$transforms, 0L)
  expect_equal(captured$transform_params, list())
})

# Error surfacing from write_lna

test_that("lna_write surfaces core errors with context", {
  arr <- array(1, dim = c(1,1,1))
  pipe <- as_pipeline(arr)
  pipe$add_step(list(type = "quant", params = list()))

  local_mocked_bindings(
    write_lna = function(...) {
      neuroarchive:::abort_lna(
        "boom",
        step_index = 0,
        transform_type = "quant",
        .subclass = "lna_error_internal"
      )
    },
    .env = asNamespace("neuroarchive")
  )

  expect_error(
    lna_write(pipe, file = "out.h5"),
    class = "lna_error_internal",
    regexp = "Pipeline failure in step 1 \(type='quant'\)"
  )
})

