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
