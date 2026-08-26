library(testthat)
  

#' Mock forward_step methods that simply add descriptors to the plan
forward_step.tA <- function(type, desc, handle) {
  handle$plan$add_descriptor(handle$plan$get_next_filename(type), desc)
  handle
}

forward_step.tB <- function(type, desc, handle) {
  handle$plan$add_descriptor(handle$plan$get_next_filename(type), desc)
  handle
}

assign("forward_step.tA", forward_step.tA, envir = .GlobalEnv)
assign("forward_step.tB", forward_step.tB, envir = .GlobalEnv)
withr::defer({
  rm(forward_step.tA, envir = .GlobalEnv)
  rm(forward_step.tB, envir = .GlobalEnv)
})


# Core test

test_that("core_write executes forward loop and merges params", {
  result <- core_write(x = array(1, dim = c(1, 1, 1)), transforms = c("tA", "tB"), transform_params = list(tB = list(foo = "bar")))
  plan <- result$plan
  expect_equal(plan$next_index, 2L)
  expect_equal(length(plan$descriptors), 2)
  expect_equal(plan$descriptors[[1]]$type, "tA")
  expect_equal(plan$descriptors[[2]]$params, list(foo = "bar"))
})

test_that("transform_params merging honors precedence and deep merge", {
  opts_env <- get(".lna_opts", envir = neuroarchive:::lna_options_env)
  rm(list = ls(envir = opts_env), envir = opts_env)
  lna_options(tB = list(a = 10, nested = list(x = 1)))

  res <- core_write(
    x = array(1, dim = c(1, 1, 1)),
    transforms = c("tB"),
    transform_params = list(tB = list(b = 20, nested = list(y = 5)))
  )

  expect_equal(
    res$plan$descriptors[[1]]$params,
    list(a = 10, nested = list(x = 1, y = 5), b = 20)
  )
})

test_that("unknown transform names in transform_params error", {
  expect_error(
    core_write(x = array(1, dim = c(1, 1, 1)), transforms = c("tA"), transform_params = list(tB = list())),
    class = "lna_error_validation"
  )
})

test_that("unnamed list input is dispatched to one forward pipeline per run", {
  captured <- list(run_ids = character(), current_run = character(), values = numeric())
  forward_step.runTest <- function(type, desc, handle) {
    captured$run_ids <<- c(captured$run_ids, handle$run_ids)
    captured$current_run <<- c(captured$current_run, handle$current_run_id)
    captured$values <<- c(captured$values, unique(as.numeric(handle$stash$input)))
    handle
  }
  assign("forward_step.runTest", forward_step.runTest, envir = .GlobalEnv)
  withr::defer(rm(forward_step.runTest, envir = .GlobalEnv))

  res <- core_write(x = list(array(1, dim = c(1,1,1)), array(2, dim = c(1,1,1))), transforms = "runTest")

  expect_equal(captured$run_ids, c("run-01", "run-02"))
  expect_equal(captured$current_run, c("run-01", "run-02"))
  expect_equal(captured$values, c(1, 2))
  expect_equal(res$handle$run_ids, c("run-01", "run-02"))
  expect_equal(res$handle$current_run_id, "run-01") # current_run_id should be the first generated ID
})

test_that("mask is validated and stored", {
  arr <- array(1, dim = c(2,2,2,3))
  msk <- array(TRUE, dim = c(2,2,2))
  res <- core_write(x = arr, transforms = "tA", mask = msk)
  expect_equal(res$handle$mask_info$active_voxels, sum(msk))
  expect_true(all(res$handle$meta$mask == msk))
})

test_that("mask voxel mismatch triggers error", {
  arr <- array(1, dim = c(2,2,2,1))
  bad <- array(TRUE, dim = c(2,2,3))
  expect_error(
    core_write(x = arr, transforms = "tA", mask = bad),
    class = "lna_error_validation"
  )
})

test_that("core_write works with progress handlers", {
  old_handlers <- progressr::handlers()
  withr::defer(progressr::handlers(old_handlers))

  progressr::handlers(progressr::handler_void())
  expect_no_error(
    progressr::with_progress(
      core_write(x = array(1, dim = c(1, 1, 1)), transforms = c("tA"))
    )
  )
})

test_that("input data is promoted to required dimensions", {
  expect_no_error(
    core_write(x = matrix(1:4, nrow = 2), transforms = "tA")
  )
})

test_that("header and plugins must be named lists", {
  arr <- array(1, dim = c(2,2,2))
  expect_error(
    core_write(x = arr, transforms = "tA", header = list(1)),
    class = "lna_error_validation"
  )
  expect_error(
    core_write(x = arr, transforms = "tA", plugins = list(1)),
    class = "lna_error_validation"
  )
})
