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
  result <- core_write(x = 1, transforms = c("tA", "tB"), transform_params = list(tB = list(foo = "bar")))
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

  with_mocked_bindings(
    default_params = function(type) {
      list(a = 1, b = 2, nested = list(x = 0, y = 0))
    }, {
      res <- core_write(
        x = 1,
        transforms = c("tB"),
        transform_params = list(tB = list(b = 20, nested = list(y = 5)))
      )
    }
  )

  expect_equal(
    res$plan$descriptors[[1]]$params,
    list(a = 10, b = 20, nested = list(x = 1, y = 5))
  )
})

test_that("unknown transform names in transform_params error", {
  expect_error(
    core_write(x = 1, transforms = c("tA"), transform_params = list(tB = list())),
    class = "lna_error_validation"
  )
})
