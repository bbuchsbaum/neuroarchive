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
