library(testthat)

# Simple chain of DSL verbs should be expressive and easy to read

test_that("pipeline verbs chain cleanly", {
  arr <- array(runif(20), dim = c(5, 4))

  pipe <- as_pipeline(arr) |> 
    pca(k = 2) |> 
    embed() |> 
    quant(bits = 6)

  step_types <- vapply(pipe$steps(), `[[`, character(1), "type")
  expect_equal(step_types, c("basis", "embed.pca", "quant"))

  captured <- list()
  local_mocked_bindings(
    write_lna = function(x, file, transforms, transform_params, run_id, checksum = "none") {
      captured$transforms <<- transforms
      captured$params <<- transform_params
      list(ok = TRUE)
    },
    .env = asNamespace("neuroarchive")
  )

  lna_write(pipe, file = "out.h5")
  expect_equal(captured$transforms, c("basis", "embed.pca", "quant"))
  expect_equal(captured$params$quant$bits, 6)
})

# Demonstrate that starting from a list of runs is equally concise

test_that("list input works seamlessly in DSL", {
  lst <- list(
    run1 = array(1:8, dim = c(2,2,2)),
    run2 = array(9:16, dim = c(2,2,2))
  )

  pipe <- as_pipeline(lst) |> delta(order = 1)

  expect_equal(pipe$runs, c("run1", "run2"))
  step <- pipe$get_last_step_spec()
  expect_equal(step$type, "delta")
  expect_equal(step$params$order, 1)
})

