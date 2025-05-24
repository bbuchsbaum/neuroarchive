library(testthat)

# simple template adding quant then pca
simple_template <- function(pipe, bits = 8, k = 2) {
  pipe <- quant(pipe, bits = bits)
  pipe <- pca(pipe, k = k)
  pipe
}

register_lna_template("simple", simple_template, force = TRUE)

arr <- array(rnorm(10), dim = c(5,2))

test_that("apply_template applies registered template", {
  pipe <- as_pipeline(arr)
  pipe2 <- apply_template(pipe, "simple")
  expect_equal(length(pipe2$steps()), 2L)
  expect_equal(pipe2$steps()[[1]]$type, "quant")
  expect_equal(pipe2$steps()[[2]]$type, "basis")
})

test_that("apply_template overrides parameters", {
  # Create fresh pipeline for this test
  pipe <- as_pipeline(arr)
  pipe2 <- apply_template(pipe, "simple", quant.bits = 4)
  expect_equal(pipe2$steps()[[1]]$params$bits, 4)

  # Create another fresh pipeline for the second part
  pipe_fresh <- as_pipeline(arr)
  pipe3 <- apply_template(pipe_fresh, "simple", quant = list(bits = 6))
  expect_equal(pipe3$steps()[[1]]$params$bits, 6)
})


