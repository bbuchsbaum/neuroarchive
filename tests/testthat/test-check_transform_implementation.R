library(testthat)
library(neuroarchive)
# Tests for check_transform_implementation

test_that("warnings emitted for collisions", {
  expect_warning(check_transform_implementation("quant"), "collides")
  expect_warning(check_transform_implementation("stats"), "collides")
})

test_that("no warning for unique name", {
  expect_warning(check_transform_implementation("myunique"), NA)
})
