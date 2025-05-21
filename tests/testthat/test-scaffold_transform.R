library(testthat)
library(withr)

# Test scaffold_transform creates files with expected content

test_that("scaffold_transform creates template files", {
  tmp <- local_tempdir()
  withr::local_dir(tmp)
  paths <- scaffold_transform("mycustom")

  expect_true(file.exists(paths$r_file))
  expect_true(file.exists(paths$schema))
  expect_true(file.exists(paths$test))

  r_lines <- readLines(paths$r_file)
  expect_true(any(grepl("forward_step.mycustom", r_lines, fixed = TRUE)))
  expect_true(any(grepl("default_params('mycustom')", r_lines, fixed = TRUE)))
})

test_that("scaffold_transform warns on namespace collisions", {
  tmp <- local_tempdir()
  withr::local_dir(tmp)
  expect_warning(scaffold_transform("delta"), "namespace")
})
