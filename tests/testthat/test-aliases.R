library(testthat)

# Tests for convenience alias functions

test_that("compress_fmri forwards to write_lna", {
  captured <- NULL
  with_mocked_bindings(
    write_lna = function(...) { captured <<- list(...); "res" },
    {
      out <- compress_fmri(x = 1, file = "foo.h5")
    }
  )
  expect_identical(captured$x, 1)
  expect_identical(captured$file, "foo.h5")
  expect_identical(out, "res")
})

test_that("open_lna forwards to read_lna", {
  captured <- NULL
  with_mocked_bindings(
    read_lna = function(...) { captured <<- list(...); "out" },
    {
      res <- open_lna(file = "bar.h5", lazy = TRUE)
    }
  )
  expect_identical(captured$file, "bar.h5")
  expect_true(captured$lazy)
  expect_identical(res, "out")
})
