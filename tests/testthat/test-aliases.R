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

test_that("open_lna is an alias of read_lna", {
  expect_identical(open_lna, read_lna)
})
