library(testthat)

# Tests for convenience alias functions

test_that("compress_fmri forwards to write_lna", {
  captured <- NULL
  local_mocked_bindings(
    write_lna = function(...) { captured <<- list(...); "res" },
    .env = asNamespace("neuroarchive")
  )
  out <- compress_fmri(x = 1, file = "foo.h5")
  expect_identical(captured$x, 1)
  expect_identical(captured$file, "foo.h5")
  expect_identical(out, "res")
})

test_that("open_lna is an alias of read_lna", {
  # Explicitly reference functions from the namespace to ensure they are found
  # This assumes devtools::load_all() has correctly loaded the package.
  expect_identical(neuroarchive::open_lna, neuroarchive::read_lna)
})
