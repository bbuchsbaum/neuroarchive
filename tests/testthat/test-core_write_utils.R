library(testthat)

# ensure_lna_array_input -------------------------------------------------

test_that("ensure_lna_array_input converts vectors to 4D with metadata", {
  vec <- 1:5
  arr <- neuroarchive:::ensure_lna_array_input(vec)
  expect_equal(dim(arr), c(5, 1, 1, 1))
  expect_false(is.null(attr(arr, "lna.orig_dims")))
  expect_identical(attr(arr, "lna.orig_dims"), 5L)
  expect_false(attr(arr, "lna.was_3d"))
})

# validate_named_list ----------------------------------------------------

test_that("validate_named_list rejects unnamed lists", {
  expect_error(
    neuroarchive:::validate_named_list(list(1, 2), "header"),
    class = "lna_error_validation"
  )

  expect_identical(
    neuroarchive:::validate_named_list(list(a = 1, b = 2), "header"),
    list(a = 1, b = 2)
  )
})
