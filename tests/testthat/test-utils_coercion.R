library(testthat)

# Tests for as_dense_mat.array

test_that("as_dense_mat handles 3D arrays", {
  arr <- array(seq_len(24), dim = c(2, 3, 4))
  mat <- neuroarchive:::as_dense_mat(arr)
  expect_equal(dim(mat), c(4, 6))
  expect_equal(mat[1, ], as.numeric(arr[, , 1]))
  expect_equal(mat[4, ], as.numeric(arr[, , 4]))
})

test_that("as_dense_mat handles 4D arrays", {
  arr <- array(seq_len(2 * 3 * 4 * 5), dim = c(2, 3, 4, 5))
  mat <- neuroarchive:::as_dense_mat(arr)
  expect_equal(dim(mat), c(5, 24))
  expect_equal(mat[1, ], as.numeric(arr[, , , 1]))
  expect_equal(mat[5, ], as.numeric(arr[, , , 5]))
})

test_that("as_vol4d passes through 4D arrays and errors otherwise", {
  arr4d <- array(seq_len(2 * 3 * 4 * 5), dim = c(2, 3, 4, 5))
  expect_identical(neuroarchive:::as_vol4d(arr4d), arr4d)

  arr3d <- array(seq_len(2 * 3 * 4), dim = c(2, 3, 4))
  expect_error(neuroarchive:::as_vol4d(arr3d), class = "lna_error_validation")
})
