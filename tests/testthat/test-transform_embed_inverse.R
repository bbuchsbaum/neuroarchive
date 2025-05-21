library(testthat)
library(hdf5r)
library(neuroarchive)
library(withr)


test_that("invert_step.embed reconstructs dense data", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- H5File$new(tmp, mode = "w")
  basis_mat <- matrix(c(1,0,0,1), nrow = 2)
  neuroarchive:::h5_write_dataset(h5, "/basis/test/matrix", basis_mat)

  desc <- list(
    type = "embed",
    params = list(basis_path = "/basis/test/matrix"),
    inputs = c("dense_mat"),
    outputs = c("coef")
  )

  coef_mat <- matrix(c(1,2,3,4), nrow = 2)
  handle <- DataHandle$new(initial_stash = list(coef = coef_mat), h5 = h5)

  h <- invert_step.embed("embed", desc, handle)

  expect_true(h$exists("dense_mat"))
  expect_false(h$exists("coef"))
  expected <- tcrossprod(coef_mat, basis_mat)
  expect_equal(h$stash$dense_mat, expected)

  h5$close_all()
})

test_that("read_lna applies roi_mask and time_idx for embed", {
  arr <- matrix(seq_len(20), nrow = 5, ncol = 4)
  tmp <- local_tempfile(fileext = ".h5")

  write_lna(arr, file = tmp, transforms = c("basis", "embed"),
            transform_params = list(
              embed = list(basis_path = "/basis/00_basis/matrix",
                            center_data_with = "/basis/00_basis/center")
            ))

  roi <- c(TRUE, FALSE, TRUE, FALSE)
  h <- read_lna(tmp, roi_mask = roi, time_idx = c(2,4))
  out <- h$stash$input

  expect_equal(dim(out), c(2, sum(roi)))
  expect_equal(out, arr[c(2,4), roi])
})
