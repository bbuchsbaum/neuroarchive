library(testthat)
library(hdf5r)
library(neuroarchive)
library(withr)


test_that("invert_step.basis reconstructs dense data", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- H5File$new(tmp, mode = "w")
  basis_mat <- matrix(c(1,0,0,1,1,1), nrow = 2)
  neuroarchive:::h5_write_dataset(h5, "/basis/test/matrix", basis_mat)

  desc <- list(
    type = "basis",
    params = list(storage_order = "component_x_voxel"),
    datasets = list(list(path = "/basis/test/matrix", role = "basis_matrix")),
    inputs = c("dense_mat"),
    outputs = c("coef")
  )

  coef_mat <- matrix(c(1,2,3,4), nrow = 2)
  handle <- DataHandle$new(initial_stash = list(coef = coef_mat), h5 = h5)

  h <- invert_step.basis("basis", desc, handle)

  expect_true(h$exists("dense_mat"))
  expect_false(h$exists("coef"))
  expected <- tcrossprod(coef_mat, basis_mat)
  expect_equal(h$stash$dense_mat, expected)

  h5$close_all()
})

