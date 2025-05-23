library(testthat)
library(hdf5r)
library(neuroarchive)

# ensure first_vals dataset dimension is stored and read correctly

test_that("first_vals handled for 1D input", {
  vec <- 1:5
  tmp <- local_tempfile(fileext = ".h5")
  res <- write_lna(vec, file = tmp, transforms = "delta",
                   transform_params = list(delta = list(axis = 1)))
  h5 <- H5File$new(tmp, mode = "r")
  first_path <- res$plan$datasets$path[res$plan$datasets$role == "first_values"]
  fv <- h5[[first_path]]$read()
  h5$close_all()
  expect_equal(dim(matrix(fv, nrow = 1)), c(1, 1))

  h <- read_lna(tmp)
  expect_equal(drop(h$stash$input), vec)
})

test_that("first_vals handled for 2D axis=1", {
  mat <- matrix(1:12, nrow = 3, ncol = 4)
  tmp <- local_tempfile(fileext = ".h5")
  res <- write_lna(mat, file = tmp, transforms = "delta",
                   transform_params = list(delta = list(axis = 1)))
  h5 <- H5File$new(tmp, mode = "r")
  first_path <- res$plan$datasets$path[res$plan$datasets$role == "first_values"]
  fv <- h5[[first_path]]$read()
  h5$close_all()
  expect_equal(dim(matrix(fv, nrow = 1)), c(1, ncol(mat)))

  h <- read_lna(tmp)
  expect_equal(drop(h$stash$input), mat)
})

test_that("first_vals handled for 3D axis>1", {
  arr <- array(seq_len(24), dim = c(2, 3, 4))
  tmp <- local_tempfile(fileext = ".h5")
  res <- write_lna(arr, file = tmp, transforms = "delta",
                   transform_params = list(delta = list(axis = 2)))
  h5 <- H5File$new(tmp, mode = "r")
  first_path <- res$plan$datasets$path[res$plan$datasets$role == "first_values"]
  fv <- h5[[first_path]]$read()
  h5$close_all()
  expect_equal(dim(matrix(fv, nrow = 1)), c(1, prod(dim(arr)[-2])))

  h <- read_lna(tmp)
  expect_equal(drop(h$stash$input), arr)
})
