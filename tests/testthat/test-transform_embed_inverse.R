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

  expect_true(h$has_key("dense_mat"))
  expect_false(h$has_key("coef"))
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


test_that("invert_step.embed errors when datasets are missing", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- H5File$new(tmp, mode = "w")
  desc <- list(
    type = "embed",
    params = list(basis_path = "/missing/matrix"),
    inputs = c("dense"),
    outputs = c("coef")
  )
  handle <- DataHandle$new(initial_stash = list(coef = matrix(0, nrow = 1, ncol = 1)), h5 = h5)

  expect_error(
    invert_step.embed("embed", desc, handle),
    class = "lna_error_contract",
    regexp = "not found"
    )
  }

test_that("invert_step.embed errors when datasets missing", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- H5File$new(tmp, mode = "w")
  desc <- list(type = "embed", params = list(basis_path = "/missing"),
               inputs = c("dense"), outputs = c("coef"))
  handle <- DataHandle$new(initial_stash = list(coef = matrix(1)), h5 = h5)
  expect_error(
    invert_step.embed("embed", desc, handle),
    class = "lna_error_contract",
    regexp = "basis"

  )
  h5$close_all()
})

test_that("invert_step.embed applies scaling and centering", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- H5File$new(tmp, mode = "w")
  basis_mat <- diag(2)
  center_vec <- c(5, 10)
  scale_vec <- c(2, 4)
  neuroarchive:::h5_write_dataset(h5, "/basis/mat", basis_mat)
  neuroarchive:::h5_write_dataset(h5, "/basis/center", center_vec)
  neuroarchive:::h5_write_dataset(h5, "/basis/scale", scale_vec)

  desc <- list(
    type = "embed",
    params = list(
      basis_path = "/basis/mat",
      center_data_with = "/basis/center",
      scale_data_with = "/basis/scale"
    ),
    inputs = c("dense_mat"),
    outputs = c("coef")
  )

  coef_mat <- matrix(c(1,2,3,4), nrow = 2)
  handle <- DataHandle$new(initial_stash = list(coef = coef_mat), h5 = h5)

  h <- invert_step.embed("embed", desc, handle)

  expected <- sweep(coef_mat %*% basis_mat, 2, scale_vec, "*")
  expected <- sweep(expected, 2, center_vec, "+")
  expect_equal(h$stash$dense_mat, expected)
})
          
test_that("invert_step.embed applies center and scale", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- H5File$new(tmp, mode = "w")
  basis_mat <- diag(2)
  neuroarchive:::h5_write_dataset(h5, "/b/mat", basis_mat)
  neuroarchive:::h5_write_dataset(h5, "/b/center", c(1,2))
  neuroarchive:::h5_write_dataset(h5, "/b/scale", c(2,2))
  desc <- list(
    type = "embed",
    params = list(basis_path = "/b/mat", center_data_with = "/b/center",
                  scale_data_with = "/b/scale"),
    inputs = c("dense"), outputs = c("coef")
  )
  coef_mat <- matrix(c(1,1,1,1), nrow = 2)
  handle <- DataHandle$new(initial_stash = list(coef = coef_mat), h5 = h5)
  h <- invert_step.embed("embed", desc, handle)
  expected <- coef_mat %*% basis_mat
  expected <- sweep(expected, 2, c(2,2), "*")
  expected <- sweep(expected, 2, c(1,2), "+")
  expect_equal(h$stash$dense, expected)

  h5$close_all()
})
