library(testthat)
library(hdf5r)
#library(neuroarchive)
library(withr)


test_that("invert_step.basis reconstructs dense data for both storage orders", {
  coef_mat <- matrix(c(1,2,3,4), nrow = 2)

  base_mat_cxv <- matrix(c(1,0,0,1,1,1), nrow = 2) # component x voxel
  base_mat_vxc <- t(base_mat_cxv)                      # voxel x component

  for (ord in c("component_x_voxel", "voxel_x_component")) {
    tmp <- tempfile(fileext = ".h5")
    on.exit(unlink(tmp), add = TRUE)
    h5 <- H5File$new(tmp, mode = "w")
    mat <- if (identical(ord, "component_x_voxel")) base_mat_cxv else base_mat_vxc
    neuroarchive:::h5_write_dataset(h5[["/"]], "/basis/test/matrix", mat)

    desc <- list(
      type = "basis",
      params = list(storage_order = ord),
      datasets = list(list(path = "/basis/test/matrix", role = "basis_matrix")),
      inputs = c("dense_mat"),
      outputs = c("coef")
    )

    handle <- DataHandle$new(initial_stash = list(coef = coef_mat), h5 = h5)
    h <- invert_step.basis("basis", desc, handle)
    expect_true(h$has_key("dense_mat"))
    expect_false(h$has_key("coef"))
    expected <- if (identical(ord, "component_x_voxel"))
      coef_mat %*% base_mat_cxv else coef_mat %*% t(base_mat_vxc)
    expect_equal(h$stash$dense_mat, expected)
    h5$close_all()
  }

  # Commenting out this potentially problematic block for now
  # handle <- DataHandle$new(initial_stash = list(coef = coef_mat), h5 = h5)
  # 
  # h <- invert_step.basis("basis", desc, handle)
  # 
  # expect_true(h$has_key("dense_mat"))
  # expect_false(h$has_key("coef"))
  # expected <- tcrossprod(coef_mat, basis_mat) # basis_mat not defined here
  # expect_equal(h$stash$dense_mat, expected)
  # 
  # h5$close_all()

})

test_that("invert_step.basis applies subset", {
  coef_mat <- matrix(1:6, nrow = 3, ncol = 2)
  subset <- list(roi_mask = c(TRUE, FALSE, TRUE), time_idx = c(1,3))
  base_mat_cxv <- matrix(c(1,0,0,1,1,1), nrow = 2)
  base_mat_vxc <- t(base_mat_cxv)

  for (ord in c("component_x_voxel", "voxel_x_component")) {
    tmp <- tempfile(fileext = ".h5")
    on.exit(unlink(tmp), add = TRUE)
    h5 <- H5File$new(tmp, mode = "w")
    mat <- if (identical(ord, "component_x_voxel")) base_mat_cxv else base_mat_vxc
    neuroarchive:::h5_write_dataset(h5[["/"]], "/basis/test/matrix", mat)
    desc <- list(
      type = "basis",
      params = list(storage_order = ord),
      datasets = list(list(path = "/basis/test/matrix", role = "basis_matrix")),
      inputs = c("dense_mat"),
      outputs = c("coef")
    )
    handle <- DataHandle$new(initial_stash = list(coef = coef_mat), h5 = h5, subset = subset)
    h <- invert_step.basis("basis", desc, handle)
    expect_equal(dim(h$stash$dense_mat),
                 c(length(subset$time_idx), sum(subset$roi_mask)))
    h5$close_all()
  }
})

test_that("invert_step.basis honours subset$roi and subset$time", {
  coef_mat <- matrix(1:6, nrow = 3, ncol = 2)
  subset <- list(roi = c(TRUE, FALSE, TRUE), time = c(1,3))
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)
  h5 <- H5File$new(tmp, mode = "w")
  neuroarchive:::h5_write_dataset(h5[["/"]], "/basis/test/matrix", matrix(c(1,0,0,1,1,1), nrow = 2))
  desc <- list(
    type = "basis",
    params = list(storage_order = "component_x_voxel"),
    datasets = list(list(path = "/basis/test/matrix", role = "basis_matrix")),
    inputs = c("dense_mat"),
    outputs = c("coef")
  )
  handle <- DataHandle$new(initial_stash = list(coef = coef_mat), h5 = h5, subset = subset)
  h <- invert_step.basis("basis", desc, handle)
  expect_equal(dim(h$stash$dense_mat), c(length(subset$time), sum(subset$roi)))
  h5$close_all()
})

test_that("invert_step.basis validates storage_order", {
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)
  h5 <- H5File$new(tmp, mode = "w")
  neuroarchive:::h5_write_dataset(h5[["/"]], "/basis/test/matrix", matrix(1))
  desc <- list(
    type = "basis",
    params = list(storage_order = "bogus"),
    datasets = list(list(path = "/basis/test/matrix", role = "basis_matrix")),
    inputs = c("dense_mat"),
    outputs = c("coef")
  )
  handle <- DataHandle$new(initial_stash = list(coef = matrix(1)), h5 = h5)
  expect_error(
    invert_step.basis("basis", desc, handle),
    class = "lna_error_validation"
  )
  h5$close_all()
})

