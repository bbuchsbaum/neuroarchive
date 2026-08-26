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

test_that("invert_step.basis passes through dense matrices in both orientations", {
  dense_full <- matrix(1:12, nrow = 3, ncol = 4)
  subset <- list(roi_mask = c(TRUE, FALSE, TRUE, FALSE), time_idx = c(1, 3))
  expected <- dense_full[subset$time_idx, c(1, 3), drop = FALSE]

  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)
  h5 <- H5File$new(tmp, mode = "w")
  neuroarchive:::h5_write_dataset(h5[["/"]], "/basis/test/matrix", matrix(1:8, nrow = 2))

  desc <- list(
    type = "basis",
    params = list(storage_order = "component_x_voxel"),
    datasets = list(list(path = "/basis/test/matrix", role = "basis_matrix")),
    inputs = c("dense_mat"),
    outputs = c("coef")
  )

  handle_txv <- DataHandle$new(initial_stash = list(coef = dense_full), h5 = h5, subset = subset)
  out_txv <- invert_step.basis("basis", desc, handle_txv)
  expect_equal(out_txv$stash$dense_mat, expected)

  handle_vxt <- DataHandle$new(initial_stash = list(coef = t(dense_full)), h5 = h5, subset = subset)
  out_vxt <- invert_step.basis("basis", desc, handle_vxt)
  expect_equal(out_vxt$stash$dense_mat, expected)

  h5$close_all()
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

test_that("read_lna applies roi_mask and time_idx for basis + quant", {
  set.seed(11)
  arr <- array(rnorm(2 * 2 * 1 * 6), dim = c(2, 2, 1, 6))
  tmp <- local_tempfile(fileext = ".h5")

  write_lna(
    arr,
    file = tmp,
    transforms = c("basis", "quant"),
    transform_params = list(
      basis = list(k = 1, allow_dense_passthrough = TRUE),
      quant = list(bits = 12)
    )
  )

  full <- as_dense_mat(read_lna(tmp)$stash$input)
  roi <- array(c(TRUE, FALSE, TRUE, FALSE), dim = c(2, 2, 1))
  idx <- c(2, 5)
  vox_idx <- which(as.logical(roi))
  expected <- full[idx, vox_idx, drop = FALSE]

  h <- read_lna(tmp, roi_mask = roi, time_idx = idx)
  out <- as_dense_mat(h$stash$input)
  expect_equal(out, expected)
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
