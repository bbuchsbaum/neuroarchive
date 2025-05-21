library(testthat)
#library(neuroarchive)
library(withr)


test_that("basis -> embed -> quant pipeline roundtrip", {
  arr <- array(runif(40), dim = c(2,2,2,5))
  tmp <- local_tempfile(fileext = ".h5")

  res <- write_lna(arr, file = tmp,
                   transforms = c("basis", "embed", "quant"),
                   transform_params = list(
                     embed = list(basis_path = "/basis/00_basis/matrix",
                                   center_data_with = "/basis/00_basis/center")
                   ))
  expect_true(file.exists(tmp))

  h <- read_lna(tmp)
  out <- h$stash$input
  expect_equal(dim(out), dim(arr))
  expect_lt(mean(abs(out - arr)), 1)
})


test_that("quant only pipeline roundtrip", {
  arr <- array(runif(12), dim = c(3,4))
  tmp <- local_tempfile(fileext = ".h5")

  res <- write_lna(arr, file = tmp, transforms = "quant")
  expect_true(file.exists(tmp))

  h <- read_lna(tmp)
  out <- h$stash$input
  expect_equal(dim(out), dim(arr))
  expect_lt(mean(abs(out - arr)), 1)
})


test_that("lna_reader works for multi-transform pipeline", {
  arr <- array(runif(40), dim = c(2,2,2,5))
  tmp <- local_tempfile(fileext = ".h5")

  res <- write_lna(arr, file = tmp,
                   transforms = c("basis", "embed", "quant"),
                   transform_params = list(
                     embed = list(basis_path = "/basis/00_basis/matrix",
                                   center_data_with = "/basis/00_basis/center")
                   ))

  reader <- read_lna(tmp, lazy = TRUE)
  out <- reader$data()$stash$input
  expect_equal(dim(out), dim(arr))
  reader$close()
})
