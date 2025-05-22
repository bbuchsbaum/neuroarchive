library(hdf5r)

test_that("lna_get_transform_report retrieves report", {
  arr <- array(runif(6), dim = c(2,3))
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(arr, file = tmp, transforms = "quant")
  rep <- lna_get_transform_report(tmp, 0)
  expect_type(rep, "list")
  expect_equal(rep$report_version, "1.0")
})
