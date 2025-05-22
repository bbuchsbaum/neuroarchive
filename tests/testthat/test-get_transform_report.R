library(hdf5r)

test_that("lna_get_transform_report retrieves report", {
  arr <- array(runif(6), dim = c(2,3))
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(arr, file = tmp, transforms = "quant")
  rep <- lna_get_transform_report(tmp, 0)
  expect_type(rep, "list")
  expect_equal(rep$report_version, "1.0")
})

test_that("lna_get_transform_report handles missing report_path and gzip", {
  arr <- array(runif(6), dim = c(2,3))
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(arr, file = tmp, transforms = "quant")
  h5 <- H5File$new(tmp, mode = "a")
  dset <- h5[["/transforms/00_quant.json"]]
  desc <- jsonlite::fromJSON(dset$read())
  desc$params$report_path <- NULL
  dset$write(jsonlite::toJSON(desc, auto_unbox = TRUE))
  dset$close(); h5$close_all()
  rep <- lna_get_transform_report(tmp, "00_quant.json")
  expect_type(rep, "list")
  expect_true("clipped_samples_count" %in% names(rep))
})

test_that("lna_get_quant_report is wrapper", {
  arr <- array(runif(6), dim = c(2,3))
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(arr, file = tmp, transforms = "quant")
  rep1 <- lna_get_transform_report(tmp, 0)
  rep2 <- lna_get_quant_report(tmp, 0)
  expect_identical(rep1, rep2)
})
