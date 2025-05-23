library(testthat)
library(hdf5r)
library(withr)

# Basic creation of empty dataset

test_that("h5_create_empty_dataset creates typed dataset", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- open_h5(tmp, mode = "w")
  root <- h5[["/"]]
  dset <- h5_create_empty_dataset(root, "grp/data", c(4,3,2), dtype = "float32", chunk_dims = c(2,2,1))
  
  expect_true(root$exists("grp/data"))
  dcpl <- dset$get_create_plist()
  
  # Use direct R6 method dset$dims to avoid S3 dispatch issues with dim()
  expect_equal(dset$dims, c(4,3,2))
  expect_equal(dcpl$get_chunk(3), c(2,2,1))
  actual_dtype_text <- dset$get_type()$to_text()
  expect_equal(actual_dtype_text, "H5T_IEEE_F32LE")
  dcpl$close(); dset$close(); close_h5_safely(h5)
})
