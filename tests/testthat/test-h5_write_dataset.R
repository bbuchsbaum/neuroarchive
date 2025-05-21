library(testthat)
library(hdf5r)
library(withr)

# Test writing a simple numeric matrix with chunking and compression

test_that("h5_write_dataset writes dataset with compression", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- neuroarchive:::open_h5(tmp, mode = "w")
  root <- h5[["/"]]

  mat <- matrix(1:9, nrow = 3)
  h5_write_dataset(root, "/group/data", mat, chunk_dims = c(2,2), compression_level = 6)

  expect_true(root$exists("group/data"))
  dset <- root[["/group/data"]]
  dcpl <- dset$get_create_plist()
  expect_equal(dcpl$get_chunk(2), c(2,2))
  expect_equal(dcpl$get_filter(0)$filter, hdf5r::h5const$H5Z_FILTER_DEFLATE)

  read_back <- dset$read()
  expect_equal(read_back, mat)

  dcpl$close()
  dset$close()
  neuroarchive:::close_h5_safely(h5)
})
