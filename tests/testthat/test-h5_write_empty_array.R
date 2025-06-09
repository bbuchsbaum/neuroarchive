library(testthat)
library(hdf5r)

# Test writing and reading an empty array using placeholder logic

test_that("h5_write_dataset preserves empty array via placeholder", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- H5File$new(tmp, mode = "w")
  root <- h5[["/"]]

  arr <- array(numeric(0), dim = c(0, 3, 2))
  neuroarchive:::h5_write_dataset(root, "empty/data", arr)

  dset <- root[["empty/data"]]
  expect_true(neuroarchive:::h5_attr_exists(dset, "lna_empty_array_placeholder"))

  read_back <- neuroarchive:::h5_read(root, "empty/data")
  expect_equal(dim(read_back), dim(arr))
  expect_equal(length(read_back), 0)

  dset$close()
  h5$close_all()
})
