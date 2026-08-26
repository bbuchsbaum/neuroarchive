library(testthat)
library(hdf5r)

test_that("list_h5_filters returns filter names and shuffle toggle works", {
  tmp <- tempfile(fileext = ".h5")
  on.exit(unlink(tmp), add = TRUE)
  f <- H5File$new(tmp, mode = "w"); root <- f[["/"]]

  # Ensure shuffle is off by default
  neuroarchive::lna_options(write.quant.nbit.shuffle = FALSE,
                            write.quant.nbit.deflate_level = 0L)
  d1 <- neuroarchive:::h5_write_qbits_nbit(root, "/qbits1", as.integer(0:31) %% 8L, bits = 3)
  fl1 <- neuroarchive::list_h5_filters(d1)
  expect_true("nbit" %in% fl1)
  expect_false("shuffle" %in% fl1)
  d1$close()

  # Turn shuffle on and verify
  neuroarchive::lna_options(write.quant.nbit.shuffle = TRUE)
  d2 <- neuroarchive:::h5_write_qbits_nbit(root, "/qbits2", as.integer(0:31) %% 8L, bits = 3)
  fl2 <- neuroarchive::list_h5_filters(d2)
  expect_true("nbit" %in% fl2)
  expect_true("shuffle" %in% fl2)
  d2$close(); f$close_all()
})
