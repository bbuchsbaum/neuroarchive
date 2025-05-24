library(testthat)
library(withr)

# Basic integration of LNAFacade

test_that("LNAFacade writes and reads", {
  fac <- LNAFacade$new()
  tmp <- local_tempfile(fileext = ".h5")
  invisible(fac$write(array(1, dim = c(1,1,1)), tmp, transforms = character()))
  expect_true(file.exists(tmp))
  expect_identical(fac$last_output, tmp)
  res <- fac$read(tmp)
  expect_true(inherits(res, "DataHandle"))
})

test_that("set_defaults replaces defaults for write", {
  fac <- LNAFacade$new()
  fac$set_defaults(list(quant = list(bits = 7L)))
  tmp <- local_tempfile(fileext = ".h5")
  fac$write(array(1, dim = c(1,1,1)), tmp, transforms = "quant")
  h5 <- hdf5r::H5File$new(tmp, mode = "r")
  dset <- h5[["/transforms/00_quant.json"]]
  desc <- jsonlite::fromJSON(dset$read())
  h5$close_all()
  expect_equal(desc$params$bits, 7L)
})
