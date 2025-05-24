library(testthat)
library(withr)

# Basic integration of LNAFacade

test_that("LNAFacade writes and reads", {
  fac <- LNAFacade$new()
  tmp <- local_tempfile(fileext = ".h5")
  result <- fac$write(array(1, dim = c(1,1,1)), tmp, transforms = character())
  expect_true(file.exists(tmp))
  expect_identical(result$file, fac$last_output)
  expect_s3_class(result, "lna_write_result")
  hdl <- fac$read(tmp)
  expect_true(inherits(hdl, "DataHandle"))
})
