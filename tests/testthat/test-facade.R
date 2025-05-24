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

test_that("LNAFacade$read uses last_output when file missing", {
  fac <- LNAFacade$new()
  tmp <- local_tempfile(fileext = ".h5")
  invisible(fac$write(array(2, dim = c(1,1,1)), tmp, transforms = character()))
  expect_identical(fac$last_output, tmp)
  res <- fac$read()
  expect_true(inherits(res, "DataHandle"))
})
