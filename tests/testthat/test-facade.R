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

test_that("LNAFacade$read uses last_output when file missing", {
  fac <- LNAFacade$new()
  tmp <- local_tempfile(fileext = ".h5")
  invisible(fac$write(array(2, dim = c(1,1,1)), tmp, transforms = character()))
  expect_identical(fac$last_output, tmp)
  res <- fac$read()
  expect_true(inherits(res, "DataHandle"))
})


test_that("LNAFacade$write_pipeline forwards to lna_write", {
  fac <- LNAFacade$new()
  arr <- array(1, dim = c(1,1,1))
  pipe <- as_pipeline(arr)

  captured <- list()
  fake_result <- list(file = "out.h5")

  local_mocked_bindings(
    lna_write = function(pipe_arg, file, ...) {
      captured$pipe <<- pipe_arg
      captured$file <<- file
      fake_result
    },
    .env = asNamespace("neuroarchive")
  )

  res <- fac$write_pipeline(pipe, file = "out.h5")

  expect_identical(captured$pipe, pipe)
  expect_equal(captured$file, "out.h5")
  expect_equal(fac$last_output, "out.h5")
  expect_identical(res, "out.h5")
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
