library(testthat)
library(hdf5r)
library(withr)

# Basic functionality test using actual file

test_that("write_lna with file=NULL uses in-memory HDF5", {
  result <- write_lna(x = 1, file = NULL, transforms = character(0))
  expect_s3_class(result, "lna_write_result")
  expect_null(result$file)
  expect_true(inherits(result$plan, "Plan"))
})

test_that("write_lna writes header attributes to file", {
  tmp <- local_tempfile(fileext = ".h5")
  result <- write_lna(x = 1, file = tmp, transforms = character(0),
                      header = list(foo = 2L))
  expect_true(file.exists(tmp))
  h5 <- neuroarchive:::open_h5(tmp, mode = "r")
  grp <- h5[["header/global"]]
  expect_identical(h5_attr_read(grp, "foo"), 2L)
  neuroarchive:::close_h5_safely(h5)
})

# Parameter forwarding for write_lna

test_that("write_lna forwards arguments to core_write and materialise_plan", {
  captured <- list()
  fake_plan <- Plan$new()
  fake_handle <- DataHandle$new()

  with_mocked_bindings(
    core_write = function(x, transforms, transform_params, mask = NULL, header = NULL) {
      captured$core <- list(x = x, transforms = transforms,
                            transform_params = transform_params,
                            header = header)
      list(handle = fake_handle, plan = fake_plan)
    },
    materialise_plan = function(h5, plan, checksum = "none", header = NULL) {
      captured$mat <- list(is_h5 = inherits(h5, "H5File"), plan = plan,
                           header = header)
    }, {
      write_lna(x = 42, file = tempfile(fileext = ".h5"),
                transforms = c("tA"),
                transform_params = list(tA = list(foo = "bar")),
                header = list(a = 1))
    }
  )

  expect_equal(captured$core$x, 42)
  expect_equal(captured$core$transforms, c("tA"))
  expect_equal(captured$core$transform_params, list(tA = list(foo = "bar")))
  expect_true(captured$mat$is_h5)
  expect_identical(captured$mat$plan, fake_plan)
  expect_equal(captured$core$header, list(a = 1))
  expect_equal(captured$mat$header, list(a = 1))
})

# Parameter forwarding for read_lna

test_that("read_lna forwards arguments to core_read", {
  captured <- list()
  with_mocked_bindings(
    core_read = function(file, allow_plugins, validate, output_dtype, lazy) {
      captured$core <- list(file = file, allow_plugins = allow_plugins,
                            validate = validate, output_dtype = output_dtype,
                            lazy = lazy)
      DataHandle$new()
    }, {
      read_lna("somefile.h5", allow_plugins = "on", validate = TRUE,
               output_dtype = "float64", lazy = TRUE)
    }
  )

  expect_equal(captured$core$file, "somefile.h5")
  expect_equal(captured$core$allow_plugins, "on")
  expect_true(captured$core$validate)
  expect_equal(captured$core$output_dtype, "float64")
  expect_true(captured$core$lazy)
})

test_that("read_lna lazy=TRUE keeps file open", {
  result <- write_lna(x = 1, file = NULL, transforms = character(0))
  handle <- read_lna(result$file, lazy = TRUE)
  expect_true(handle$h5$is_valid())
  neuroarchive:::close_h5_safely(handle$h5)
})
