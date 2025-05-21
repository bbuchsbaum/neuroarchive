library(testthat)
library(hdf5r)
library(withr)

# Basic functionality test using actual file

test_that("write_lna with file=NULL writes stub structure and read_lna works", {
  result <- write_lna(x = 1, file = NULL, transforms = character(0))
  expect_s3_class(result, "lna_write_result")
  expect_true(file.exists(result$file))
  expect_true(inherits(result$plan, "Plan"))

  handle <- read_lna(result$file)
  expect_true(inherits(handle, "DataHandle"))
  expect_false(handle$h5$is_valid())
})

# Parameter forwarding for write_lna

test_that("write_lna forwards arguments to core_write and materialise_plan", {
  captured <- list()
  fake_plan <- Plan$new()
  fake_handle <- DataHandle$new()

  with_mocked_bindings(
    core_write = function(x, transforms, transform_params) {
      captured$core <- list(x = x, transforms = transforms,
                            transform_params = transform_params)
      list(handle = fake_handle, plan = fake_plan)
    },
    materialise_plan = function(h5, plan) {
      captured$mat <- list(is_h5 = inherits(h5, "H5File"), plan = plan)
    }, {
      write_lna(x = 42, file = tempfile(fileext = ".h5"),
                transforms = c("tA"),
                transform_params = list(tA = list(foo = "bar")))
    }
  )

  expect_equal(captured$core$x, 42)
  expect_equal(captured$core$transforms, c("tA"))
  expect_equal(captured$core$transform_params, list(tA = list(foo = "bar")))
  expect_true(captured$mat$is_h5)
  expect_identical(captured$mat$plan, fake_plan)
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
  handle$h5$close_all()
})
