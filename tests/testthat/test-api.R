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

test_that("write_lna plugins list is written to /plugins", {
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(x = 1, file = tmp, transforms = character(0),
            plugins = list(myplugin = list(a = 1)))
  h5 <- neuroarchive:::open_h5(tmp, mode = "r")
  expect_true(h5$exists("plugins/myplugin.json"))
  grp <- h5[["plugins"]]
  desc <- read_json_descriptor(grp, "myplugin.json")
  expect_identical(desc, list(a = 1))
  neuroarchive:::close_h5_safely(h5)
})

# Parameter forwarding for write_lna

test_that("write_lna forwards arguments to core_write and materialise_plan", {
  captured <- list()
  fake_plan <- Plan$new()
  fake_handle <- DataHandle$new()

  with_mocked_bindings(
    core_write = function(x, transforms, transform_params, mask = NULL,
                          header = NULL, plugins = NULL) {
      captured$core <- list(x = x, transforms = transforms,
                            transform_params = transform_params,
                            header = header, plugins = plugins)
      list(handle = fake_handle, plan = fake_plan)
    },
    materialise_plan = function(h5, plan, checksum = "none", header = NULL,
                                plugins = NULL) {
      captured$mat <- list(is_h5 = inherits(h5, "H5File"), plan = plan,
                           header = header, plugins = plugins)
    }, {
      write_lna(x = 42, file = tempfile(fileext = ".h5"),
                transforms = c("tA"),
                transform_params = list(tA = list(foo = "bar")),
                header = list(a = 1),
                plugins = list(p = list(val = 2)))
    }
  )

  expect_equal(captured$core$x, 42)
  expect_equal(captured$core$transforms, c("tA"))
  expect_equal(captured$core$transform_params, list(tA = list(foo = "bar")))
  expect_true(captured$mat$is_h5)
  expect_identical(captured$mat$plan, fake_plan)
  expect_equal(captured$core$header, list(a = 1))
  expect_equal(captured$mat$header, list(a = 1))
  expect_equal(captured$core$plugins, list(p = list(val = 2)))
  expect_equal(captured$mat$plugins, list(p = list(val = 2)))
})

# Parameter forwarding for read_lna

test_that("read_lna forwards arguments to core_read", {
  captured <- list()
  with_mocked_bindings(
    core_read = function(file, run_id, allow_plugins, validate, output_dtype, lazy) {
      captured$core <- list(file = file, run_id = run_id, allow_plugins = allow_plugins,
                            validate = validate, output_dtype = output_dtype,
                            lazy = lazy)
      DataHandle$new()
    }, {
      read_lna("somefile.h5", run_id = "run-*", allow_plugins = "prompt", validate = TRUE,
               output_dtype = "float64", lazy = TRUE)
    }
  )

  expect_equal(captured$core$file, "somefile.h5")
  expect_equal(captured$core$run_id, "run-*")
  expect_equal(captured$core$allow_plugins, "prompt")
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

test_that("read_lna validates file argument", {
  expect_error(read_lna(1), class = "lna_error_validation")
  expect_error(read_lna(c("a", "b")), class = "lna_error_validation")
})

test_that("write_lna writes block_table dataset", {
  tmp <- local_tempfile(fileext = ".h5")
  arr <- array(1, dim = c(1, 1, 1))
  msk <- array(TRUE, dim = c(1, 1, 1))
  bt <- data.frame(start = 1L, end = 1L)
  write_lna(x = arr, file = tmp, transforms = character(0), mask = msk,
            block_table = bt)
  h5 <- neuroarchive:::open_h5(tmp, mode = "r")
  expect_true(h5$exists("spatial/block_table"))
  val <- h5[["spatial/block_table"]]$read()
  expect_equal(val, as.matrix(bt))
  neuroarchive:::close_h5_safely(h5)
})

test_that("write_lna validates block_table ranges", {
  arr <- array(1, dim = c(1, 1, 1))
  msk <- array(TRUE, dim = c(1, 1, 1))
  bt_bad <- data.frame(idx = 2L)
  expect_error(
    write_lna(x = arr, file = tempfile(fileext = ".h5"),
              transforms = character(0), mask = msk, block_table = bt_bad),
    class = "lna_error_validation"
  )
})
