library(testthat)
library(withr)

# Helper to create simple LNA file with no transforms
create_empty_lna <- function(path) {
  h5 <- neuroarchive:::open_h5(path, mode = "w")
  h5$create_group("transforms")
  neuroarchive:::close_h5_safely(h5)
}

# Helper to create lna with one dummy descriptor
create_dummy_lna <- function(path) {
  h5 <- neuroarchive:::open_h5(path, mode = "w")
  tf <- h5$create_group("transforms")
  write_json_descriptor(tf, "00_dummy.json", list(type = "dummy"))
  neuroarchive:::close_h5_safely(h5)
}


test_that("read_lna(lazy=TRUE) returns lna_reader", {
  tmp <- local_tempfile(fileext = ".h5")
  create_empty_lna(tmp)

  reader <- read_lna(tmp, lazy = TRUE)
  expect_s3_class(reader, "lna_reader")
  expect_true(reader$h5$is_valid())
  reader$close()
})

test_that("lna_reader initialize closes file on failure", {
  tmp <- local_tempfile(fileext = ".h5")
  create_empty_lna(tmp)
  orig_open_h5 <- getFromNamespace("open_h5", "neuroarchive")
  captured_h5 <- NULL
  expect_error(
    with_mocked_bindings(
      open_h5 = function(file, mode = "r") {
        captured_h5 <<- orig_open_h5(file, mode)
        captured_h5
      },
      read_lna(tmp, run_id = "run-01", lazy = TRUE)
    ),
    class = "lna_error_run_id"
  )
  expect_true(inherits(captured_h5, "H5File"))
  expect_false(captured_h5$is_valid())
})


test_that("lna_reader close is idempotent", {
  tmp <- local_tempfile(fileext = ".h5")
  create_empty_lna(tmp)

  reader <- read_lna(tmp, lazy = TRUE)
  expect_true(reader$h5$is_valid())
  reader$close()
  expect_null(reader$h5)
  expect_silent(reader$close())
})


test_that("lna_reader data caches result and respects subset", {
  tmp <- local_tempfile(fileext = ".h5")
  create_empty_lna(tmp)

  reader <- read_lna(tmp, lazy = TRUE)
  h1 <- reader$data()
  expect_s3_class(h1, "DataHandle")
  h2 <- reader$data()
  expect_identical(h1, h2)

  reader$subset(roi_mask = 1)
  h3 <- reader$data()
  expect_false(identical(h1, h3))
  expect_identical(h3$subset$roi_mask, 1)

  h4 <- reader$data()
  expect_identical(h3, h4)

  reader$close()
})

test_that("read_lna lazy passes subset params", {
  tmp <- local_tempfile(fileext = ".h5")
  create_empty_lna(tmp)
  msk <- array(TRUE, dim = c(1,1,1))
  reader <- read_lna(tmp, lazy = TRUE, roi_mask = msk, time_idx = 2)
  expect_identical(reader$subset_params$roi_mask, msk)
  expect_identical(reader$subset_params$time_idx, 2L)
  reader$close()
})

test_that("lna_reader$data allow_plugins='none' errors on unknown transform", {
  tmp <- local_tempfile(fileext = ".h5")
  create_dummy_lna(tmp)
  reader <- read_lna(tmp, lazy = TRUE, allow_plugins = "none")
  expect_error(reader$data(), class = "lna_error_no_method")
  reader$close()
})

test_that("lna_reader$data allow_plugins='prompt' falls back when non-interactive", {
  tmp <- local_tempfile(fileext = ".h5")
  create_dummy_lna(tmp)
  reader <- read_lna(tmp, lazy = TRUE, allow_plugins = "prompt")
  with_mocked_bindings(
    rlang::is_interactive = function() FALSE,
    { expect_warning(reader$data(), "Missing invert_step") }
  )
  reader$close()
})

test_that("lna_reader$data allow_plugins='prompt' interactive respects choice", {
  tmp <- local_tempfile(fileext = ".h5")
  create_dummy_lna(tmp)
  reader <- read_lna(tmp, lazy = TRUE, allow_plugins = "prompt")
  with_mocked_bindings(
    rlang::is_interactive = function() TRUE,
    readline = function(prompt = "") "n",
    { expect_error(reader$data(), class = "lna_error_no_method") }
  )
  reader$close()

  reader <- read_lna(tmp, lazy = TRUE, allow_plugins = "prompt")
  with_mocked_bindings(
    rlang::is_interactive = function() TRUE,
    readline = function(prompt = "") "y",
    { expect_warning(reader$data(), "Missing invert_step") }
  )
  reader$close()
})
