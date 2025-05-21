library(testthat)

# forward_step error provenance

forward_step.fail <- function(type, desc, handle) {
  stop("boom")
}
assign("forward_step.fail", forward_step.fail, envir = .GlobalEnv)
withr::defer(rm("forward_step.fail", envir = .GlobalEnv))

test_that("core_write reports step provenance", {
  err <- expect_error(core_write(x = array(1, dim = c(1, 1, 1)), transforms = "fail"))
  expect_true(grepl("forward_step.fail\[0\]", err$location))
})

# invert_step error provenance
invert_step.fail <- function(type, desc, handle) {
  stop("oops")
}
assign("invert_step.fail", invert_step.fail, envir = .GlobalEnv)
withr::defer(rm("invert_step.fail", envir = .GlobalEnv))

create_fail_file <- function(path) {
  h5 <- neuroarchive:::open_h5(path, mode = "w")
  tf <- h5$create_group("transforms")
  write_json_descriptor(tf, "00_fail.json", list(type = "fail"))
  neuroarchive:::close_h5_safely(h5)
}

test_that("core_read reports step provenance", {
  tmp <- local_tempfile(fileext = ".h5")
  create_fail_file(tmp)
  err <- expect_error(core_read(tmp))
  expect_true(grepl("invert_step.fail\[0\]", err$location))
})

create_double_fail_file <- function(path) {
  h5 <- neuroarchive:::open_h5(path, mode = "w")
  tf <- h5$create_group("transforms")
  write_json_descriptor(tf, "00_dummy.json", list(type = "dummy"))
  write_json_descriptor(tf, "01_fail.json", list(type = "fail"))
  neuroarchive:::close_h5_safely(h5)
}

test_that("core_read error location uses transform index", {
  tmp <- local_tempfile(fileext = ".h5")
  create_double_fail_file(tmp)
  invert_step.dummy <- function(type, desc, handle) handle
  assign("invert_step.dummy", invert_step.dummy, envir = .GlobalEnv)
  withr::defer(rm("invert_step.dummy", envir = .GlobalEnv))

  err <- expect_error(core_read(tmp))
  expect_true(grepl("invert_step.fail\[1\]", err$location))
})
