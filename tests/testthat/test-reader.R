library(testthat)
library(withr)

# Helper to create simple LNA file with no transforms
create_empty_lna <- function(path) {
  h5 <- neuroarchive:::open_h5(path, mode = "w")
  h5$create_group("transforms")
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
