library(testthat)
library(hdf5r)
library(withr)

create_valid_lna <- function(path, checksum = TRUE) {
  h5 <- neuroarchive:::open_h5(path, mode = "w")
  plan <- Plan$new()
  plan$add_descriptor("00_dummy.json", list(type = "dummy"))
  plan$add_payload("payload", matrix(1:4, nrow = 2))
  plan$add_dataset_def("/scans/run-01/data", "data", "dummy", "run-01", 0L, "{}",
                      "payload", "eager")
  materialise_plan(h5, plan, checksum = if (checksum) "sha256" else "none")
}



#' validate_lna works in a forked worker when the schema cache is cleared
#'
#' This test uses the future package with multicore plan. The worker clears the
#' internal schema cache before calling validate_lna. We expect validation to
#' succeed and return TRUE.

skip_if_not_installed("future")
skip_on_cran()

test_that("validate_lna works with schema_cache_clear in forked worker", {
  tmp <- local_tempfile(fileext = ".h5")
  create_valid_lna(tmp)
  future::plan("multicore")
  on.exit(future::plan("sequential"))
  fut <- future::future({
    library(neuroarchive)
    schema_cache_clear()
    validate_lna(tmp, checksum = FALSE)
  })
  expect_true(future::value(fut))
})
