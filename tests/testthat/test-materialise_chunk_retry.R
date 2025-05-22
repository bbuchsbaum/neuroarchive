library(testthat)
library(hdf5r)
library(withr)

# Simulate chunk size failures to test retry heuristics

test_that("materialise_plan retries with chunk heuristics", {
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- neuroarchive:::open_h5(tmp, mode = "w")
  plan <- Plan$new()
  plan$add_payload("p", matrix(1:10, nrow = 2))
  plan$add_dataset_def("/scans/run-01/data", "data", "dummy", "run-01", 0L, "{}", "p", "eager")

  calls <- list()
  fake_write <- function(h5_group, path, data, chunk_dims = NULL, compression_level = 0) {
    calls[[length(calls) + 1]] <<- chunk_dims
    if (length(calls) < 3) {
      stop("chunk too large")
    }
    parts <- strsplit(path, "/")[[1]]
    parts <- parts[nzchar(parts)]
    grp <- h5_group
    if (length(parts) > 1) {
      for (g in parts[-length(parts)]) {
        grp <- if (!grp$exists(g)) grp$create_group(g) else grp[[g]]
      }
    }
    grp$create_dataset(tail(parts, 1), data)
  }

  local_mocked_bindings(
    h5_write_dataset = fake_write,
    .env = asNamespace("neuroarchive")
  )

  expect_warning(
    expect_warning(
      materialise_plan(h5, plan),
      "<1 GiB"
    ),
    "256 MiB"
  )
  expect_equal(length(calls), 3)
  neuroarchive:::close_h5_safely(h5)
})

