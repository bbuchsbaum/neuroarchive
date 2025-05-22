library(testthat)
library(hdf5r)
library(withr)

check_roundtrip <- function(arr, axis = -1, coding = "none") {
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(arr, file = tmp, transforms = "delta",
            transform_params = list(delta = list(axis = axis,
                                               coding_method = coding)))
  h <- read_lna(tmp)
  drop(h$stash$input)
}

# dims[axis] == 1 should yield empty deltas and reconstruct from first_vals

test_that("delta handles axis dimension of length 1", {
  arr <- matrix(1:5, nrow = 1)
  for (coding in c("none", "rle")) {
    tmp <- local_tempfile(fileext = ".h5")
    res <- write_lna(arr, file = tmp, transforms = "delta",
                     transform_params = list(delta = list(axis = 1,
                                                        coding_method = coding)))
    h5 <- H5File$new(tmp, mode = "r")
    ds <- h5[["/scans/run-01/deltas/00_delta/delta_stream"]][]
    h5$close_all()
    expect_equal(nrow(ds), 0)
    h <- read_lna(tmp)
    expect_equal(drop(h$stash$input), arr)
  }
})

# purely 1D input roundtrips

test_that("delta handles 1D input", {
  vec <- 1:5
  for (coding in c("none", "rle")) {
    expect_equal(check_roundtrip(vec, axis = 1, coding = coding), vec)
  }
})

# single effective series in multi-dim array

test_that("delta handles prod(dims[-axis]) == 1", {
  arr <- array(1:5, dim = c(5, 1, 1))
  for (coding in c("none", "rle")) {
    expect_equal(check_roundtrip(arr, axis = 1, coding = coding), arr)
  }
})

# constant along the differencing axis

test_that("delta handles constant input along diff axis", {
  arr <- matrix(5, nrow = 4, ncol = 3)
  for (coding in c("none", "rle")) {
    expect_equal(check_roundtrip(arr, axis = 1, coding = coding), arr)
  }
})

# very short series (2 timepoints)

test_that("delta handles very short series", {
  arr <- matrix(c(1, 2, 1, 2), nrow = 2, ncol = 2)
  for (coding in c("none", "rle")) {
    expect_equal(check_roundtrip(arr, axis = 1, coding = coding), arr)
  }
})

