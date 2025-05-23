library(testthat)
library(hdf5r)
library(withr)
library(neuroarchive)

check_roundtrip <- function(arr, axis = -1, coding = "none") {
  tmp <- local_tempfile(fileext = ".h5")
  write_lna(arr, file = tmp, transforms = "delta",
            transform_params = list(delta = list(axis = axis,
                                               coding_method = coding)))
  h <- read_lna(tmp)
  h # Return the handle itself
}

# dims[axis] == 1 should yield empty deltas and reconstruct from first_vals

test_that("delta handles axis dimension of length 1", {
  arr <- matrix(1:5, nrow = 1)

  for (coding in c("none", "rle")) {
    # message(paste0("\n[test-debug] Testing 'delta handles axis dimension of length 1' with coding: ", coding))
    tmp <- local_tempfile(fileext = ".h5")
    res <- write_lna(arr, file = tmp, transforms = "delta",
                     transform_params = list(delta = list(axis = 1,
                                                        coding_method = coding)))
    h5 <- H5File$new(tmp, mode = "r")
    dset <- h5[["/scans/run-01/deltas/00_delta/delta_stream"]]
    ds <- dset$read()
    dset$close()
    h5$close_all()
    expect_equal(nrow(ds), 0)

    h <- read_lna(tmp) # h here is the DataHandle
    
    # message("[test-debug] ---- Details for h$stash$input ----")
    # message(paste0("[test-debug] Class(h$stash$input): ", class(h$stash$input)))
    # message(paste0("[test-debug] Dim(h$stash$input): ", paste(dim(h$stash$input), collapse=",")))
    # message("[test-debug] dput(h$stash$input):")
    # print(capture.output(dput(h$stash$input)))
    # 
    # message("[test-debug] ---- Details for arr ----")
    # message(paste0("[test-debug] Class(arr): ", class(arr)))
    # message(paste0("[test-debug] Dim(arr): ", paste(dim(arr), collapse=",")))
    # message("[test-debug] dput(arr):")
    # print(capture.output(dput(arr)))
    # 
    # comparison_result <- all.equal(h$stash$input, arr)
    # message(paste0("[test-debug] all.equal(h$stash$input, arr): ", comparison_result))
    
    expect_equal(h$stash$input, arr)
  }
})

# purely 1D input roundtrips

test_that("delta handles 1D input", {
  vec <- 1:5
  for (coding in c("none", "rle")) {
    h <- check_roundtrip(vec, axis = 1, coding = coding)
    expect_equal(h$stash$input, vec)
  }
})

# single effective series in multi-dim array

test_that("delta handles prod(dims[-axis]) == 1", {
  arr <- array(1:5, dim = c(5, 1, 1))
  for (coding in c("none", "rle")) {
    h <- check_roundtrip(arr, axis = 1, coding = coding)
    expect_equal(h$stash$input, arr)
  }
})

# constant along the differencing axis

test_that("delta handles constant input along diff axis", {
  arr <- matrix(5, nrow = 4, ncol = 3)
  for (coding in c("none", "rle")) {
    h <- check_roundtrip(arr, axis = 1, coding = coding)
    expect_equal(h$stash$input, arr)
  }
})

# very short series (2 timepoints)

test_that("delta handles very short series", {
  arr <- matrix(c(1, 2, 1, 2), nrow = 2, ncol = 2)
  for (coding in c("none", "rle")) {
    h <- check_roundtrip(arr, axis = 1, coding = coding)
    expect_equal(h$stash$input, arr)
  }
})

# Removed the simplified (0x5 only) test, uncommented original
test_that("delta handles axis dimension of length 0", {
  for (coding in c("none", "rle")) {
    arr <- matrix(numeric(0), nrow = 0, ncol = 5)
    h <- check_roundtrip(arr, axis = 1, coding = coding)
    expect_equal(h$stash$input, arr)

    arr <- matrix(numeric(0), nrow = 5, ncol = 0)
    # For 5x0 matrix, axis=1 means diffing along rows. There are 5 rows, each of length 0.
    # p$orig_dims = c(5,0). p$axis=1. dims[p$axis] = 5. dim_xp_col should be 0 (prod of empty set? or 1 if length(dims[-p$axis])==0)
    # Let's test with axis = 2 (diff along columns) as well to see behavior for 0 columns if axis=1 fails.
    # However, the current setup should handle 0 columns correctly due to dim_xp_col logic. Defaulting to axis = 1.
    h <- check_roundtrip(arr, axis = 1, coding = coding) 
    expect_equal(h$stash$input, arr) 
  }
})

