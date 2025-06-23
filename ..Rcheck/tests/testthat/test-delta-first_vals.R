library(testthat)
library(hdf5r)
library(withr)


# ensure first_vals dataset dimension is stored and read correctly

test_that("first_vals handled for 1D input", {
  vec <- 1:5
  tmp <- local_tempfile(fileext = ".h5")
  
  # Capture any warnings or messages during write_lna
  result <- tryCatch({
    write_lna(vec, file = tmp, transforms = "delta",
              transform_params = list(delta = list(axis = 1)))
  }, warning = function(w) {
    cat("WARNING during write_lna:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
    write_lna(vec, file = tmp, transforms = "delta",
              transform_params = list(delta = list(axis = 1)))
  }, error = function(e) {
    cat("ERROR during write_lna:", conditionMessage(e), "\n")
    stop(e)
  })
  
  cat("DEBUG: write_lna completed successfully\n")
  cat("DEBUG: Result class:", class(result), "\n")
  cat("DEBUG: Result structure:\n")
  str(result, max.level = 1)
  
  h5 <- H5File$new(tmp, mode = "r")
  
  # Get the datasets tibble and find first_values
  datasets <- result$plan$datasets
  cat("DEBUG: Number of datasets:", nrow(datasets), "\n")
  cat("DEBUG: Roles found:", paste(datasets$role, collapse = ", "), "\n")
  cat("DEBUG: Paths found:", paste(datasets$path, collapse = ", "), "\n")
  
  first_values_idx <- which(datasets$role == "first_values")
  expect_length(first_values_idx, 1)
  first_path <- datasets$path[first_values_idx[1]]
  
  fv <- h5[[first_path]]$read()
  h5$close_all()
  expect_equal(dim(matrix(fv, nrow = 1)), c(1, 1))

  h <- read_lna(tmp)
  expect_equal(drop(h$stash$input), vec)
})

test_that("first_vals handled for 2D axis=1", {
  mat <- matrix(1:12, nrow = 3, ncol = 4)
  tmp <- local_tempfile(fileext = ".h5")
  res <- write_lna(mat, file = tmp, transforms = "delta",
                   transform_params = list(delta = list(axis = 1)))
  h5 <- H5File$new(tmp, mode = "r")
  
  # Get the datasets tibble and find first_values
  datasets <- res$plan$datasets
  first_values_idx <- which(datasets$role == "first_values")
  expect_length(first_values_idx, 1)
  first_path <- datasets$path[first_values_idx[1]]
  
  fv <- h5[[first_path]]$read()
  h5$close_all()
  expect_equal(dim(matrix(fv, nrow = 1)), c(1, ncol(mat)))

  h <- read_lna(tmp)
  expect_equal(drop(h$stash$input), mat)
})

test_that("first_vals handled for 3D axis>1", {
  arr <- array(seq_len(24), dim = c(2, 3, 4))
  tmp <- local_tempfile(fileext = ".h5")
  res <- write_lna(arr, file = tmp, transforms = "delta",
                   transform_params = list(delta = list(axis = 2)))
  h5 <- H5File$new(tmp, mode = "r")
  
  # Get the datasets tibble and find first_values
  datasets <- res$plan$datasets
  first_values_idx <- which(datasets$role == "first_values")
  expect_length(first_values_idx, 1)
  first_path <- datasets$path[first_values_idx[1]]
  
  fv <- h5[[first_path]]$read()
  h5$close_all()
  expect_equal(dim(matrix(fv, nrow = 1)), c(1, prod(dim(arr)[-2])))

  h <- read_lna(tmp)
  expect_equal(drop(h$stash$input), arr)
})
