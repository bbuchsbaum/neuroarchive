library(withr)
library(hdf5r)
library(neuroarchive)

# Test 1: Basic local_tempfile behavior
test_basic_tempfile <- function() {
  cat("Test 1: Basic local_tempfile behavior\n")
  
  # Create temp file
  tmp <- local_tempfile(fileext = ".h5")
  cat("Created temp file:", tmp, "\n")
  cat("File exists after creation:", file.exists(tmp), "\n")
  
  # Write something to it
  h5 <- H5File$new(tmp, mode = "w")
  h5$create_group("test")
  h5$close_all()
  
  cat("File exists after writing:", file.exists(tmp), "\n")
  
  # Try to read it
  h5_read <- H5File$new(tmp, mode = "r")
  cat("Successfully opened for reading\n")
  h5_read$close_all()
  
  cat("File exists after reading:", file.exists(tmp), "\n")
}

# Test 2: Simulate the test scenario
test_lna_scenario <- function() {
  cat("\nTest 2: LNA scenario\n")
  
  tmp <- local_tempfile(fileext = ".h5")
  cat("Created temp file:", tmp, "\n")
  cat("File exists after creation:", file.exists(tmp), "\n")
  
  # Write using write_lna
  arr <- matrix(1:5, nrow = 1)
  res <- write_lna(arr, file = tmp, transforms = "delta",
                   transform_params = list(delta = list(axis = 1,
                                                      coding_method = "none")))
  
  cat("File exists after write_lna:", file.exists(tmp), "\n")
  
  # Try to read it back
  tryCatch({
    h <- read_lna(tmp)
    cat("Successfully read with read_lna\n")
  }, error = function(e) {
    cat("Error reading with read_lna:", conditionMessage(e), "\n")
  })
}

# Test 3: Check if it's a scoping issue
test_scoping <- function() {
  cat("\nTest 3: Scoping test\n")
  
  # Store the temp file path outside
  tmp_path <- NULL
  
  local({
    tmp <- local_tempfile(fileext = ".h5")
    tmp_path <<- tmp
    cat("Inside local block - temp file:", tmp, "\n")
    cat("Inside local block - file exists:", file.exists(tmp), "\n")
  })
  
  cat("Outside local block - file exists:", file.exists(tmp_path), "\n")
}

# Run tests
test_basic_tempfile()
test_lna_scenario()
test_scoping()