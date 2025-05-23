library(devtools)
devtools::load_all(".")

# Let's test the sparsepca descriptor and examine what it generates
tryCatch({
  library(testthat)
  library(withr)
  
  # Simple aggregator plugin used for testing
  .forward_step.myorg.aggregate_runs <- function(type, desc, handle) {
    # Respect the pipeline-provided input key
    input_key <- desc$inputs[[1]]
    lst <- handle$stash[[input_key]]
    stopifnot(is.list(lst))
    mats <- lapply(lst, function(x) {
      if (is.matrix(x)) x else as.matrix(x)
    })
    aggregated <- do.call(rbind, mats)
    desc$version <- "1.0"
    # Keep the pipeline-provided inputs and outputs
    handle$plan$add_descriptor(handle$plan$get_next_filename(type), desc)
    
    # Use the pipeline-provided output key
    output_key <- desc$outputs[[1]]
    new_values <- setNames(list(aggregated), output_key)
    handle <- handle$update_stash(keys = input_key, new_values = new_values)
  }
  
  assign("forward_step.myorg.aggregate_runs", .forward_step.myorg.aggregate_runs, envir = .GlobalEnv)
  
  testthat::local_mocked_bindings(
    default_params = function(type) {
      if (type == "myorg.aggregate_runs") return(list())
      if (type == "myorg.sparsepca") return(list(k=10))
      if (type == "delta") return(list())
      neuroarchive:::default_params(type)
    },
    .package = "neuroarchive"
  )
  
  set.seed(1)
  run1_data <- matrix(rnorm(50), nrow = 10, ncol = 5)
  dim(run1_data) <- c(dim(run1_data), 1)
  run2_data <- matrix(rnorm(50), nrow = 10, ncol = 5)
  dim(run2_data) <- c(dim(run2_data), 1)
  
  # Create a test file with both transforms
  tmp <- tempfile(fileext = ".h5")
  cat(sprintf("Writing to: %s\n", tmp))
  
  write_lna(list(`run-01` = run1_data, `run-02` = run2_data), 
            file = tmp,
            transforms = c("myorg.aggregate_runs", "myorg.sparsepca"),
            transform_params = list(myorg.sparsepca = list(k = 3)))
  
  # Now examine what's in the HDF5 file
  h5_file <- hdf5r::H5File$new(tmp, mode = "r")
  tf_group <- h5_file[["transforms"]]
  
  cat("\n=== Content of 01_myorg.sparsepca.json ===\n")
  desc_dataset <- tf_group[["01_myorg.sparsepca.json"]]
  desc_content <- desc_dataset[]
  desc_dataset$close()
  cat(desc_content)
  cat("\n")
  
  h5_file$close()
  
}, error = function(e) {
  cat("Error occurred:\n")
  print(e)
  cat("\nError message:\n")
  cat(conditionMessage(e))
  cat("\n")
}) 