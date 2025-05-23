library(hdf5r)
message("[neuroarchive-test-setup] Setting HDF5_PLUGIN_PATH to empty string.")
original_path <- Sys.getenv("HDF5_PLUGIN_PATH", unset = "<UNSET>")
message(paste0("[neuroarchive-test-setup] Original HDF5_PLUGIN_PATH was: '", original_path, "'"))
Sys.setenv(HDF5_PLUGIN_PATH = "")
message(paste0("[neuroarchive-test-setup] HDF5_PLUGIN_PATH after set to empty: '", Sys.getenv("HDF5_PLUGIN_PATH", unset = "<UNSET>"), "'"))

# Add a check for hdf5r namespace and h5test availability for debugging
tryCatch({
  if (requireNamespace("hdf5r", quietly = TRUE)) {
    message("[neuroarchive-test-setup] hdf5r namespace is available.")
    if (exists("h5test", where = "package:hdf5r")) {
        can_core <- hdf5r::h5test(type = "core")
        message(paste0("[neuroarchive-test-setup] hdf5r::h5test(type=\"core\") result: ", can_core))
    } else {
        message("[neuroarchive-test-setup] hdf5r::h5test not found in hdf5r.")
    }
  } else {
    message("[neuroarchive-test-setup] hdf5r namespace NOT available.")
  }
}, error = function(e) {
  message(paste0("[neuroarchive-test-setup] Error during HDF5 check: ", e$message))
})

message("--- HDF5R debug --- ")
tryCatch({
  message(paste0("hdf5r::H5File class: ", class(hdf5r::H5File)))
  message(paste0("is.function(hdf5r::H5File$new): ", is.function(hdf5r::H5File$new)))
  message("ls(\"package:hdf5r\"):")
  print(ls("package:hdf5r"))
}, error = function(e) {
  message(paste0("Error during hdf5r debug prints: ", e$message))
})
message("--- End HDF5R debug --- ") 