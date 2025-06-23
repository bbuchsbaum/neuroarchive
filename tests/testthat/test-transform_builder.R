test_that("TransformBuilder creates and builds transforms correctly", {
  # Create mock plan with simple tracking
  datasets_added <- list()
  descriptors_added <- list()
  payloads_added <- list()
  
  plan <- list(
    next_index = 1,
    get_next_filename = function(type) "00_test.json",
    add_descriptor = function(fname, desc) {
      descriptors_added[[fname]] <<- desc
    },
    add_payload = function(path, data) {
      payloads_added[[path]] <<- data
    },
    add_dataset_def = function(...) {
      args <- list(...)
      datasets_added[[args[[1]]]] <<- args
    }
  )
  
  stash_update_capture <- NULL
  
  handle <- list(
    plan = plan,
    current_run_id = "run-01",
    meta = list(),
    update_stash = function(keys, values) {
      stash_update_capture <<- list(keys = keys, values = values)
      handle
    }
  )
  
  desc <- list(
    params = list(test_param = "value")
  )
  
  # Create builder
  builder <- TransformBuilder$new(handle, "test_transform", desc)
  
  # Test basic properties
  expect_equal(builder$get_run_id(), "run-01")
  expect_equal(builder$get_base_name(), "00_test")
  
  # Add datasets
  builder$add_dataset("/scans/run-01/data", c(1, 2, 3), "test_data", "float32")
  builder$add_standard_dataset("output", c(4, 5, 6), dtype = "int32")
  
  # Add report
  builder$add_report(list(status = "success", value = 42))
  
  # Set IO
  builder$set_io(inputs = "input", outputs = "output")
  
  # Build
  result <- builder$build(stash_keys = "old_key", 
                         stash_values = list(new_key = "new_value"))
  
  # Check descriptor was added
  expect_length(descriptors_added, 1)
  added_desc <- descriptors_added[["00_test.json"]]
  expect_equal(added_desc$type, "test_transform")
  expect_equal(added_desc$version, "1.0")
  expect_equal(added_desc$inputs, "input")
  expect_equal(added_desc$outputs, "output")
  
  # Check datasets were added
  expect_length(datasets_added, 3)  # 2 datasets + 1 report
  
  # Check stash update
  expect_equal(stash_update_capture$keys, "old_key")
  expect_equal(stash_update_capture$values, list(new_key = "new_value"))
})

test_that("TransformBuilder generates correct standard paths", {
  # Mock objects
  captured_calls <- list()
  
  plan <- list(
    next_index = 1,
    get_next_filename = function(type) "00_test.json",
    add_descriptor = function(...) {},
    add_payload = function(...) {},
    add_dataset_def = function(...) {
      args <- list(...)
      captured_calls <<- append(captured_calls, list(args))
    }
  )
  
  handle <- list(
    plan = plan,
    current_run_id = "run-01",
    update_stash = function(...) handle
  )
  
  builder <- TransformBuilder$new(handle, "quant", list())
  
  builder$add_standard_dataset("quantized", array(1:10), dtype = "uint8")
  builder$build()
  
  # Check the path was generated correctly
  expect_length(captured_calls, 1)
  expect_match(captured_calls[[1]][[1]], "/scans/run-01/00_test/quantized")
})

test_that("TransformBuilder handles reports correctly", {
  # Mock objects
  payloads_added <- list()
  
  plan <- list(
    next_index = 1,
    get_next_filename = function(type) "00_test.json",
    add_descriptor = function(...) {},
    add_payload = function(path, data) {
      payloads_added[[path]] <<- data
    },
    add_dataset_def = function(...) {}
  )
  
  handle <- list(
    plan = plan,
    current_run_id = "run-01",
    update_stash = function(...) handle
  )
  
  builder <- TransformBuilder$new(handle, "test", list())
  
  # Add report
  report_data <- list(
    version = "1.0",
    metrics = list(accuracy = 0.95),
    status = "complete"
  )
  
  builder$add_report(report_data)
  builder$build()
  
  # Check report was gzipped and added
  expect_length(payloads_added, 1)
  report_path <- names(payloads_added)[1]
  expect_match(report_path, "_report.json$")
  
  # Check it's gzipped
  report_payload <- payloads_added[[report_path]]
  expect_type(report_payload, "raw")
  
  # Decompress and verify
  decompressed <- memDecompress(report_payload, type = "gzip")
  json_str <- rawToChar(decompressed)
  recovered <- jsonlite::fromJSON(json_str)
  
  expect_equal(recovered$version, "1.0")
  expect_equal(recovered$metrics$accuracy, 0.95)
})

test_that("TransformBuilder parameter methods work", {
  plan <- list(
    next_index = 1,
    get_next_filename = function(type) "00_test.json",
    add_descriptor = function(...) {},
    add_payload = function(...) {},
    add_dataset_def = function(...) {}
  )
  
  handle <- list(
    plan = plan,
    current_run_id = "run-01",
    update_stash = function(...) handle
  )
  
  builder <- TransformBuilder$new(handle, "test", list())
  
  # Test getting params
  expect_null(builder$params())
  
  # Test setting params
  new_params <- list(a = 1, b = "test")
  builder$params(new_params)
  expect_equal(builder$params(), new_params)
  
  # Test chaining
  result <- builder$params(list(c = 3))$set_version("2.0")
  expect_s3_class(result, "TransformBuilder")
})