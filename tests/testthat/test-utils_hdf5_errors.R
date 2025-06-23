test_that("h5_attr_write validates inputs correctly", {
  temp_h5_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_h5_file))
  
  h5_file <- H5File$new(temp_h5_file, mode = "w")
  on.exit(h5_file$close_all(), add = TRUE)
  
  root_group <- h5_file[["/"]]
  
  # Test invalid h5_obj
  expect_error(
    h5_attr_write("not_an_h5_object", "attr", 123),
    class = "lna_error_validation"
  )
  
  # Test invalid attribute name
  expect_error(
    h5_attr_write(root_group, 123, "value"),
    class = "lna_error_validation"
  )
  
  expect_error(
    h5_attr_write(root_group, c("attr1", "attr2"), "value"),
    class = "lna_error_validation"
  )
})

test_that("h5_attr_read validates inputs and handles missing attributes", {
  temp_h5_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_h5_file))
  
  h5_file <- H5File$new(temp_h5_file, mode = "w")
  on.exit(h5_file$close_all(), add = TRUE)
  
  root_group <- h5_file[["/"]]
  
  # Test invalid h5_obj
  expect_error(
    h5_attr_read("not_an_h5_object", "attr"),
    class = "lna_error_validation"
  )
  
  # Test invalid attribute name
  expect_error(
    h5_attr_read(root_group, 123),
    class = "lna_error_validation"
  )
  
  # Test missing attribute
  expect_error(
    h5_attr_read(root_group, "nonexistent_attr"),
    class = "lna_error_missing_data"
  )
})

test_that("h5_create_empty_dataset validates inputs", {
  temp_h5_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_h5_file))
  
  h5_file <- H5File$new(temp_h5_file, mode = "w")
  on.exit(h5_file$close_all(), add = TRUE)
  
  root_group <- h5_file[["/"]]
  
  # Test invalid h5_group
  expect_error(
    h5_create_empty_dataset("not_a_group", "dataset", c(10, 10), "float32"),
    class = "lna_error_validation"
  )
  
  # Test invalid path
  expect_error(
    h5_create_empty_dataset(root_group, 123, c(10, 10), "float32"),
    class = "lna_error_validation"
  )
  
  expect_error(
    h5_create_empty_dataset(root_group, "", c(10, 10), "float32"),
    class = "lna_error_validation"
  )
  
  # Test invalid dimensions
  expect_error(
    h5_create_empty_dataset(root_group, "dataset", "not_numeric", "float32"),
    class = "lna_error_validation"
  )
  
  # Test invalid dtype
  expect_error(
    h5_create_empty_dataset(root_group, "dataset", c(10, 10), 123),
    class = "lna_error_validation"
  )
})

test_that("h5_write_dataset validates inputs", {
  temp_h5_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_h5_file))
  
  h5_file <- H5File$new(temp_h5_file, mode = "w")
  on.exit(h5_file$close_all(), add = TRUE)
  
  root_group <- h5_file[["/"]]
  
  # Test invalid h5_group
  expect_error(
    h5_write_dataset("not_a_group", "dataset", matrix(1:4, 2, 2)),
    class = "lna_error_validation"
  )
  
  # Test invalid path
  expect_error(
    h5_write_dataset(root_group, 123, matrix(1:4, 2, 2)),
    class = "lna_error_validation"
  )
  
  # Test invalid compression level
  expect_error(
    h5_write_dataset(root_group, "dataset", matrix(1:4, 2, 2), compression_level = "high"),
    class = "lna_error_validation"
  )
  
  # Test invalid data type
  expect_error(
    h5_write_dataset(root_group, "dataset", list(a = 1, b = 2)),
    class = "lna_error_validation"
  )
})

test_that("guess_chunk_dims handles edge cases", {
  # Test with empty dimensions
  chunks <- guess_chunk_dims(c(0, 100, 100), 4)
  expect_equal(chunks, c(1L, 1L, 1L))
  
  # Test with negative dimensions should error
  expect_error(
    guess_chunk_dims(c(-10, 100, 100), 4),
    class = "lna_error_validation"
  )
  
  # Test with empty dimension vector should error
  expect_error(
    guess_chunk_dims(numeric(0), 4),
    class = "lna_error_validation"
  )
})

test_that("reduce_chunk_dims validates inputs", {
  # Test with non-numeric chunk
  expect_error(
    reduce_chunk_dims("not_numeric", 4, 1024),
    class = "lna_error_validation"
  )
  
  # Test normal operation
  reduced <- reduce_chunk_dims(c(1000, 1000, 10), 4, 1024)
  expect_true(all(reduced <= c(1000, 1000, 10)))
  expect_true(prod(reduced) * 4 <= 1024)
})

test_that("open_h5 validates inputs and handles failures", {
  # Test invalid path type
  expect_error(
    open_h5(123),
    class = "lna_error_validation"
  )
  
  expect_error(
    open_h5(c("file1.h5", "file2.h5")),
    class = "lna_error_validation"
  )
  
  # Test invalid mode
  expect_error(
    open_h5("test.h5", mode = 123),
    class = "lna_error_validation"
  )
  
  # Test opening non-existent file in read mode
  expect_error(
    open_h5("/nonexistent/path/file.h5", mode = "r"),
    class = "lna_error_io"
  )
})

test_that("assert_h5_path validates inputs", {
  temp_h5_file <- tempfile(fileext = ".h5")
  on.exit(unlink(temp_h5_file))
  
  h5_file <- H5File$new(temp_h5_file, mode = "w")
  on.exit(h5_file$close_all(), add = TRUE)
  
  # Test invalid h5 object
  expect_error(
    assert_h5_path("not_h5", "path"),
    class = "lna_error_validation"
  )
  
  # Test invalid path
  expect_error(
    assert_h5_path(h5_file, 123),
    class = "lna_error_validation"
  )
  
  # Test missing path
  expect_error(
    assert_h5_path(h5_file, "/nonexistent/path"),
    class = "lna_error_missing_path"
  )
})

test_that("map_dtype handles unknown types", {
  expect_error(
    map_dtype("unknown_type"),
    class = "lna_error_validation"
  )
  
  expect_error(
    map_dtype(123),
    class = "lna_error_validation"
  )
})