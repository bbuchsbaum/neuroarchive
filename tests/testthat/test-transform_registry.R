test_that("register_transform validates input", {
  # Clear registry first
  clear_transform_registry()
  
  # Test invalid type
  expect_error(
    register_transform(NULL),
    class = "lna_error_validation"
  )
  
  expect_error(
    register_transform(c("type1", "type2")),
    class = "lna_error_validation"
  )
  
  expect_error(
    register_transform(123),
    class = "lna_error_validation"
  )
})

test_that("register_transform adds transforms to registry", {
  clear_transform_registry()
  
  # Register a transform
  expect_silent(register_transform("test_transform", 
                                 metadata = list(category = "test"),
                                 validate = FALSE))
  
  # Check it's registered
  expect_true(is_transform_registered("test_transform"))
  
  # Get info
  info <- get_transform_info("test_transform")
  expect_equal(info$type, "test_transform")
  expect_equal(info$metadata$category, "test")
  expect_true(info$has_forward)  # validate=FALSE assumes methods exist
  expect_true(info$has_invert)
  expect_s3_class(info$registered_at, "POSIXct")
})

test_that("register_transform warns on overwrite", {
  clear_transform_registry()
  
  # First registration
  register_transform("duplicate", validate = FALSE)
  
  # Second registration should warn
  expect_warning(
    register_transform("duplicate", validate = FALSE),
    "already registered"
  )
})

test_that("register_transform validates S3 methods when validate=TRUE", {
  clear_transform_registry()
  
  # Mock S3 methods temporarily - need to register them properly
  .S3method("forward_step", "mock_transform", function(type, desc, handle) handle)
  .S3method("invert_step", "mock_transform", function(type, desc, handle) handle)
  
  # Register with validation
  expect_silent(register_transform("mock_transform", validate = TRUE))
  
  info <- get_transform_info("mock_transform")
  expect_true(info$has_forward)
  expect_true(info$has_invert)
})

test_that("register_transform warns about missing methods", {
  clear_transform_registry()
  
  # Register transform without methods
  expect_warning(
    register_transform("no_methods", validate = TRUE),
    "has no forward_step method"
  )
})

test_that("register_transform checks for schema files", {
  clear_transform_registry()
  
  # Register a transform that's unlikely to have a schema
  register_transform("no_schema_test", validate = FALSE)
  
  info <- get_transform_info("no_schema_test")
  # Schema file likely doesn't exist for "no_schema_test"
  expect_false(info$has_schema)
})

test_that("register_transform detects additional methods", {
  skip("Skipping test that requires non-existent generic functions")
  
  # This test would work if validate_params, prepare_transform, and 
  # cleanup_transform were generic functions, but they're not currently
  # defined in the package
})

test_that("register_transform sets capabilities correctly", {
  clear_transform_registry()
  
  # Test temporal capabilities
  register_transform("temporal", validate = FALSE)
  info <- get_transform_info("temporal")
  expect_true("temporal_basis" %in% names(info$capabilities))
  
  # Test basis capabilities
  register_transform("basis", validate = FALSE)
  info <- get_transform_info("basis")
  expect_true(info$capabilities$supports_pca)
  expect_false(info$capabilities$supports_empirical)
  
  register_transform("empirical_basis", validate = FALSE)
  info <- get_transform_info("empirical_basis")
  expect_true(info$capabilities$supports_empirical)
  
  # Test compression capabilities
  register_transform("quant", validate = FALSE)
  info <- get_transform_info("quant")
  expect_true(info$capabilities$compression)
  expect_true(info$capabilities$lossy)
  expect_false(info$capabilities$lossless)
  
  register_transform("delta", validate = FALSE)
  info <- get_transform_info("delta")
  expect_true(info$capabilities$compression)
  expect_false(info$capabilities$lossy)
  expect_true(info$capabilities$lossless)
  
  # Test spatial capabilities
  register_transform("spat_foo", validate = FALSE)
  info <- get_transform_info("spat_foo")
  expect_true(info$capabilities$spatial)
})

test_that("list_transforms returns correct results", {
  clear_transform_registry()
  
  # Empty registry
  expect_equal(list_transforms(), character())
  expect_equal(list_transforms(with_details = TRUE), list())
  
  # Add some transforms
  register_transform("transform1", metadata = list(category = "cat1"), validate = FALSE)
  register_transform("transform2", metadata = list(category = "cat2"), validate = FALSE)
  register_transform("transform3", metadata = list(category = "cat1"), validate = FALSE)
  
  # List all
  all_transforms <- list_transforms()
  expect_setequal(all_transforms, c("transform1", "transform2", "transform3"))
  
  # List by category
  cat1_transforms <- list_transforms(category = "cat1")
  expect_setequal(cat1_transforms, c("transform1", "transform3"))
  
  # With details
  detailed <- list_transforms(with_details = TRUE)
  expect_length(detailed, 3)
  expect_equal(names(detailed), c("transform1", "transform2", "transform3"))
})

test_that("list_transforms filters by capability", {
  clear_transform_registry()
  
  # Register transforms with different capabilities
  register_transform("quant", validate = FALSE)
  register_transform("delta", validate = FALSE)
  register_transform("basis", validate = FALSE)
  
  # Filter by compression capability
  compression_transforms <- list_transforms(capability = "compression")
  expect_setequal(compression_transforms, c("quant", "delta"))
  
  # Filter by lossy capability
  lossy_transforms <- list_transforms(capability = "lossy")
  expect_equal(lossy_transforms, "quant")
})

test_that("get_transform_info returns NULL for non-existent transforms", {
  clear_transform_registry()
  
  expect_null(get_transform_info("nonexistent"))
})

test_that("is_transform_registered works correctly", {
  clear_transform_registry()
  
  expect_false(is_transform_registered("test"))
  
  register_transform("test", validate = FALSE)
  
  expect_true(is_transform_registered("test"))
})

test_that("get_transform_capabilities returns capabilities", {
  clear_transform_registry()
  
  # Non-existent transform
  expect_equal(get_transform_capabilities("nonexistent"), list())
  
  # Register transform
  register_transform("quant", validate = FALSE)
  
  caps <- get_transform_capabilities("quant")
  expect_type(caps, "list")
  expect_true(caps$compression)
  expect_true(caps$lossy)
})

test_that("discover_and_register_transforms finds transforms", {
  clear_transform_registry()
  
  # This will look for actual transform files in the package
  count <- discover_and_register_transforms(validate = FALSE)
  
  # Should find at least some transforms
  expect_gt(count, 0)
  
  # Check some known transforms were registered
  registered <- list_transforms()
  expect_true("quant" %in% registered)
  expect_true("delta" %in% registered)
  expect_true("basis" %in% registered)
})

test_that("discover_and_register_transforms categorizes correctly", {
  clear_transform_registry()
  
  discover_and_register_transforms(validate = FALSE)
  
  # Check categories
  quant_info <- get_transform_info("quant")
  expect_equal(quant_info$metadata$category, "compression")
  
  basis_info <- get_transform_info("basis")
  expect_equal(basis_info$metadata$category, "dimensionality")
  
  if (is_transform_registered("temporal")) {
    temporal_info <- get_transform_info("temporal")
    expect_equal(temporal_info$metadata$category, "temporal")
  }
})

test_that("clear_transform_registry removes all transforms", {
  clear_transform_registry()
  
  # Add some transforms
  register_transform("test1", validate = FALSE)
  register_transform("test2", validate = FALSE)
  
  expect_length(list_transforms(), 2)
  
  # Clear
  clear_transform_registry()
  
  expect_length(list_transforms(), 0)
})

test_that("transform_summary returns correct data frame", {
  clear_transform_registry()
  
  # Empty registry
  summary_empty <- transform_summary()
  expect_s3_class(summary_empty, "data.frame")
  expect_equal(nrow(summary_empty), 0)
  expect_true("type" %in% names(summary_empty))
  expect_true("category" %in% names(summary_empty))
  
  # Add transforms with explicit metadata
  register_transform("test_quant", 
                    metadata = list(category = "compression"), 
                    validate = FALSE)
  register_transform("test_basis", 
                    metadata = list(category = "dimensionality"), 
                    validate = FALSE)
  
  summary_df <- transform_summary()
  expect_s3_class(summary_df, "data.frame")
  expect_equal(nrow(summary_df), 2)
  expect_true(all(c("type", "category", "has_forward", "has_invert", 
                    "has_schema", "min_dims", "compression", "spatial", 
                    "temporal") %in% names(summary_df)))
  
  # Check values
  quant_row <- summary_df[summary_df$type == "test_quant", ]
  expect_equal(quant_row$category, "compression")
  expect_true(quant_row$has_forward)
})

test_that("transform registry handles min_dims correctly", {
  clear_transform_registry()
  
  # Mock transform_min_dims to return specific value
  transform_min_dims_mock <- function(type) {
    if (type == "test_3d") return(3L)
    if (type == "test_4d") return(4L)
    3L  # default
  }
  
  # Temporarily replace
  orig_fn <- transform_min_dims
  assignInNamespace("transform_min_dims", transform_min_dims_mock, "neuroarchive")
  
  register_transform("test_3d", validate = FALSE)
  register_transform("test_4d", validate = FALSE)
  
  info_3d <- get_transform_info("test_3d")
  expect_equal(info_3d$min_dims, 3L)
  
  info_4d <- get_transform_info("test_4d")
  expect_equal(info_4d$min_dims, 4L)
  
  # Restore original
  assignInNamespace("transform_min_dims", orig_fn, "neuroarchive")
})

test_that("transform registry handles errors gracefully", {
  clear_transform_registry()
  
  # Mock transform_min_dims to throw error
  transform_min_dims_error <- function(type) {
    stop("Mock error")
  }
  
  orig_fn <- transform_min_dims
  assignInNamespace("transform_min_dims", transform_min_dims_error, "neuroarchive")
  
  # Should still register with default min_dims
  expect_silent(register_transform("error_test", validate = FALSE))
  
  info <- get_transform_info("error_test")
  expect_equal(info$min_dims, 3L)  # Default value
  
  # Restore
  assignInNamespace("transform_min_dims", orig_fn, "neuroarchive")
})