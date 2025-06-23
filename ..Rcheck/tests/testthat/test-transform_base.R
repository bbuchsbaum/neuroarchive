test_that("extract_transform_params validates parameters correctly", {
  # Test successful extraction with defaults
  desc <- list(params = list(bits = 10))
  param_spec <- list(
    bits = list(default = 8, type = "integer", min = 1, max = 16),
    method = list(default = "range", type = "character", enum = c("range", "sd"))
  )
  
  result <- extract_transform_params(desc, param_spec)
  expect_equal(result$bits, 10)
  expect_equal(result$method, "range")
  
  # Test type validation
  desc_bad <- list(params = list(bits = "not_a_number"))
  expect_error(
    extract_transform_params(desc_bad, param_spec),
    class = "lna_error_validation"
  )
  
  # Test range validation
  desc_bad <- list(params = list(bits = 20))
  expect_error(
    extract_transform_params(desc_bad, param_spec),
    class = "lna_error_validation"
  )
  
  # Test enum validation
  desc_bad <- list(params = list(method = "invalid"))
  expect_error(
    extract_transform_params(desc_bad, param_spec),
    class = "lna_error_validation"
  )
})

test_that("validator functions work correctly", {
  # Test integer range validator
  validator <- validate_integer_range(1, 10)
  expect_true(validator(5, "test"))
  expect_match(validator(0, "test"), "must be between")
  expect_match(validator(11, "test"), "must be between")
  expect_match(validator(5.5, "test"), "must be an integer")
  
  # Test enum validator
  validator <- validate_enum(c("a", "b", "c"))
  expect_true(validator("b", "test"))
  expect_match(validator("d", "test"), "must be one of")
  
  # Test probability validator
  expect_true(validate_probability(0.5, "test"))
  expect_match(validate_probability(1.5, "test"), "must be between")
  expect_match(validate_probability(-0.1, "test"), "must be between")
})

test_that("generate_transform_path creates correct paths", {
  # Save current options
  old_scans <- lna_options("paths.scans_root")[[1]]
  old_transforms <- lna_options("paths.transforms_root")[[1]]
  old_temporal <- lna_options("paths.temporal_root")[[1]]
  
  # Test standard dataset path
  path <- generate_transform_path("quant", "run-01", "00_quant", "quantized")
  expect_equal(path, "/scans/run-01/00_quant/quantized")
  
  # Test temporal basis path
  path <- generate_transform_path("temporal", "run-01", "01_temporal", "basis")
  expect_equal(path, "/temporal/01_temporal/basis")
  
  # Test report path
  path <- generate_transform_path("quant", "run-01", "00_quant", "report")
  expect_equal(path, "/transforms/00_quant_report.json")
  
  # Test with custom paths
  lna_options(paths.scans_root = "/custom/scans/")
  path <- generate_transform_path("delta", "run-02", "02_delta", "delta_stream")
  expect_equal(path, "/custom/scans/run-02/02_delta/delta_stream")
  
  # Restore options
  lna_options(
    paths.scans_root = old_scans,
    paths.transforms_root = old_transforms,
    paths.temporal_root = old_temporal
  )
})

test_that("resolve_input_key finds correct input", {
  # Create mock handle with stash
  handle <- list(
    stash = list(
      input = 1,
      temporal_coefficients = 2,
      delta_stream = 3
    ),
    has_key = function(key) key %in% names(handle$stash)
  )
  
  # Test preferred keys
  desc <- list(inputs = c("other_key"))
  key <- resolve_input_key(handle, desc, c("temporal_coefficients", "delta_stream"))
  expect_equal(key, "temporal_coefficients")
  
  # Test descriptor inputs
  desc <- list(inputs = c("delta_stream"))
  key <- resolve_input_key(handle, desc, character())
  expect_equal(key, "delta_stream")
  
  # Test fallback
  desc <- list()
  key <- resolve_input_key(handle, desc, c("nonexistent"))
  expect_equal(key, "input")
})

test_that("apply_subset handles different data types", {
  # Test matrix subsetting
  mat <- matrix(1:20, nrow = 5, ncol = 4)
  subset <- list(
    roi_mask = c(TRUE, FALSE, TRUE, FALSE),
    time_idx = c(1, 3, 5)
  )
  
  result <- apply_subset(mat, subset, "matrix")
  expect_equal(dim(result), c(3, 2))
  expect_equal(result[1, 1], 1)  # First time point, first ROI
  expect_equal(result[2, 1], 3)  # Third time point, first ROI
  
  # Test with NULL subset
  result <- apply_subset(mat, NULL, "matrix")
  expect_identical(result, mat)
  
  # Test 4D array subsetting
  arr <- array(1:120, dim = c(3, 4, 5, 2))
  subset <- list(time_idx = c(1))
  
  result <- apply_subset(arr, subset, "array")
  expect_equal(dim(result), c(3, 4, 5, 1))
})

test_that("with_transform_error_handling provides context", {
  # Test validation error gets context
  err <- expect_error(
    with_transform_error_handling(
      abort_lna("test error", .subclass = "lna_error_validation"),
      type = "test_transform",
      step = "forward"
    ),
    class = "lna_error_validation"
  )
  
  expect_equal(err$transform_type, "test_transform")
  expect_equal(err$transform_step, "forward")
  
  # Test generic error gets wrapped
  err <- expect_error(
    with_transform_error_handling(
      stop("generic error"),
      type = "test_transform",
      step = "forward"
    ),
    class = "lna_error_transform"
  )
  
  expect_match(conditionMessage(err), "Transform 'test_transform' failed")
})