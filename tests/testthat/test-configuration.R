test_that("lna_options path configuration works", {
  # Test default path values
  expect_equal(lna_options("paths.scans_root")[[1]], "/scans/")
  expect_equal(lna_options("paths.transforms_root")[[1]], "/transforms/")
  expect_equal(lna_options("paths.temporal_root")[[1]], "/temporal/")
  expect_equal(lna_options("paths.metadata_root")[[1]], "/metadata/")
  
  # Test setting path values
  old_val <- lna_options("paths.scans_root")[[1]]
  lna_options(paths.scans_root = "/custom/scans/")
  expect_equal(lna_options("paths.scans_root")[[1]], "/custom/scans/")
  
  # Restore original value
  lna_options(paths.scans_root = old_val)
  expect_equal(lna_options("paths.scans_root")[[1]], old_val)
})

test_that("lna_options memory configuration works", {
  # Test default memory values
  expect_equal(lna_options("memory.target_slab_bytes")[[1]], 64e6)
  
  # Test setting memory values
  old_val <- lna_options("memory.target_slab_bytes")[[1]]
  lna_options(memory.target_slab_bytes = 128e6)
  expect_equal(lna_options("memory.target_slab_bytes")[[1]], 128e6)
  
  # Restore
  lna_options(memory.target_slab_bytes = old_val)
})

test_that("lna_options quantization configuration works", {
  # Test default quant values
  expect_equal(lna_options("quant.clip_warn_pct")[[1]], 0.5)
  expect_equal(lna_options("quant.clip_abort_pct")[[1]], 5.0)
  expect_equal(lna_options("quant.snr_sample_frac")[[1]], 0.01)
  expect_equal(lna_options("quant.bits_min")[[1]], 1L)
  expect_equal(lna_options("quant.bits_max")[[1]], 16L)
  
  # Test setting multiple values at once
  old_warn <- lna_options("quant.clip_warn_pct")[[1]]
  old_abort <- lna_options("quant.clip_abort_pct")[[1]]
  
  lna_options(
    quant.clip_warn_pct = 1.0,
    quant.clip_abort_pct = 10.0
  )
  
  expect_equal(lna_options("quant.clip_warn_pct")[[1]], 1.0)
  expect_equal(lna_options("quant.clip_abort_pct")[[1]], 10.0)
  
  # Restore
  lna_options(
    quant.clip_warn_pct = old_warn,
    quant.clip_abort_pct = old_abort
  )
})

test_that("lna_options temporal configuration works", {
  # Test temporal defaults
  expect_equal(lna_options("temporal.polynomial_order")[[1]], 3L)
  expect_equal(lna_options("temporal.bspline_order")[[1]], 3L)
  expect_equal(lna_options("temporal.wavelet_type")[[1]], "db4")
  expect_equal(lna_options("temporal.dpss_nw")[[1]], 4)
})

test_that("lna_options default identifiers work", {
  # Test default identifiers
  expect_equal(lna_options("defaults.run_id")[[1]], "run-01")
  expect_equal(lna_options("defaults.input_key")[[1]], "input")
  expect_equal(lna_options("defaults.origin_label")[[1]], "global")
})

test_that("lna_options returns NULL for non-existent options", {
  expect_null(lna_options("non.existent.option")[[1]])
  expect_null(lna_options("fake.setting")[[1]])
})

test_that("lna_options returns all options when called without arguments", {
  all_opts <- lna_options()
  expect_type(all_opts, "list")
  expect_true(length(all_opts) >= 20)  # We have at least 20 options
  expect_true("paths.scans_root" %in% names(all_opts))
  expect_true("memory.target_slab_bytes" %in% names(all_opts))
  expect_true("quant.clip_warn_pct" %in% names(all_opts))
})

test_that("lna_options can retrieve multiple values", {
  vals <- lna_options("paths.scans_root", "paths.transforms_root", "memory.target_slab_bytes")
  expect_length(vals, 3)
  expect_equal(vals$paths.scans_root, "/scans/")
  expect_equal(vals$paths.transforms_root, "/transforms/")
  expect_equal(vals$memory.target_slab_bytes, 64e6)
})

test_that(".get_lna_constant helper function works", {
  # Test getting from constants
  expect_equal(.get_lna_constant("PATHS", "SCANS_ROOT"), "/scans/")
  expect_equal(.get_lna_constant("LIMITS", "QUANT_BITS_MIN"), 1L)
  expect_equal(.get_lna_constant("LIMITS", "QUANT_BITS_MAX"), 16L)
  
  # Test with default fallback
  expect_equal(.get_lna_constant("FAKE", "NONEXISTENT", default = "fallback"), "fallback")
  expect_null(.get_lna_constant("FAKE", "NONEXISTENT"))
  
  # Test option override
  old_val <- lna_options("paths.scans_root")[[1]]
  lna_options(paths.scans_root = "/override/")
  expect_equal(.get_lna_constant("PATHS", "SCANS_ROOT"), "/override/")
  lna_options(paths.scans_root = old_val)
})

test_that("constants are properly defined", {
  # Test that key constant lists exist
  expect_true(exists(".LNA_PATHS"))
  expect_true(exists(".LNA_PATTERNS"))
  expect_true(exists(".LNA_LIMITS"))
  expect_true(exists(".LNA_ALGORITHM_DEFAULTS"))
  expect_true(exists(".LNA_TRANSFORM_TYPES"))
  expect_true(exists(".LNA_DEFAULTS"))
  expect_true(exists(".LNA_STORAGE_TYPES"))
  
  # Test some specific values
  expect_equal(.LNA_PATHS$SCANS_ROOT, "/scans/")
  expect_equal(.LNA_PATTERNS$DEFAULT_RUN_ID, "run-01")
  expect_equal(.LNA_LIMITS$QUANT_BITS_MIN, 1L)
  expect_equal(.LNA_STORAGE_TYPES$get_dtype_for_bits(8), "uint8")
  expect_equal(.LNA_STORAGE_TYPES$get_dtype_for_bits(16), "uint16")
})