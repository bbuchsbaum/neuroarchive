#' LNA Constants
#'
#' @description
#' Central location for package-wide constants, magic numbers, and default values.
#' These values are used throughout the package and can be referenced instead of
#' hardcoding values in multiple locations.
#'
#' @keywords internal

# HDF5 Path Constants
.LNA_PATHS <- list(
  SCANS_ROOT = "/scans/",
  TRANSFORMS_ROOT = "/transforms/",
  TEMPORAL_ROOT = "/temporal/",
  METADATA_ROOT = "/metadata/",
  DATA_SUBPATH = "data",
  VALUES_DATASET = "values",
  QUANTIZED_DATASET = "quantized",
  QUANT_SCALE_DATASET = "quant_scale",
  QUANT_OFFSET_DATASET = "quant_offset"
)

# File Patterns and Extensions
.LNA_PATTERNS <- list(
  TRANSFORM_DESC_PATTERN = "%02d_%s.json",  # e.g., "00_quant.json"
  REPORT_SUFFIX = "_report.json",
  RUN_ID_PATTERN = "^run-[0-9]+$",
  DEFAULT_RUN_ID = "run-01",
  HDF5_EXTENSIONS = c(".h5", ".hdf5"),
  NIFTI_EXTENSIONS = c(".nii", ".nii.gz"),
  SCHEMA_EXTENSION = ".schema.json"
)

# Numeric Limits and Thresholds
.LNA_LIMITS <- list(
  # Quantization
  QUANT_BITS_MIN = 1L,
  QUANT_BITS_MAX = 16L,
  QUANT_CLIP_WARN_PCT = 0.5,
  QUANT_CLIP_ABORT_PCT = 5.0,
  QUANT_SNR_SAMPLE_FRAC = 0.01,
  
  # Memory and Performance
  TARGET_SLAB_BYTES = 64e6,  # 64 MB
  CHUNK_TARGET_MIB = 1,      # 1 MiB
  FILE_READ_LINE_LIMIT = 1000L,
  FILE_WRITE_LINE_LIMIT = 50L,
  
  # Numeric Tolerances
  NUMERIC_TOLERANCE = 1e-7,
  
  # Compression
  DEFAULT_COMPRESSION_LEVEL = 0L,
  MAX_COMPRESSION_LEVEL = 9L
)

# Algorithm-Specific Defaults
.LNA_ALGORITHM_DEFAULTS <- list(
  # Poisson Sampling
  POISSON_RADIUS_FACTOR = 2.5,
  POISSON_DENSITY_FACTOR = 1.5,
  
  # Edge Detection
  EDGE_THRESH_K = 3.0,
  
  # Energy Preservation
  DEFAULT_KEEP_ENERGY = 0.99,
  
  # Temporal Basis Defaults
  TEMPORAL_BASIS = list(
    DCT = list(
      default_n_basis = function(n_timepoints) min(20L, n_timepoints %/% 5L)
    ),
    POLYNOMIAL = list(
      default_order = 3L,
      default_n_basis = function(n_timepoints) min(5L, n_timepoints %/% 10L)
    ),
    BSPLINE = list(
      default_order = 3L,
      default_n_basis = function(n_timepoints) min(n_timepoints %/% 5L, 30L)
    ),
    WAVELET = list(
      default_type = "db4",
      default_n_basis = function(n_timepoints) min(n_timepoints %/% 4L, 50L)
    ),
    DPSS = list(
      default_nw = 4
    )
  ),
  
  # fMRI Frequency Ranges (Hz)
  FMRI_FREQ_RANGES = list(
    resting_state = c(0, 0.08),
    task_based = c(0, 0.12),
    event_related = c(0, 0.16)
  )
)

# Transform Type Identifiers
.LNA_TRANSFORM_TYPES <- list(
  QUANT = "quant",
  BASIS = "basis",
  TEMPORAL = "temporal",
  DELTA = "delta",
  EMBED = "embed",
  AGGREGATE_RUNS = "aggregate_runs",
  SPARSEPCA = "sparsepca",
  HRBF = "hrbf",
  HRBF_PROJECT = "hrbf_project"
)

# Default Keys and Identifiers
.LNA_DEFAULTS <- list(
  INPUT_KEY = "input",
  OUTPUT_KEY = "output",
  AGGREGATED_MATRIX_KEY = "aggregated_matrix",
  GLOBAL_ORIGIN_LABEL = "global"
)

# Storage Type Mappings
.LNA_STORAGE_TYPES <- list(
  UINT8_MAX_BITS = 8L,
  UINT16_MAX_BITS = 16L,
  get_dtype_for_bits = function(bits) {
    if (bits <= 8L) "uint8" else "uint16"
  }
)

# Helper function to get a constant value with optional override from options
.get_lna_constant <- function(category, name, default = NULL) {
  # Try to get from options first (allows runtime override)
  opt_key <- paste(tolower(category), tolower(name), sep = ".")
  opt_val <- lna_options(opt_key)[[1]]
  if (!is.null(opt_val)) {
    return(opt_val)
  }
  
  # Otherwise get from constants
  const_name <- paste0(".LNA_", toupper(category))
  if (exists(const_name)) {
    const_list <- get(const_name)
    if (!is.null(const_list[[name]])) {
      return(const_list[[name]])
    }
  }
  
  # Return default if provided
  default
}