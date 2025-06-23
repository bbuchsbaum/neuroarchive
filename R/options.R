#' Package Options for LNA
#'
#' Provides a lightweight mechanism for storing global package defaults.
#' Options are kept in an internal environment and can be retrieved or
#' updated via this helper.  Typical options include
#' `write.compression_level`, `write.chunk_target_mib` and per-transform
#' defaults such as `quant` or `delta` lists.
#'
#' @param ... Named options to set, or character names of options to
#'   retrieve.  If no arguments are provided, the full option list is
#'   returned.
#' @return A list of current options or the requested subset.  When setting
#'   values the updated option list is returned invisibly.
#' @export
lna_options <- function(...) {
  .lna_opts <- get(".lna_opts", envir = lna_options_env)
  args <- list(...)
  if (length(args) == 0) {
    return(as.list(.lna_opts))
  }
  if (is.null(names(args))) {
    return(mget(unlist(args), envir = .lna_opts, ifnotfound = list(NULL)))
  }
  for (nm in names(args)) {
    assign(nm, args[[nm]], envir = .lna_opts)
  }
  invisible(as.list(.lna_opts))
}

lna_options_env <- new.env(parent = emptyenv())
default_opts <- list(
  # I/O Configuration
  write.compression_level = 0L,
  write.chunk_target_mib = 1,
  write.file_line_limit = 50L,
  read.file_line_limit = 1000L,
  read.strict_mask_hash_validation = FALSE,
  
  # Memory and Performance
  memory.target_slab_bytes = 64e6,  # 64 MB
  
  # Quantization
  quant.clip_warn_pct = 0.5,
  quant.clip_abort_pct = 5.0,
  quant.snr_sample_frac = 0.01,
  quant.bits_min = 1L,
  quant.bits_max = 16L,
  quant = list(),
  
  # Delta Transform
  delta = list(),
  
  # Basis Transforms
  `basis.pca` = list(),
  basis.default_keep_energy = 0.99,
  
  # Temporal Basis
  temporal.polynomial_order = 3L,
  temporal.bspline_order = 3L,
  temporal.wavelet_type = "db4",
  temporal.dpss_nw = 4,
  
  # Poisson Sampling
  poisson.radius_factor = 2.5,
  poisson.density_factor = 1.5,
  
  # Edge Detection
  edge.thresh_k = 3.0,
  
  # Numeric Tolerances
  numeric.tolerance = 1e-7,
  
  # Path Configuration
  paths.scans_root = "/scans/",
  paths.transforms_root = "/transforms/",
  paths.temporal_root = "/temporal/",
  paths.metadata_root = "/metadata/",
  
  # Default Identifiers
  defaults.run_id = "run-01",
  defaults.input_key = "input",
  defaults.origin_label = "global"
)
assign(".lna_opts", list2env(default_opts, parent = emptyenv()),
       envir = lna_options_env)
