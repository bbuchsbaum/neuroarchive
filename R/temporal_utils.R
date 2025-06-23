#' Temporal Utility Functions
#'
#' This file contains utility functions for temporal transformations,
#' including parameter suggestion helpers and other convenience functions.

#' Suggest DPSS parameters for fMRI applications
#'
#' Provides recommended `time_bandwidth_product` values for DPSS temporal 
#' basis functions based on TR and study type. This is a convenience function
#' to help users select appropriate parameters without deep knowledge of 
#' spectral analysis.
#'
#' @param TR Numeric repetition time in seconds.
#' @param n_time Integer number of time points (TRs).
#' @param study_type Character scalar indicating study design:
#' \itemize{
#'   \item \code{"resting"}: Resting-state fMRI (preserves 0-0.08 Hz)
#'   \item \code{"task"}: Task-based fMRI (preserves 0-0.12 Hz)  
#'   \item \code{"event"}: Event-related fMRI (preserves 0-0.16 Hz)
#'   \item \code{"custom"}: Use `max_freq` parameter
#' }
#' @param max_freq Numeric maximum frequency to preserve (Hz). 
#'   Only used when `study_type = "custom"`.
#' @param conservative Logical. If `TRUE`, use slightly narrower bandwidth
#'   for more aggressive noise removal.
#'
#' @return A list with suggested parameters:
#' \itemize{
#'   \item \code{time_bandwidth_product}: Recommended NW value
#'   \item \code{n_basis}: Suggested number of basis functions
#'   \item \code{preserved_freq}: Approximate frequency range preserved (Hz)
#'   \item \code{notes}: Additional guidance
#' }
#'
#' @examples
#' # Resting-state fMRI with TR = 2s, 300 timepoints
#' params_rest <- suggest_dpss_fmri(TR = 2.0, n_time = 300, study_type = "resting")
#' basis <- temporal_basis("dpss", n_time = 300, 
#'                         n_basis = params_rest$n_basis,
#'                         time_bandwidth_product = params_rest$time_bandwidth_product)
#' 
#' # Task fMRI with faster sampling
#' params_task <- suggest_dpss_fmri(TR = 1.0, n_time = 400, study_type = "task")
#' 
#' # Custom frequency range
#' params_custom <- suggest_dpss_fmri(TR = 1.5, n_time = 350, 
#'                                    study_type = "custom", max_freq = 0.1)
#' 
#' @export
suggest_dpss_fmri <- function(TR, n_time, study_type = c("resting", "task", "event", "custom"),
                              max_freq = NULL, conservative = FALSE) {
  study_type <- match.arg(study_type)
  
  stopifnot(is.numeric(TR), length(TR) == 1, TR > 0)
  stopifnot(is.numeric(n_time), length(n_time) == 1, n_time > 0)
  
  nyquist_freq <- 1 / (2 * TR)
  
  # Define target frequencies by study type
  target_freq <- switch(study_type,
    "resting" = 0.08,   # Typical BOLD bandwidth
    "task" = 0.12,      # Include task frequencies
    "event" = 0.16,     # Broader for event-related
    "custom" = {
      if (is.null(max_freq)) {
        stop("max_freq must be specified when study_type = 'custom'")
      }
      max_freq
    }
  )
  
  if (target_freq >= nyquist_freq) {
    warning(sprintf("Target frequency (%.3f Hz) approaches Nyquist (%.3f Hz). Consider faster sampling.", 
                    target_freq, nyquist_freq))
  }
  
  # Calculate time-bandwidth product
  # For DPSS: bandwidth ≈ 2*NW / (n_time * TR) Hz
  # To preserve frequencies up to target_freq, we want: 2*NW/(n_time*TR) ≥ target_freq
  # Therefore: NW ≥ target_freq * n_time * TR / 2
  
  # Simplified approach: empirically tested values for fMRI
  if (study_type == "resting") {
    # Conservative: typical values that preserve 0-0.08 Hz well
    time_bandwidth_product <- if (TR >= 2.0) 2.0 else 2.5
  } else if (study_type == "task") {
    # Moderate: preserve task-related frequencies up to ~0.12 Hz
    time_bandwidth_product <- if (TR >= 2.0) 3.0 else 3.5
  } else if (study_type == "event") {
    # Broader: preserve event-related frequencies up to ~0.16 Hz
    time_bandwidth_product <- if (TR >= 1.5) 4.0 else 4.5
  } else {
    # Custom: use theoretical calculation but with practical bounds
    # For DPSS: to preserve up to f_max, need NW ≈ f_max * T_total / 2
    # where T_total = n_time * TR is the total duration
    total_duration <- n_time * TR
    time_bandwidth_product <- target_freq * total_duration / 2
    time_bandwidth_product <- pmax(1.5, pmin(time_bandwidth_product, 6.0))
  }
  
  # Apply conservative adjustment if requested
  if (conservative) {
    time_bandwidth_product <- time_bandwidth_product * 0.8
  }
  
  # Practical bounds for fMRI
  time_bandwidth_product <- pmax(1.5, pmin(time_bandwidth_product, n_time / 4))
  
  # Suggest number of basis functions (must satisfy n_basis <= 2*NW constraint)
  # Use about 80% of maximum allowed to ensure good concentration
  max_allowed_basis <- floor(2 * time_bandwidth_product)
  n_basis_suggested <- round(0.8 * max_allowed_basis)
  n_basis_suggested <- pmax(1, pmin(n_basis_suggested, max_allowed_basis, n_time / 2))
  
  # Calculate actual preserved frequency
  # For DPSS: bandwidth in Hz = 2 * NW / (n_time * TR)
  # This is the half-bandwidth; total bandwidth is twice this
  preserved_freq <- 2 * time_bandwidth_product / (n_time * TR)
  
  # Generate notes
  notes <- sprintf(
    "For %s fMRI with TR=%.1fs: preserves 0-%.3f Hz, attenuates >%.3f Hz", 
    study_type, TR, preserved_freq, preserved_freq * 1.5
  )
  
  if (time_bandwidth_product < 2) {
    notes <- paste(notes, "Very narrow bandwidth - may over-smooth.")
  } else if (time_bandwidth_product > 4) {
    notes <- paste(notes, "Broad bandwidth - less noise reduction.")
  }
  
  list(
    time_bandwidth_product = round(time_bandwidth_product, 2),
    n_basis = as.integer(n_basis_suggested),
    preserved_freq = round(preserved_freq, 4),
    notes = notes
  )
} 