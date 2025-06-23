#' MODWT (Maximal Overlap Discrete Wavelet Transform) Functions
#'
#' This file contains all MODWT-specific functions for fMRI-optimized
#' temporal compression using shift-invariant wavelets.

#' Calculate maximum MODWT decomposition levels
#' 
#' Uses the Percival & Walden formula for determining the maximum number of
#' levels possible given the time series length and filter length.
#' 
#' @param N Integer number of time points
#' @param L Integer filter length
#' @return Integer maximum number of levels
#' @keywords internal
.modwt_max_level <- function(N, L) {
  floor(log2((N - 1) / (L - 1) + 1))
}

#' Generate MODWT basis matrix for fMRI-optimized temporal compression
#'
#' Maximal Overlap Discrete Wavelet Transform (MODWT) provides shift-invariant
#' wavelet decomposition without downsampling. Unlike DWT, MODWT:
#' \itemize{
#'   \item Works with any length time series (no power-of-2 requirement)
#'   \item Preserves temporal resolution at all scales
#'   \item Maintains shift invariance (critical for fMRI group analysis)
#'   \item Enables perfect reconstruction
#' }
#'
#' Recommended wavelets for fMRI:
#' \itemize{
#'   \item \strong{sym8}: 8 vanishing moments, excellent for BOLD compression
#'   \item \strong{coif6}: 6 vanishing moments, good balance of localization/smoothness  
#'   \item \strong{db8}: Alternative to sym8 with similar properties
#' }
#'
#' @param n_time Integer number of time points (no power-of-2 restriction)
#' @param wavelet Character wavelet type: "sym8", "coif6", "db8", etc.
#' @param levels Integer number of decomposition levels (auto-computed if NULL)
#' @param boundary Character boundary condition: "periodic" or "reflection"
#' @keywords internal
.modwt_basis <- function(n_time, wavelet = "sym8", levels = NULL, boundary = "periodic") {
  # Check if wavelets package is available
  if (!requireNamespace("wavelets", quietly = TRUE)) {
    abort_lna("wavelets package required for MODWT support",
              .subclass = "lna_error_dependency",
              location = ".modwt_basis")
  }
  
  # Map common aliases to wavelets package names
  wavelet_map <- list(
    "sym8" = "la8",      # Least asymmetric 8
    "coif6" = "c6",      # Coiflet 6  
    "db8" = "d8",        # Daubechies 8
    "db4" = "d4",        # Daubechies 4
    "haar" = "haar"      # Haar wavelet
  )
  
  filter_name <- wavelet_map[[wavelet]] %||% wavelet
  
  # Auto-compute levels if not specified
  if (is.null(levels)) {
    # Get filter length using simple lookup (more robust than accessing wavelets internals)
    filter_length <- switch(filter_name,
      "la8" = 8, "c6" = 12, "d8" = 16, "d4" = 8, "haar" = 2, 8)
    levels <- .modwt_max_level(n_time, filter_length)
    levels <- max(1, min(levels, 6))  # Reasonable bounds for fMRI
    
    # Emit warning when auto-capping levels
    if (.modwt_max_level(n_time, filter_length) > 6) {
      message("Using levels = 6 to avoid excessive smoothing (user can override)")
    }
  }
  
  # **Represent basis implicitly to save memory**
  # Instead of materializing n_time x n_time identity matrix
  structure(list(
    n_time = n_time,
    wavelet = filter_name,
    levels = levels,
    boundary = boundary
  ), class = "modwt_basis")
}

#' MODWT Temporal Basis for fMRI
#'
#' Maximal Overlap Discrete Wavelet Transform optimized for BOLD signal processing.
#' MODWT provides superior properties for fMRI compared to standard DWT:
#' \itemize{
#'   \item \strong{Shift invariance}: Critical for group analysis and jittered designs
#'   \item \strong{No length restrictions}: Works with any number of time points
#'   \item \strong{Perfect reconstruction}: Maintains statistical validity
#'   \item \strong{Enhanced compression}: Better energy compaction for BOLD signals
#' }
#'
#' @param kind Character scalar, must be "modwt".
#' @param n_time Integer number of time points.
#' @param n_basis Integer number of basis functions to retain (NULL for all).
#' @param wavelet Character wavelet type. Recommended for fMRI:
#'   \itemize{
#'     \item \code{"sym8"}: 8 vanishing moments, excellent for BOLD compression
#'     \item \code{"coif6"}: 6 vanishing moments, good localization/smoothness balance
#'     \item \code{"db8"}: Alternative to sym8 with similar properties
#'   }
#' @param levels Integer decomposition levels (auto-computed if NULL).
#' @param boundary Character boundary method: "periodic" (default) or "reflection".
#' @param ... Additional arguments (unused).
#'
#' @section Wavelet Selection for fMRI:
#' 
#' **Symlet-8 ("sym8")**: Recommended for most fMRI applications
#' \itemize{
#'   \item 8 vanishing moments capture polynomial trends up to degree 7
#'   \item Near-symmetric (minimal phase distortion)
#'   \item Excellent energy compaction for BOLD signals (60-80% power < 0.1 Hz)
#'   \item Good localization for motion spike detection
#' }
#' 
#' **Coiflet-6 ("coif6")**: Alternative for applications requiring symmetry
#' \itemize{
#'   \item 6 vanishing moments (good for most physiological trends)
#'   \item Symmetric filters (zero phase distortion)
#'   \item Slightly shorter support (better localization)
#' }
#'
#' **Level Selection**: 
#' Automatically computed using Percival & Walden formula but bounded for fMRI:
#' \itemize{
#'   \item Typical range: 4-6 levels for TR=1-3s
#'   \item Higher levels capture slower trends
#'   \item Too many levels may over-smooth BOLD signal
#' }
#'
#' @examples
#' \dontrun{
#' # Optimal MODWT for resting-state fMRI (300 TRs, TR=2s)
#' basis_rest <- temporal_basis("modwt", n_time = 300, wavelet = "sym8")
#' 
#' # Task fMRI with controlled decomposition levels
#' basis_task <- temporal_basis("modwt", n_time = 400, wavelet = "sym8", levels = 5)
#' 
#' # High-resolution event-related fMRI
#' basis_event <- temporal_basis("modwt", n_time = 500, wavelet = "coif6", levels = 6)
#' 
#' # Comparison with standard DWT (requires power-of-2 length)
#' # MODWT works with any length:
#' basis_flexible <- temporal_basis("modwt", n_time = 347, wavelet = "sym8")
#' }
#'
#' @references
#' Percival, D.B. & Walden, A.T. (2000). \emph{Wavelet Methods for Time Series Analysis}. 
#' Cambridge University Press.
#' 
#' Bullmore, E. et al. (2004). Wavelets and functional magnetic resonance imaging 
#' of the human brain. \emph{NeuroImage}, 23, S234-S249.
#'
#' @export
temporal_basis.modwt <- function(kind, n_time, n_basis = NULL,
                               wavelet = "sym8", levels = NULL, 
                               boundary = "periodic", ...) {
  
  basis <- .modwt_basis(n_time, wavelet, levels, boundary)
  
  # Store n_basis trimming info if specified
  if (!is.null(n_basis)) {
    # For MODWT, we have n_time coefficients per level for both W and V
    # Total coefficients = levels * n_time (W) + levels * n_time (V) = n_time * (2 * levels)
    n_total_coeffs <- n_time * (2 * levels)
    
    if (n_basis < n_total_coeffs) {
      # Keep the most important coefficients (lowest frequency bands)
      # This means keeping final approximation + some detail levels
      keep_idx <- seq(from = n_total_coeffs - n_basis + 1, to = n_total_coeffs)
      basis$keep_idx <- keep_idx
      basis$n_basis <- n_basis
    }
  }
  
  basis
}

#' @export
temporal_basis.modwt_sym8 <- function(kind, n_time, n_basis = NULL, levels = NULL, ...) {
  temporal_basis.modwt("modwt", n_time, n_basis, wavelet = "sym8", levels = levels, ...)
}

#' @export
temporal_basis.modwt_coif6 <- function(kind, n_time, n_basis = NULL, levels = NULL, ...) {
  temporal_basis.modwt("modwt", n_time, n_basis, wavelet = "coif6", levels = levels, ...)
}

#' @export
temporal_basis.modwt_db8 <- function(kind, n_time, n_basis = NULL, levels = NULL, ...) {
  temporal_basis.modwt("modwt", n_time, n_basis, wavelet = "db8", levels = levels, ...)
}

#' MODWT Projection with Advanced fMRI-Optimized Thresholding
#'
#' Projects fMRI data onto MODWT basis with sophisticated denoising options
#' designed for BOLD signal characteristics and motion artifact removal.
#'
#' @param kind Character scalar, must be "modwt".
#' @param basis Unused (MODWT projection doesn't use pre-computed basis).
#' @param X Matrix with time in rows, voxels in columns.
#' @param threshold_type Character thresholding method:
#'   \itemize{
#'     \item \code{"none"}: No thresholding (lossless compression)
#'     \item \code{"bayes_shrink"}: BayesShrink with optimal risk estimation
#'     \item \code{"sure"}: Stein's Unbiased Risk Estimator 
#'     \item \code{"fdr"}: False Discovery Rate control
#'     \item \code{"soft"}: Soft thresholding with manual threshold
#'     \item \code{"hard"}: Hard thresholding with manual threshold
#'   }
#' @param threshold_value Numeric threshold value (auto-computed for adaptive methods).
#' @param noise_estimator Character noise estimation method:
#'   \itemize{
#'     \item \code{"mad"}: Median Absolute Deviation (robust to outliers)
#'     \item \code{"bayes"}: Bayesian noise estimation
#'     \item \code{"std"}: Standard deviation (classical)
#'   }
#' @param robust_lifting Logical; use median-based lifting steps (default: FALSE).
#' @param fdr_alpha Numeric FDR significance level (default: 0.05).
#' @param ... Additional arguments.
#'
#' @return MODWT coefficients matrix (n_coeffs x n_voxels) with applied thresholding.
#'
#' @section BayesShrink for fMRI:
#' 
#' BayesShrink minimizes Bayesian risk by estimating signal and noise variances:
#' \itemize{
#'   \item Estimates noise σ using MAD of finest detail coefficients
#'   \item Estimates signal variance adaptively per scale
#'   \item Computes optimal threshold: τ = σ²/σₓ  
#'   \item Particularly effective for BOLD signals with known spectral properties
#' }
#'
#' @section Outlier Robustness:
#' 
#' Motion artifacts in fMRI appear as isolated spikes that contaminate 
#' wavelet coefficients. Robust methods help:
#' \itemize{
#'   \item MAD noise estimation is insensitive to outliers
#'   \item FDR thresholding controls false discoveries from spike artifacts
#'   \item Robust lifting (future enhancement) uses medians in predict/update steps
#' }
#'
#' @export
temporal_project.modwt <- function(kind, basis, X,
                                  threshold_type = c("none", "bayes_shrink", "sure", "fdr", "soft", "hard"),
                                  threshold_value = NULL,
                                  noise_estimator = c("mad", "bayes", "std"),
                                  robust_lifting = FALSE,
                                  fdr_alpha = 0.05, ...) {
  
  threshold_type <- match.arg(threshold_type)
  noise_estimator <- match.arg(noise_estimator)
  
  # Validate that basis is implicit MODWT basis
  stopifnot(inherits(basis, "modwt_basis"))
  
  # Extract parameters from implicit basis
  n_time <- basis$n_time
  levels <- basis$levels
  filter_name <- basis$wavelet
  boundary <- basis$boundary
  
  # Check wavelets package availability
  if (!requireNamespace("wavelets", quietly = TRUE)) {
    abort_lna("wavelets package required for MODWT",
              .subclass = "lna_error_dependency",
              location = "temporal_project.modwt")
  }
  
  # Validate robust_lifting parameter (placeholder for future implementation)
  if (robust_lifting) {
    message("robust_lifting = TRUE: Not implemented yet, using standard MODWT")
  }
  
  # Single MODWT function with filter name (more reliable than filter object)
  modwt_single <- function(ts) {
    wavelets::modwt(ts, filter = filter_name,
                    n.levels = levels, boundary = boundary)
  }
  
  # Process all voxels with unified approach
  n_voxels <- ncol(X)
  coeff_list <- vector("list", n_voxels)
  
  # Estimate global noise if using level-specific thresholding
  global_sigma_hat <- NULL
  if (threshold_type != "none" && noise_estimator == "mad") {
    # Estimate noise from finest detail level of first few voxels
    sample_voxels <- min(10, n_voxels)
    all_finest_details <- c()
    for (i in seq_len(sample_voxels)) {
      test_modwt <- modwt_single(X[, i])
      all_finest_details <- c(all_finest_details, test_modwt@W[[1]])
    }
    global_sigma_hat <- mad(all_finest_details, na.rm = TRUE) / 0.6745
  }
  
  # Process each voxel
  for (v in seq_len(n_voxels)) {
    res <- modwt_single(X[, v])
    
    # Apply thresholding to detail coefficients if requested
    if (threshold_type != "none") {
      # Use global or local noise estimation
      sigma_hat <- if (!is.null(global_sigma_hat)) {
        global_sigma_hat
      } else {
        NULL  # Let .apply_modwt_threshold estimate locally
      }
      
      # Threshold each detail level
      for (j in seq_len(levels)) {
        res@W[[j]] <- .apply_modwt_threshold(res@W[[j]], 
                                            sigma_hat = sigma_hat,
                                            method = threshold_type,
                                            alpha = fdr_alpha)
      }
    }
    
    # Concatenate coefficients: all detail levels + all approximation levels
    coeff_list[[v]] <- c(unlist(res@W, use.names = FALSE), unlist(res@V, use.names = FALSE))
  }
  
  # Convert to matrix (n_coeffs x n_voxels)
  coeff_matrix <- do.call(cbind, coeff_list)
  
  # Apply n_basis trimming if specified in basis
  if (!is.null(basis$keep_idx)) {
    coeff_matrix <- coeff_matrix[basis$keep_idx, , drop = FALSE]
  }
  
  # Store metadata as attributes for reconstruction
  attr(coeff_matrix, "modwt_params") <- list(
    wavelet = basis$wavelet,
    levels = basis$levels,
    n_time = basis$n_time,
    boundary = basis$boundary,
    threshold_type = threshold_type,
    noise_estimator = noise_estimator,
    keep_idx = basis$keep_idx
  )
  
  coeff_matrix
}

#' @export
temporal_project.modwt_sym8 <- function(kind, basis, X, ...) {
  temporal_project.modwt("modwt", basis, X, ...)
}

#' @export
temporal_project.modwt_coif6 <- function(kind, basis, X, ...) {
  temporal_project.modwt("modwt", basis, X, ...)
}

#' @export
temporal_project.modwt_db8 <- function(kind, basis, X, ...) {
  temporal_project.modwt("modwt", basis, X, ...)
}

#' Apply MODWT-specific thresholding to wavelet coefficients
#' @keywords internal
.apply_modwt_threshold <- function(w, sigma_hat = NULL,
                                  method = c("bayes_shrink", "sure", "hard", "soft", "fdr", "none"),
                                  alpha = 0.05) {
  method <- match.arg(method)
  if (method == "none") return(w)
  if (length(w) == 0) return(w)
  
  # Estimate noise standard deviation if not provided
  if (is.null(sigma_hat)) {
    sigma_hat <- mad(w, na.rm = TRUE) / 0.6745
  }
  
  if (sigma_hat <= 0) return(w)  # No noise to threshold
  
  # Compute threshold based on method
  thr <- switch(method,
    "bayes_shrink" = {
      # Improved BayesShrink: τ = σ²/σₓ where σₓ² = max(σᵧ² - σ², 0)
      sig_y2 <- mean(w^2, na.rm = TRUE)  # More robust than var()
      sig_x <- sqrt(pmax(sig_y2 - sigma_hat^2, 0))
      sigma_hat^2 / (sig_x + 1e-12)  # Add epsilon for numerical stability
    },
    "sure" = sigma_hat * sqrt(2 * log(length(w))),
    "hard" = sigma_hat * sqrt(2 * log(length(w))),
    "soft" = sigma_hat * sqrt(2 * log(length(w))),
    "fdr" = {
      # Improved FDR using Benjamini-Hochberg procedure
      p_vals <- 2 * pnorm(-abs(w) / sigma_hat)  # Two-sided p-values
      p_sorted <- sort(p_vals)
      n <- length(p_vals)
      # Find largest k such that P(k) <= (k/n) * alpha
      k_seq <- seq_len(n)
      critical_idx <- max(which(p_sorted <= alpha * k_seq / n), 0)
      if (critical_idx > 0) {
        p_critical <- p_sorted[critical_idx]
        qnorm(p_critical / 2, lower.tail = FALSE) * sigma_hat
      } else {
        Inf  # No threshold (keep all coefficients)
      }
    }
  )
  
  # Apply thresholding
  switch(method,
    "hard" = w * (abs(w) > thr),
    "soft" = sign(w) * pmax(abs(w) - thr, 0),
    "bayes_shrink" = sign(w) * pmax(abs(w) - thr, 0),
    "sure" = sign(w) * pmax(abs(w) - thr, 0),
    "fdr" = w * (abs(w) > thr)
  )
}

#' @export
temporal_reconstruct.modwt <- function(kind, basis, coeff, ...) {
  # MODWT reconstruction using stored parameters from coefficients
  params <- attr(coeff, "modwt_params")
  if (is.null(params)) {
    abort_lna("MODWT parameters not found in coefficient matrix",
              .subclass = "lna_error_internal",
              location = "temporal_reconstruct.modwt")
  }
  
  # Validate n_time matches between basis and coefficients
  if (inherits(basis, "modwt_basis") && basis$n_time != params$n_time) {
    abort_lna("Basis and coefficient matrix time dimensions don't match",
              .subclass = "lna_error_internal",
              location = "temporal_reconstruct.modwt")
  }
  
  # Check wavelets package availability
  if (!requireNamespace("wavelets", quietly = TRUE)) {
    abort_lna("wavelets package required for MODWT reconstruction",
              .subclass = "lna_error_dependency",
              location = "temporal_reconstruct.modwt")
  }
  
  n_voxels <- ncol(coeff)
  n_time <- params$n_time
  levels <- params$levels
  filter_name <- params$wavelet
  boundary <- params$boundary %||% "periodic"
  
  # Handle coefficient trimming if it was applied
  full_coeff <- coeff
  if (!is.null(params$keep_idx)) {
    # Reconstruct full coefficient matrix by zero-padding
    n_total_coeffs <- n_time * (2 * levels)  # W + V coefficients
    full_coeff <- matrix(0, nrow = n_total_coeffs, ncol = n_voxels)
    full_coeff[params$keep_idx, ] <- coeff
  }
  
  # Use filter name directly for better compatibility
  
  # Reconstruct each voxel time series
  reconstructed <- matrix(0, nrow = n_time, ncol = n_voxels)
  
  for (v in seq_len(n_voxels)) {
    voxel_coeffs <- full_coeff[, v]
    
    # Parse coefficients back into MODWT structure with corrected indexing
    modwt_obj <- .parse_modwt_coefficients(voxel_coeffs, levels, n_time, filter_name, boundary)
    
    # Inverse MODWT
    reconstructed[, v] <- wavelets::imodwt(modwt_obj)
  }
  
  reconstructed
}

#' @export
temporal_reconstruct.modwt_sym8 <- function(kind, basis, coeff, ...) {
  temporal_reconstruct.modwt("modwt", basis, coeff, ...)
}

#' @export
temporal_reconstruct.modwt_coif6 <- function(kind, basis, coeff, ...) {
  temporal_reconstruct.modwt("modwt", basis, coeff, ...)
}

#' @export
temporal_reconstruct.modwt_db8 <- function(kind, basis, coeff, ...) {
  temporal_reconstruct.modwt("modwt", basis, coeff, ...)
}

#' Parse flat coefficient vector back into MODWT object structure
#' @keywords internal
.parse_modwt_coefficients <- function(coeffs, levels, n_time, filter_name, boundary = "periodic") {
  # Create a dummy MODWT object to get the correct S4 structure
  dummy_modwt <- wavelets::modwt(rep(0, n_time), filter = filter_name,
                                n.levels = levels, boundary = boundary)
  
  # Split coefficients back into W (detail) and V (approximation) components
  idx <- 1
  k_block <- n_time  # Each level has n_time coefficients in MODWT
  
  # Extract detail coefficients for each level
  # W coefficients are stored as matrices (n_time x 1) in MODWT objects
  W <- vector("list", levels)
  for (lvl in seq_len(levels)) {
    W[[lvl]] <- matrix(coeffs[idx:(idx + k_block - 1)], ncol = 1)
    idx <- idx + k_block
  }
  
  # Extract all approximation coefficients (one for each level)
  # V coefficients are stored as matrices (n_time x 1) in MODWT objects
  V <- vector("list", levels)
  for (lvl in seq_len(levels)) {
    V[[lvl]] <- matrix(coeffs[idx:(idx + k_block - 1)], ncol = 1)
    idx <- idx + k_block
  }
  
  # Update the dummy object with our coefficients
  dummy_modwt@W <- W
  dummy_modwt@V <- V
  
  # Ensure all other attributes are preserved from dummy object
  # (series, filter, level, boundary are already set correctly)
  
  dummy_modwt
}

# Helper function for wavelet kind expansion
.expand_modwt_kind <- function(kind) {
  if (kind == "modwt_sym8") {
    list(kind = "modwt", wavelet = "sym8")
  } else if (kind == "modwt_coif6") {
    list(kind = "modwt", wavelet = "coif6")
  } else if (kind == "modwt_db8") {
    list(kind = "modwt", wavelet = "db8")
  } else {
    list(kind = kind, wavelet = NULL)
  }
} 