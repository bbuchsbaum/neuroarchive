#' Temporal Projection Functions
#'
#' This file contains all temporal projection functions that project data
#' onto temporal bases with optional thresholding/denoising.

#' Project data onto a temporal basis
#'
#' Each temporal basis kind can implement customised projection logic.
#' The default assumes an orthonormal basis and uses `crossprod`.
#' @keywords internal
temporal_project <- function(kind, basis, X, ...) {
  stopifnot(is.character(kind), length(kind) == 1)
  obj <- structure(kind, class = c(kind, "character"))
  UseMethod("temporal_project", obj)
}

#' @export
temporal_project.default <- function(kind, basis, X, ...) {
  crossprod(basis, X)
}

#' @export
temporal_project.dct <- function(kind, basis, X, 
                                threshold_type = c("none", "energy", "hard"), 
                                threshold_value = NULL,
                                keep_energy = NULL, ...) {
  threshold_type <- match.arg(threshold_type)
  coeffs <- crossprod(basis, X)
  
  if (threshold_type != "none") {
    if (threshold_type == "energy" && !is.null(keep_energy)) {
      # Energy-based thresholding: keep coefficients that preserve X% of energy
      energy_per_coeff <- rowSums(coeffs^2)
      total_energy <- sum(energy_per_coeff)
      
      if (total_energy > 0) {
        cumsum_energy <- cumsum(sort(energy_per_coeff, decreasing = TRUE))
        n_keep <- which(cumsum_energy >= keep_energy * total_energy)[1]
        
        if (!is.na(n_keep) && n_keep < length(energy_per_coeff)) {
          threshold <- sort(energy_per_coeff, decreasing = TRUE)[n_keep]
          coeffs[energy_per_coeff < threshold, ] <- 0
        }
      }
    } else if (threshold_type == "hard" && !is.null(threshold_value)) {
      # Hard thresholding
      coeffs[abs(coeffs) < threshold_value] <- 0
    }
  }
  
  coeffs
}

#' @export
temporal_project.wavelet <- function(kind, basis, X, 
                                   threshold_type = c("none", "soft", "hard", "adaptive"),
                                   threshold_value = NULL, ...) {
  threshold_type <- match.arg(threshold_type)
  coeffs <- crossprod(basis, X)
  
  if (threshold_type != "none") {
    if (is.null(threshold_value)) {
      # Adaptive threshold using Donoho-Johnstone estimator
      sigma_est <- mad(as.vector(coeffs), na.rm = TRUE) / 0.6745
      threshold_value <- sigma_est * sqrt(2 * log(length(coeffs)))
    }
    
    coeffs <- switch(threshold_type,
      "soft" = sign(coeffs) * pmax(abs(coeffs) - threshold_value, 0),
      "hard" = coeffs * (abs(coeffs) > threshold_value),
      "adaptive" = sign(coeffs) * pmax(abs(coeffs) - threshold_value, 0)
    )
  }
  
  coeffs
}

#' @export
temporal_project.bspline <- function(kind, basis, X, ...) {
  qrB <- qr(basis)
  if (qrB$rank < ncol(basis)) {
    message("[temporal_project.bspline WARN] B-spline basis is rank deficient. Projection may be unstable.")
  }
  qr.coef(qrB, X)
} 