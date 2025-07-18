% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/experimental_api.R
\name{temporal_project_denoise}
\alias{temporal_project_denoise}
\title{Project data onto temporal basis with optional denoising}
\usage{
temporal_project_denoise(
  data,
  basis_type = c("dct", "polynomial", "bspline", "dpss", "wavelet", "modwt"),
  n_components = NULL,
  threshold_type = c("none", "energy", "hard", "soft", "adaptive", "bayes_shrink",
    "sure", "fdr"),
  threshold_value = NULL,
  keep_energy = 0.99,
  return_basis = FALSE,
  ...
)
}
\arguments{
\item{data}{Matrix with time in rows, variables (voxels) in columns}

\item{basis_type}{Type of temporal basis: "dct", "polynomial", "bspline", "dpss", "wavelet"}

\item{n_components}{Number of basis components (NULL for automatic selection)}

\item{threshold_type}{Type of coefficient thresholding: "none", "energy", "hard", "soft", "adaptive"}

\item{threshold_value}{Threshold value (automatic if NULL for adaptive methods)}

\item{keep_energy}{For energy thresholding: fraction of energy to preserve (0-1)}

\item{return_basis}{Logical; return the basis matrix for inspection}

\item{...}{Additional parameters passed to basis generation}
}
\value{
List containing:
\itemize{
  \item \code{coefficients}: Projected coefficients (n_components x n_variables)
  \item \code{basis}: Temporal basis matrix (if \code{return_basis = TRUE})
  \item \code{compression_ratio}: Actual compression achieved
  \item \code{denoising_stats}: Statistics on thresholding effects (if applied)
}
}
\description{
High-level function that combines temporal basis generation and projection
with optional coefficient thresholding for denoising. This provides a
complete workflow for temporal dimensionality reduction.
}
\examples{
\donttest{
# Create example fMRI time series (200 TRs, 1000 voxels)
set.seed(123)
fmri_data <- matrix(rnorm(200 * 1000), nrow = 200, ncol = 1000)
# Add some temporal structure
trend <- outer(1:200, rep(1, 1000)) * 0.01
fmri_data <- fmri_data + trend

# Basic DCT compression without thresholding
result_basic <- temporal_project_denoise(fmri_data, "dct", n_components = 50)

# DCT with energy-based denoising (preserve 95\% of energy)
result_denoised <- temporal_project_denoise(
  fmri_data, "dct", 
  n_components = 50,
  threshold_type = "energy", 
  keep_energy = 0.95
)

# DPSS for resting-state fMRI with optimal parameters
dpss_params <- suggest_dpss_fmri(TR = 2.0, n_time = 200, study_type = "resting")
result_dpss <- temporal_project_denoise(
  fmri_data, "dpss",
  n_components = dpss_params$n_basis,
  time_bandwidth_product = dpss_params$time_bandwidth_product,
  return_basis = TRUE
)

# Wavelet denoising with adaptive thresholding
# (requires power-of-2 length, so subsample)
fmri_subset <- fmri_data[1:128, ]
result_wavelet <- temporal_project_denoise(
  fmri_subset, "wavelet",
  wavelet = "db4",
  threshold_type = "adaptive"
)

# Reconstruct denoised data
if (result_denoised$compression_ratio > 1) {
  reconstructed <- result_dpss$basis \%*\% result_dpss$coefficients
  cat("Reconstruction error:", mean((fmri_data - reconstructed)^2), "\n")
}
}

}
