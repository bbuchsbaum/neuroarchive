% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/temporal_modwt.R
\name{.apply_modwt_threshold}
\alias{.apply_modwt_threshold}
\title{Apply MODWT-specific thresholding to wavelet coefficients}
\usage{
.apply_modwt_threshold(
  w,
  sigma_hat = NULL,
  method = c("bayes_shrink", "sure", "hard", "soft", "fdr", "none"),
  alpha = 0.05
)
}
\description{
Apply MODWT-specific thresholding to wavelet coefficients
}
\keyword{internal}
