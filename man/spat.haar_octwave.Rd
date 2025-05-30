\name{spat.haar_octwave}
\alias{spat.haar_octwave}
\title{Mask-Adaptive Haar Wavelet Pyramid Transform}
\description{
The \code{spat.haar_octwave} transform applies a Haar lifting scheme on a Morton-ordered implicit octree. Only the root and per-level detail coefficients are stored, keeping the basis implicit while preserving perfect reconstruction before quantisation.
}
\section{Parameters}{
\itemize{
  \item{\code{levels}}{Number of decomposition levels.}
  \item{\code{detail_threshold_type}}{Strategy for sparsifying detail coefficients.}
  \item{\code{detail_threshold_value}}{Threshold used when sparsifying.}
}}
\details{
A full description of the algorithm and layout is provided in \file{raw-data/Haar_Pyramind.md}. See also the vignette section \dQuote{Appendix: Mask-Adaptive Haar Octwave}.
}
