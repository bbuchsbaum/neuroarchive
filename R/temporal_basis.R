#' Temporal Basis Generation Functions
#'
#' This file contains all functions related to generating temporal basis matrices
#' for different basis types (DCT, B-spline, DPSS, polynomial, wavelet).

#' Generate an orthonormal DCT basis matrix
#' @keywords internal
.dct_basis <- function(n_time, n_basis) {
  if (n_basis <= 0) {
    return(matrix(0, nrow = n_time, ncol = 0))
  }
  
  t <- seq_len(n_time) - 0.5
  k <- seq_len(n_basis) - 1
  B <- sqrt(2 / n_time) * cos(outer(t, k, function(ti, ki) pi * ti * ki / n_time))
  B[,1] <- B[,1] / sqrt(2)
  B
}

#' Generate a B-spline basis matrix
#' @keywords internal
.bspline_basis <- function(n_time, n_basis, order) {
  x <- seq_len(n_time)
  splines::bs(x, df = n_basis, degree = order, intercept = TRUE)
}

#' Generate DPSS basis matrix
#'
#' Uses \code{multitaper::dpss()} when available for optimal performance.
#' Falls back to eigendecomposition of the sinc kernel matrix when multitaper
#' is not available.
#' @keywords internal
.dpss_basis <- function(n_time, n_basis, NW) {
  stopifnot(is.numeric(NW), length(NW) == 1)
  NW <- as.numeric(NW)
  stopifnot(NW > 0, NW < n_time / 2)
  stopifnot(is.numeric(n_basis), length(n_basis) == 1)
  n_basis <- as.integer(n_basis)
  stopifnot(n_basis <= 2 * NW)

  if (requireNamespace("multitaper", quietly = TRUE)) {
    res <- multitaper::dpss(n = n_time, k = n_basis, nw = NW)
    V <- res$v
    V <- V[, seq_len(min(n_basis, ncol(V))), drop = FALSE]
    return(V)
  }

  # Fallback implementation
  W <- NW / n_time
  m <- as.double(seq_len(n_time) - 1)
  diff <- outer(m, m, "-")

  S <- sin(2 * pi * W * diff) / (pi * diff)
  diag(S) <- 2 * W

  eig <- eigen(S, symmetric = TRUE)
  V <- eig$vectors[, seq_len(n_basis), drop = FALSE]

  V <- sweep(V, 2, sqrt(colSums(V^2)), "/")
  for (j in seq_len(ncol(V))) if (V[1, j] < 0) V[, j] <- -V[, j]

  V
}

#' Generate orthogonal polynomial basis matrix
#' @keywords internal
.polynomial_basis <- function(n_time, n_basis) {
  if (n_basis <= 0) {
    return(matrix(0, nrow = n_time, ncol = 0))
  }
  
  # First column: constant (0-th order polynomial, scaled to norm 1)
  col_const <- matrix(1 / sqrt(n_time), nrow = n_time, ncol = 1)
  
  if (n_basis == 1) {
    return(col_const)
  } else {
    # Higher order orthogonal polynomials (degrees 1 to n_basis-1)
    # stats::poly(..., degree = k) gives k columns for degrees 1...k
    # These are orthogonal to each other and to a constant intercept.
    degree_for_poly <- n_basis - 1 
    P_ortho <- stats::poly(seq_len(n_time), degree = degree_for_poly, raw = FALSE)
    
    # Combine the constant term with the higher-order orthogonal polynomials
    # The columns of P_ortho are already orthogonal to an intercept.
    return(cbind(col_const, P_ortho))
  }
}

#' Generate wavelet basis matrix using the `wavelets` package
#'
#' Daubechies 4 ("db4") tends to balance temporal resolution and smoothness
#' for fMRI applications, so it is used as the default.
#' Any wavelet supported by the `wavelets` package may be supplied. The
#' computation is vectorised with `vapply` to avoid creating an explicit
#' identity matrix.
#' @keywords internal
.wavelet_basis <- function(n_time, wavelet = "db4") {
  if (log2(n_time) %% 1 != 0) {
    abort_lna("wavelet basis requires power-of-two length",
              .subclass = "lna_error_validation",
              location = ".wavelet_basis")
  }
  # Map common aliases (e.g., "db4") to names expected by `wavelets::dwt`
  if (is.character(wavelet)) {
    wl <- tolower(wavelet)
    if (wl == "db1") {
      wl <- "haar"
    } else if (grepl("^db[0-9]+$", wl)) {
      wl <- sub("^db", "d", wl)
    }
    wavelet <- wl
  }
  J <- log2(n_time)
  filter <- wavelets::wavelet(wavelet)
  basis <- matrix(0, n_time, n_time)
  for (k in seq_len(n_time)) {
    e_k <- numeric(n_time)
    e_k[k] <- 1
    w <- wavelets::dwt(e_k, filter = filter, n.levels = J, boundary = "periodic")
    basis[, k] <- c(unlist(w@W), w@V[[w@level]])
  }
  basis
}

#' Generate temporal basis matrix
#'
#' Dispatches on \code{kind} to create a temporal basis. Package authors can
#' extend this generic by defining methods named \code{temporal_basis.<kind>}.
#'
#' @param kind Character scalar identifying the basis type. Options include:
#' \itemize{
#'   \item \code{"dct"}: Discrete Cosine Transform (fast, general-purpose)
#'   \item \code{"polynomial"}: Orthogonal polynomials (simple detrending)
#'   \item \code{"bspline"}: B-spline basis (smooth, flexible trends)
#'   \item \code{"dpss"}: Discrete Prolate Spheroidal Sequences (optimal spectral concentration)
#'   \item \code{"wavelet"}: Wavelet basis (time-frequency localization)
#' }
#' @param n_time Integer number of time points.
#' @param n_basis Integer number of basis functions.
#' @param ... Additional arguments passed to methods.
#' @return A basis matrix with dimensions \code{n_time x n_basis}.
#' 
#' @section Choosing a temporal basis for fMRI:
#' 
#' \itemize{
#'   \item \strong{DCT}: Fast, good for general compression. Use when computational speed matters.
#'   \item \strong{Polynomial}: Simple linear/quadratic detrending. Use for basic drift removal.
#'   \item \strong{B-spline}: Smooth, flexible trends. Good for slow drift with irregular patterns.
#'   \item \strong{DPSS}: Optimal frequency concentration. Use when preserving specific BOLD frequencies is critical.
#'   \item \strong{Wavelet}: Time-frequency localization. Use for non-stationary signals.
#' }
#' 
#' For most fMRI preprocessing, \code{"dct"} or \code{"polynomial"} suffice. 
#' Use \code{"dpss"} when you need precise control over frequency content 
#' (e.g., removing physiological artifacts while preserving BOLD signal).
#'
#' @examples
#' # General-purpose DCT basis for compression
#' dct_basis <- temporal_basis("dct", n_time = 200, n_basis = 50)
#' 
#' # Simple polynomial detrending (linear + quadratic)
#' poly_basis <- temporal_basis("polynomial", n_time = 200, n_basis = 3)
#' 
#' # DPSS for preserving BOLD frequencies (0-0.08 Hz with TR=2s)
#' dpss_basis <- temporal_basis("dpss", n_time = 300, n_basis = 12, 
#'                              time_bandwidth_product = 2.0)
#' 
#' @export
temporal_basis <- function(kind, n_time, n_basis, ...) {
  stopifnot(is.character(kind), length(kind) == 1)
  obj <- structure(kind, class = c(kind, "character"))
  UseMethod("temporal_basis", obj)
}

#' @export
temporal_basis.dct <- function(kind, n_time, n_basis, ...) {
  .dct_basis(n_time, n_basis)
}

#' @export
temporal_basis.bspline <- function(kind, n_time, n_basis, order = 3, ...) {
  .bspline_basis(n_time, n_basis, order)
}

#' DPSS Temporal Basis for fMRI
#'
#' Discrete Prolate Spheroidal Sequences (DPSS) provide optimal time-frequency 
#' localization for spectral analysis. In fMRI applications, DPSS bases can 
#' preserve BOLD signal frequencies while attenuating scanner artifacts.
#'
#' @param kind Character scalar, must be "dpss".
#' @param n_time Integer number of time points (TRs).
#' @param n_basis Integer number of DPSS basis functions to retain.
#' @param time_bandwidth_product Numeric controlling frequency concentration.
#'   Lower values (1-2) provide narrower frequency bands but less flexibility.
#'   Higher values (3-4) allow broader frequency content but less concentration.
#' @param n_tapers Integer number of tapers to compute (defaults to n_basis).
#' @param ... Additional arguments (unused).
#'
#' @section Choosing time_bandwidth_product for fMRI:
#' 
#' The `time_bandwidth_product` (NW) controls the trade-off between frequency
#' resolution and concentration. For fMRI applications:
#' 
#' **For TR = 2s (0.5 Hz Nyquist frequency):**
#' \itemize{
#'   \item NW = 1.5: Preserves 0-0.05 Hz (ultra-low frequency BOLD)
#'   \item NW = 2.0: Preserves 0-0.08 Hz (typical BOLD bandwidth)  
#'   \item NW = 3.0: Preserves 0-0.12 Hz (includes some task frequencies)
#' }
#' 
#' **For TR = 1s (1.0 Hz Nyquist frequency):**
#' \itemize{
#'   \item NW = 2.0: Preserves 0-0.08 Hz (BOLD signal)
#'   \item NW = 3.0: Preserves 0-0.12 Hz (BOLD + low task frequencies)
#'   \item NW = 4.0: Preserves 0-0.16 Hz (includes more task-related signals)
#' }
#' 
#' **Rule of thumb:** Frequency bandwidth ≈ NW / (n_time * TR)
#'
#' @examples
#' \dontrun{
#' # Example 1: Conservative denoising for resting-state fMRI
#' # TR = 2s, 300 TRs (10 minutes), preserve 0-0.08 Hz
#' basis_rest <- temporal_basis("dpss", n_time = 300, n_basis = 12, 
#'                              time_bandwidth_product = 2.0)
#' 
#' # Example 2: Task fMRI with faster sampling  
#' # TR = 1s, 400 TRs, preserve 0-0.12 Hz for task frequencies
#' basis_task <- temporal_basis("dpss", n_time = 400, n_basis = 20,
#'                              time_bandwidth_product = 3.0)
#' 
#' # Example 3: High-resolution temporal filtering
#' # TR = 0.8s, 500 TRs, broader bandwidth for event-related designs
#' basis_event <- temporal_basis("dpss", n_time = 500, n_basis = 25,
#'                               time_bandwidth_product = 4.0)
#' 
#' # Example 4: Artifact removal while preserving BOLD
#' # Remove respiratory (~0.3 Hz) and cardiac (~1 Hz) with TR = 1s
#' # Use NW = 2 to concentrate energy below 0.08 Hz
#' basis_clean <- temporal_basis("dpss", n_time = 600, n_basis = 15,
#'                               time_bandwidth_product = 2.0)
#' }
#'
#' @references
#' Thomson, D.J. (1982). Spectrum estimation and harmonic analysis. 
#' \emph{Proceedings of the IEEE}, 70(9), 1055-1096.
#' 
#' Percival, D.B. & Walden, A.T. (1993). \emph{Spectral Analysis for Physical 
#' Applications}. Cambridge University Press.
#' 
#' @export
temporal_basis.dpss <- function(kind, n_time, n_basis,
                               time_bandwidth_product = 3,
                               n_tapers = n_basis, ...) {
  n_tapers <- n_tapers %||% n_basis
  stopifnot(is.numeric(n_tapers), length(n_tapers) == 1)
  n_tapers <- as.integer(n_tapers)
  n_basis <- min(n_basis, n_tapers, n_time)
  .dpss_basis(n_time, n_basis, NW = time_bandwidth_product)
}

#' @export
temporal_basis.polynomial <- function(kind, n_time, n_basis, ...) {
  .polynomial_basis(n_time, n_basis)
}

#' @export
temporal_basis.wavelet <- function(kind, n_time, n_basis, wavelet = "db4", ...) {
  basis <- .wavelet_basis(n_time, wavelet)
  if (!is.null(n_basis)) basis <- basis[, seq_len(min(n_basis, ncol(basis))), drop = FALSE]
  basis
}

#' @export
temporal_basis.default <- function(kind, n_time, n_basis, ...) {
  abort_lna(
    sprintf("Unsupported temporal kind '%s'", kind),
    .subclass = "lna_error_validation",
    location = "temporal_basis:kind"
  )
} 