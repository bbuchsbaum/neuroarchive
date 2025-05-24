#' Temporal Transform - Forward Step
#'
#' Projects data onto a temporal basis (DCT, B-spline, DPSS, polynomial, or wavelet).
#' Debug messages are controlled by the `lna.debug.temporal` option.
#' @keywords internal
#' @export
forward_step.temporal <- function(type, desc, handle) {
  dbg <- isTRUE(getOption("lna.debug.temporal", FALSE))
  p <- desc$params %||% list()
  # Extract temporal-specific parameters and remove them from p to avoid duplication
  kind <- p$kind %||% "dct"
  n_basis <- p$n_basis
  p$kind <- NULL
  p$n_basis <- NULL
  order <- p$order %||% 3
  p$order <- NULL

  if (!is.null(n_basis)) {
    if (!is.numeric(n_basis) || length(n_basis) != 1 ||
        n_basis <= 0 || n_basis %% 1 != 0) {
      abort_lna(
        "n_basis must be a positive integer",
        .subclass = "lna_error_validation",
        location = "forward_step.temporal:n_basis"
      )
    }
    n_basis <- as.integer(n_basis)
  }

  if (!is.null(order)) {
    if (!is.numeric(order) || length(order) != 1 ||
        order <= 0 || order %% 1 != 0) {
      abort_lna(
        "order must be a positive integer",
        .subclass = "lna_error_validation",
        location = "forward_step.temporal:order"
      )
    }
    order <- as.integer(order)
  }
  # Determine input key based on previous transform's output.
  # When temporal coefficients are already present we treat them as
  # the input so additional temporal steps operate on the projected
  # coefficients rather than reusing the raw matrix.
  if (handle$has_key("temporal_coefficients")) {
    input_key <- "temporal_coefficients"
  } else if (handle$has_key("delta_stream")) {
    input_key <- "delta_stream"
  } else if (handle$has_key("sparsepca_embedding")) {
    input_key <- "sparsepca_embedding"
  } else if (handle$has_key("aggregated_matrix")) {
    input_key <- "aggregated_matrix"
  } else {
    input_key <- if (!is.null(desc$inputs)) desc$inputs[[1]] else "input"
  }

  X <- handle$get_inputs(input_key)[[1]]
  X <- as_dense_mat(X)

  n_time <- nrow(X)
  if (is.null(n_basis)) n_basis <- n_time
  n_basis <- min(n_basis, n_time)

  # After resolving defaults, store parameters back in desc$params
  p_final <- c(list(kind = kind, n_basis = n_basis, order = order), p)
  desc$params <- p_final

  args <- c(list(kind = kind, n_time = n_time, n_basis = n_basis, order = order),
            p)
  
  basis <- do.call(temporal_basis, args)

  # Delegate projection logic to per-kind methods for extensibility
  coeff <- temporal_project(kind, basis, X)

  if (dbg) {
    # DEBUG: Check reconstruction locally
    if (is.matrix(basis) && is.matrix(coeff) && ncol(basis) == nrow(coeff)) {
      if (identical(kind, "polynomial")) {
        message("[forward_step.temporal POLY DEBUG] Checking orthogonality of basis (t(basis) %*% basis):")
        # Ensure it's a plain matrix for printing, and round for clarity
        t_basis_basis <- as.matrix(crossprod(basis))
        print(round(t_basis_basis, 5))
      }
      X_reconstructed_debug <- basis %*% coeff # Should be time x features
      if (!isTRUE(all.equal(X, X_reconstructed_debug, tolerance = 1e-7))) {
        message("[forward_step.temporal DEBUG] Local reconstruction MISMATCH.")
        if (identical(kind, "polynomial")) {
           message("Sum of squared differences: ", sum((X - X_reconstructed_debug)^2))
        }
      } else {
        message("[forward_step.temporal DEBUG] Local reconstruction MATCHES.")
      }
    } else {
      message("[forward_step.temporal DEBUG] Could not perform local reconstruction check due to matrix non-conformance.")
    }
    # END DEBUG
  }

  run_id <- handle$current_run_id %||% "run-01"
  run_id <- sanitize_run_id(run_id)
  plan <- handle$plan
  fname <- plan$get_next_filename(type)
  base_name <- tools::file_path_sans_ext(fname)
  basis_path <- paste0("/temporal/", base_name, "/basis")
  coef_path <- paste0("/scans/", run_id, "/", base_name, "/coefficients")
  knots_path <- paste0("/temporal/", base_name, "/knots")
  params_json <- as.character(jsonlite::toJSON(desc$params, auto_unbox = TRUE))

  desc$version <- "1.0"
  desc$inputs <- c(input_key)
  desc$outputs <- c("temporal_coefficients")
  datasets <- list(
    list(path = basis_path, role = "temporal_basis"),
    list(path = coef_path, role = "temporal_coefficients")
  )
  knots_data <- attr(basis, "knots")
  if (!is.null(knots_data)) {
    datasets[[length(datasets) + 1]] <- list(path = knots_path, role = "knots")
  }
  desc$datasets <- datasets

  plan$add_descriptor(fname, desc)
  
  basis_payload <- basis
  plan$add_payload(basis_path, basis_payload)
  


  plan$add_dataset_def(basis_path, "temporal_basis", as.character(type), run_id,
                       as.integer(plan$next_index), params_json,
                       basis_path, "eager", dtype = "float32")

  if (!is.null(knots_data)) {
    knots_payload <- knots_data
    plan$add_payload(knots_path, knots_payload)
    plan$add_dataset_def(knots_path, "knots", as.character(type), run_id,
                         as.integer(plan$next_index), params_json,
                         knots_path, "eager", dtype = "float32")
  }
  
  coeff_payload <- coeff
  plan$add_payload(coef_path, coeff_payload)

  plan$add_dataset_def(coef_path, "temporal_coefficients", as.character(type), run_id,
                       as.integer(plan$next_index), params_json,
                       coef_path, "eager", dtype = "float32")
  handle$plan <- plan

  handle$update_stash(keys = c(input_key),
                      new_values = list(temporal_coefficients = coeff))
}

#' Temporal Transform - Inverse Step
#'
#' Reconstructs data from stored temporal basis coefficients.
#' Debug messages are controlled by the `lna.debug.temporal` option.
#' @keywords internal
#' @export
invert_step.temporal <- function(type, desc, handle) {
  dbg <- isTRUE(getOption("lna.debug.temporal", FALSE))
  if (dbg) message(sprintf("[invert_step.temporal ENTRY] Incoming handle stash keys: %s. Is input NULL? %s", paste(names(handle$stash), collapse=", "), is.null(handle$stash$input)))
  basis_path <- NULL
  coeff_path <- NULL
  
  if (!is.null(desc$datasets)) {
    roles <- vapply(desc$datasets, function(d) d$role, character(1))
    idx_b <- which(roles == "temporal_basis")
    if (length(idx_b) > 0) basis_path <- desc$datasets[[idx_b[1]]]$path
    idx_c <- which(roles == "temporal_coefficients")
    if (length(idx_c) > 0) coeff_path <- desc$datasets[[idx_c[1]]]$path
  }

  if (is.null(basis_path)) {
    abort_lna(
      "temporal_basis path not found in descriptor",
      .subclass = "lna_error_descriptor",
      location = "invert_step.temporal:basis_path"
    )
  }
  if (is.null(coeff_path)) {
    abort_lna("temporal_coefficients path not found in descriptor datasets", .subclass = "lna_error_descriptor", location = "invert_step.temporal")
  }

  output_stash_key  <- desc$inputs[[1]] %||% "input"

  root <- handle$h5[["/"]]
  basis <- h5_read(root, basis_path)
  coeff <- h5_read(root, coeff_path)

  if (dbg) {
    message(sprintf("[invert_step.temporal] Basis dims: %s", paste(dim(basis), collapse = "x")))
    message(sprintf("[invert_step.temporal] Coeff dims: %s", paste(dim(coeff), collapse = "x")))
  }
  

  
  if (dbg) {
    message("--- Invert Step Pre-Dense Calculation Debug ---")
    if (nrow(basis) >= 2 && ncol(basis) >= 2) {
      message("basis_loaded[1:2, 1:2]:")
      print(basis[1:2, 1:2, drop = FALSE])
    }
    if (nrow(coeff) >= 2 && ncol(coeff) >= 2) {
      message("coeff_loaded[1:2, 1:2]:")
      print(coeff[1:2, 1:2, drop = FALSE])
    }
  }
  
  # Check for valid matrix dimensions before multiplication
  if (!is.matrix(basis) || !is.matrix(coeff)) {
    # Handle case where basis/coeff are stored dimension vectors from empty arrays
    if (length(basis) == 2 && length(coeff) == 2 && 
        all(basis >= 0) && all(coeff >= 0) && 
        all(basis == as.integer(basis)) && all(coeff == as.integer(coeff))) {
      # These look like stored dimensions - reconstruct the original empty matrices
      basis_dims <- as.integer(basis)
      coeff_dims <- as.integer(coeff)
      
      
      
      # Reconstruct the matrices
      basis <- array(numeric(0), dim = basis_dims)
      coeff <- array(numeric(0), dim = coeff_dims)
      
      # Now check if we can do matrix multiplication
      if (basis_dims[2] != coeff_dims[1]) {
        # Dimensions don't match for multiplication - create empty result
        dense <- matrix(numeric(0), nrow = basis_dims[1], ncol = coeff_dims[2])
        
      } else {
        # Dimensions match - do the multiplication (which will result in empty matrix)
        dense <- basis %*% coeff
        
      }
    } else if (length(basis) == 0 && length(coeff) == 0) {
      # Both are empty - create empty result matrix
      dense <- matrix(numeric(0), nrow = 0, ncol = 0)
      
    } else {
      abort_lna("Invalid matrix dimensions for multiplication", .subclass="lna_error_internal", location="invert_step.temporal")
    }
  } else {
    # Special case: if basis is 0x0 (n_time=0, n_basis=0), create empty result with correct dimensions
    if (nrow(basis) == 0 && ncol(basis) == 0) {
      # When n_time=0, we need to reconstruct to match the expected output dimensions
      # The output should have 0 rows (time) and the same number of columns as the original data
      # We can infer the number of columns from the coefficient matrix
      n_features <- if (is.matrix(coeff) && ncol(coeff) > 0) ncol(coeff) else 1
      dense <- matrix(numeric(0), nrow = 0, ncol = n_features)
      
    } else {
      # Check that matrix dimensions are compatible for multiplication
      if (ncol(basis) != nrow(coeff)) {
        abort_lna(
          sprintf("Matrix dimension mismatch: basis has %d columns but coeff has %d rows", 
                  ncol(basis), nrow(coeff)),
          .subclass = "lna_error_internal", 
          location = "invert_step.temporal"
        )
      }

      dense <- temporal_reconstruct(desc$params$kind %||% "dct", basis, coeff)
    }
  }
  
  if (dbg) message(sprintf("[invert_step.temporal] Dense dims after matmult: %s", paste(dim(dense), collapse="x")))
  if (dbg && nrow(dense) >= 2 && ncol(dense) >= 2) {
    message("dense[1:2, 1:2]:")
    print(dense[1:2, 1:2, drop = FALSE])
  }
  
  subset <- handle$subset
  if (!is.null(subset$roi_mask)) {
    roi <- as.logical(subset$roi_mask)
    if (length(roi) == ncol(dense)) { 
      dense <- dense[, roi, drop = FALSE]
    }
  }
  if (!is.null(subset$time_idx)) {
    idx <- as.integer(subset$time_idx)
    # Ensure that idx is not empty and all indices are within bounds
    if (length(idx) > 0 && nrow(dense) >= max(idx) && min(idx) > 0) { 
        dense <- dense[idx, , drop = FALSE]
    } else if (length(idx) > 0) {
        # Handle out-of-bounds or empty idx if necessary, or let it error if that's desired.
        warning("time_idx for temporal subsetting is invalid or out of bounds.")
    }
  }
  if (dbg) message(sprintf("[invert_step.temporal] Dense dims after subsetting: %s", paste(dim(dense), collapse="x")))
  
  if (is.null(dense)) {
    abort_lna("Reconstructed data (dense) is NULL before stashing", .subclass="lna_error_internal", location="invert_step.temporal")
  }
  if (dbg) message(sprintf("[invert_step.temporal] Stashing to key: '%s'. Is dense NULL? %s", output_stash_key, is.null(dense)))
  new_values_list <- setNames(list(dense), output_stash_key)

  handle <- handle$update_stash(keys = character(),
                                new_values = new_values_list)
  if (dbg) message(sprintf("[invert_step.temporal] invert_step.temporal IS RETURNING handle with Stash keys: %s. Is input NULL? %s", paste(names(handle$stash), collapse=", "), is.null(handle$stash$input)))
  return(handle)
}

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
temporal_project.bspline <- function(kind, basis, X, ...) {
  qrB <- qr(basis)
  if (qrB$rank < ncol(basis)) {
    message("[temporal_project.bspline WARN] B-spline basis is rank deficient. Projection may be unstable.")
  }
  qr.coef(qrB, X)
}

#' Reconstruct data from temporal coefficients
#'
#' Mirrors `temporal_project` for the inverse operation. The default
#' simply multiplies the basis by the coefficient matrix.
#' @keywords internal
temporal_reconstruct <- function(kind, basis, coeff, ...) {
  stopifnot(is.character(kind), length(kind) == 1)
  obj <- structure(kind, class = c(kind, "character"))
  UseMethod("temporal_reconstruct", obj)
}

#' @export
temporal_reconstruct.default <- function(kind, basis, coeff, ...) {
  basis %*% coeff
}

#' @export
temporal_reconstruct.bspline <- function(kind, basis, coeff, ...) {
  basis %*% coeff
}

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
