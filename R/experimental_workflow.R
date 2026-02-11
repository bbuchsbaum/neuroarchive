#' Experimental Workflow and Utility Functions
#'
#' @description
#' High-level workflow functions and utilities for data preprocessing,
#' quantization, dimensionality reduction, and Rcpp acceleration control.
#'
#' @useDynLib neuroarchive, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @name experimental_workflow
NULL

#' Control Rcpp acceleration usage
#'
#' @description
#' Enable or disable Rcpp acceleration for neuroarchive functions.
#' When disabled, pure R implementations are used as fallbacks.
#'
#' @param enable Logical; if TRUE, enable Rcpp acceleration, if FALSE disable it.
#'   If NULL (default), just return current status without changing settings.
#' @param components Character vector specifying which components to control:
#'   "hrbf" (HRBF basis functions), "edge" (edge detection), "hwt" (Haar wavelets),
#'   or "all" (default)
#' @param verbose Logical; if TRUE, print current acceleration status
#'
#' @return Invisibly returns a list of current option values
#'
#' @examples
#' \donttest{
#' # Check current status
#' rcpp_control(verbose = TRUE)
#'
#' # Disable all Rcpp acceleration
#' rcpp_control(FALSE)
#'
#' # Enable only HRBF acceleration
#' rcpp_control(TRUE, components = "hrbf")
#'
#' # Re-enable all
#' rcpp_control(TRUE)
#' }
#'
#' @export
rcpp_control <- function(enable = NULL, components = "all", verbose = FALSE) {
  components <- match.arg(components, c("all", "hrbf", "edge", "hwt"), several.ok = TRUE)

  if ("all" %in% components) {
    components <- c("hrbf", "edge", "hwt")
  }

  option_map <- list(
    hrbf = "lna.hrbf.use_rcpp_helpers",
    edge = "lna.edge_adaptive.use_rcpp",
    hwt = "lna.hwt.use_rcpp"
  )

  # Set options only if enable is explicitly provided
  if (!is.null(enable)) {
    for (comp in components) {
      if (comp %in% names(option_map)) {
        options(setNames(list(enable), option_map[[comp]]))
      }
    }
  }

  # Get current status AFTER setting options (don't use defaults to see actual set values)
  current_status <- list(
    hrbf = getOption("lna.hrbf.use_rcpp_helpers"),
    edge = getOption("lna.edge_adaptive.use_rcpp"),
    hwt = getOption("lna.hwt.use_rcpp")
  )

  # Handle NULL values (not set) as TRUE (default behavior)
  current_status$hrbf <- if (is.null(current_status$hrbf)) TRUE else current_status$hrbf
  current_status$edge <- if (is.null(current_status$edge)) TRUE else current_status$edge
  current_status$hwt <- if (is.null(current_status$hwt)) TRUE else current_status$hwt

  if (verbose) {
    cat("Rcpp acceleration status:\n")
    cat(sprintf("  HRBF functions: %s\n",
                if(current_status$hrbf) "enabled" else "disabled"))
    cat(sprintf("  Edge detection: %s\n",
                if(current_status$edge) "enabled" else "disabled"))
    cat(sprintf("  Haar wavelets:  %s\n",
                if(current_status$hwt) "enabled" else "disabled"))
  }

  invisible(current_status)
}

#' Quantize data matrix with advanced compression options
#'
#' @description
#' Apply quantization to reduce data precision for compression.
#' Supports both global and per-voxel quantization strategies with
#' intelligent error monitoring and automatic parameter adjustment.
#'
#' @param data Matrix to quantize (time x voxels or any numeric matrix)
#' @param bits Number of quantization bits (1-16, default: 8)
#' @param method Scaling method:
#'   \itemize{
#'     \item \code{"range"}: Scale by min/max values
#'     \item \code{"sd"}: Scale by standard deviation (robust to outliers)
#'     \item \code{"robust"}: Scale by median and MAD
#'   }
#' @param scope Quantization scope:
#'   \itemize{
#'     \item \code{"global"}: Single scale factor for entire matrix
#'     \item \code{"voxel"}: Per-column scaling (better precision)
#'     \item \code{"adaptive"}: Choose scope based on data characteristics
#'   }
#' @param center Center data before quantization (default: TRUE)
#' @param clip_threshold Abort if clipping exceeds this percentage (default: 5)
#' @param auto_adjust Automatically increase bits if clipping is excessive (default: TRUE)
#' @param return_stats Return detailed compression statistics (default: FALSE)
#'
#' @return List with quantized data and reconstruction parameters:
#' \itemize{
#'   \item \code{quantized_data}: Quantized matrix (integer values)
#'   \item \code{params}: Parameters needed for reconstruction
#'   \item \code{compression_ratio}: Achieved compression ratio
#'   \item \code{clip_report}: Statistics on value clipping
#'   \item \code{mse}: Mean squared error (if \code{return_stats = TRUE})
#' }
#'
#' @examples
#' \donttest{
#' # Create example fMRI data matrix
#' set.seed(123)
#' fmri_data <- matrix(rnorm(200 * 1000, mean = 1000, sd = 50),
#'                     nrow = 200, ncol = 1000)
#'
#' # Basic 8-bit quantization with global scaling
#' quant_basic <- quantize_data(fmri_data, bits = 8, scope = "global")
#'
#' # High-precision per-voxel quantization
#' quant_precise <- quantize_data(
#'   fmri_data,
#'   bits = 12,
#'   scope = "voxel",
#'   method = "robust",
#'   return_stats = TRUE
#' )
#'
#' # Adaptive quantization (automatically chooses best parameters)
#' quant_adaptive <- quantize_data(
#'   fmri_data,
#'   scope = "adaptive",
#'   auto_adjust = TRUE
#' )
#'
#' # Reconstruct and check error
#' reconstructed <- dequantize_data(quant_precise)
#' mse <- mean((fmri_data - reconstructed)^2)
#' cat("Compression ratio:", quant_precise$compression_ratio, "\n")
#' cat("MSE:", mse, "\n")
#' }
#'
#' @export
quantize_data <- function(data,
                         bits = 8L,
                         method = c("range", "sd", "robust"),
                         scope = c("global", "voxel", "adaptive"),
                         center = TRUE,
                         clip_threshold = 5.0,
                         auto_adjust = TRUE,
                         return_stats = FALSE) {

  method <- match.arg(method)
  scope <- match.arg(scope)

  if (!is.array(data)) {
    if (is.vector(data) && is.atomic(data)) {
      dim(data) <- length(data)
    } else {
      abort_lna(
        "data must be a numeric array/matrix",
        .subclass = "lna_error_validation",
        location = "experimental_workflow::quantize_data:data"
      )
    }
  }
  if (!is.numeric(data)) {
    abort_lna(
      "data must be numeric",
      .subclass = "lna_error_validation",
      location = "experimental_workflow::quantize_data:data"
    )
  }
  if (!(is.numeric(bits) && length(bits) == 1 && !is.na(bits) &&
        bits == as.integer(bits) && bits >= 1 && bits <= 16)) {
    abort_lna(
      "bits must be an integer between 1 and 16",
      .subclass = "lna_error_validation",
      location = "experimental_workflow::quantize_data:bits"
    )
  }
  bits <- as.integer(bits)
  if (!(is.logical(center) && length(center) == 1 && !is.na(center))) {
    abort_lna(
      "center must be a single logical value",
      .subclass = "lna_error_validation",
      location = "experimental_workflow::quantize_data:center"
    )
  }
  if (!(is.numeric(clip_threshold) && length(clip_threshold) == 1 &&
        !is.na(clip_threshold) && clip_threshold >= 0)) {
    abort_lna(
      "clip_threshold must be a single non-negative number",
      .subclass = "lna_error_validation",
      location = "experimental_workflow::quantize_data:clip_threshold"
    )
  }
  if (!(is.logical(auto_adjust) && length(auto_adjust) == 1 && !is.na(auto_adjust))) {
    abort_lna(
      "auto_adjust must be a single logical value",
      .subclass = "lna_error_validation",
      location = "experimental_workflow::quantize_data:auto_adjust"
    )
  }
  if (!(is.logical(return_stats) && length(return_stats) == 1 && !is.na(return_stats))) {
    abort_lna(
      "return_stats must be a single logical value",
      .subclass = "lna_error_validation",
      location = "experimental_workflow::quantize_data:return_stats"
    )
  }
  if (any(!is.finite(data))) {
    abort_lna(
      "quantize_data cannot handle non-finite values",
      .subclass = "lna_error_validation",
      location = "experimental_workflow::quantize_data:finite"
    )
  }

  # Auto-select scope if adaptive
  if (scope == "adaptive") {
    d <- dim(data)
    if (length(d) == 2) {
      n_voxels <- ncol(data)
      col_ranges <- apply(data, 2, function(x) diff(range(x, na.rm = TRUE)))
      range_cv <- stats::sd(col_ranges, na.rm = TRUE) / mean(col_ranges, na.rm = TRUE)
      if (is.na(range_cv)) range_cv <- 0
      scope <- if (n_voxels < 100 || range_cv < 0.5) "global" else "voxel"
    } else if (length(d) == 4) {
      n_voxels <- prod(d[1:3])
      if (n_voxels < 100) {
        scope <- "global"
      } else {
        # Quick heuristic: estimate range variability across a small sample of voxels
        mat <- matrix(as.numeric(data), nrow = n_voxels, ncol = d[4])
        sample_n <- min(n_voxels, 512L)
        idx <- sample.int(n_voxels, sample_n)
        rng <- matrixStats::rowRanges(mat[idx, , drop = FALSE])
        col_ranges <- rng[, 2] - rng[, 1]
        range_cv <- stats::sd(col_ranges, na.rm = TRUE) / mean(col_ranges, na.rm = TRUE)
        if (is.na(range_cv)) range_cv <- 0
        scope <- if (range_cv < 0.5) "global" else "voxel"
      }
    } else {
      scope <- "global"
    }
  }

  robust_global <- function(x, bits, center) {
    med <- stats::median(x)
    s <- stats::mad(x, center = med, constant = 1.4826)
    if (center) {
      max_abs <- 3 * s
      scale <- (2 * max_abs) / (2^bits - 1)
      offset <- med - max_abs
    } else {
      lo <- med - 3 * s
      hi <- med + 3 * s
      scale <- (hi - lo) / (2^bits - 1)
      offset <- lo
    }
    if (scale == 0) {
      q_raw <- rep.int(0L, length(x))
      nclip <- 0L
      scale <- 1
    } else {
      q_raw <- round((x - offset) / scale)
      nclip <- sum(q_raw < 0 | q_raw > 2^bits - 1)
    }
    q <- q_raw
    q[q < 0] <- 0L
    q[q > 2^bits - 1] <- 2^bits - 1L
    clip_pct <- if (length(x) > 0) 100 * nclip / length(x) else 0
    list(q = q, scale = scale, offset = offset,
         n_clipped_total = as.integer(nclip), clip_pct = clip_pct)
  }

  robust_voxel_4d <- function(x, bits, center) {
    dims <- dim(x)
    vox <- prod(dims[1:3])
    time <- dims[4]
    mat <- matrix(as.numeric(x), nrow = vox, ncol = time)
    med <- matrixStats::rowMedians(mat)
    s <- matrixStats::rowMads(mat, center = med, constant = 1.4826)
    if (center) {
      max_abs <- 3 * s
      scale <- (2 * max_abs) / (2^bits - 1)
      offset <- med - max_abs
    } else {
      lo <- med - 3 * s
      hi <- med + 3 * s
      scale <- (hi - lo) / (2^bits - 1)
      offset <- lo
    }
    zero_idx <- scale == 0
    scale[zero_idx] <- 1
    q_raw <- round((mat - offset) / scale)
    nclip <- sum(q_raw < 0 | q_raw > 2^bits - 1)
    q <- q_raw
    q[q < 0] <- 0L
    q[q > 2^bits - 1] <- 2^bits - 1L
    clip_pct <- if (length(x) > 0) 100 * nclip / length(x) else 0
    list(
      q = array(as.integer(q), dim = dims),
      scale = array(scale, dim = dims[1:3]),
      offset = array(offset, dim = dims[1:3]),
      n_clipped_total = as.integer(nclip),
      clip_pct = clip_pct
    )
  }

  voxel_matrix <- function(x, bits, method, center) {
    X <- as.matrix(x)
    if (method == "robust") {
      med <- matrixStats::colMedians(X)
      s <- matrixStats::colMads(X, center = med, constant = 1.4826)
      if (center) {
        max_abs <- 3 * s
        scale <- (2 * max_abs) / (2^bits - 1)
        offset <- med - max_abs
      } else {
        lo <- med - 3 * s
        hi <- med + 3 * s
        scale <- (hi - lo) / (2^bits - 1)
        offset <- lo
      }
      zero_idx <- scale == 0
      scale[zero_idx] <- 1
      q_raw <- round((X - rep(offset, each = nrow(X))) / rep(scale, each = nrow(X)))
      nclip <- sum(q_raw < 0 | q_raw > 2^bits - 1)
      q <- q_raw
      q[q < 0] <- 0L
      q[q > 2^bits - 1] <- 2^bits - 1L
      clip_pct <- if (length(X) > 0) 100 * nclip / length(X) else 0
      list(q = matrix(as.integer(q), nrow = nrow(X), ncol = ncol(X)),
           scale = scale, offset = offset,
           n_clipped_total = as.integer(nclip), clip_pct = clip_pct)
    } else {
      # Fall back to per-column global quantization.
      q <- matrix(0L, nrow = nrow(X), ncol = ncol(X))
      scale <- numeric(ncol(X))
      offset <- numeric(ncol(X))
      nclip <- 0L
      for (j in seq_len(ncol(X))) {
        res <- .quantize_global(X[, j], bits, method, center)
        q[, j] <- as.integer(res$q)
        scale[j] <- as.numeric(res$scale)
        offset[j] <- as.numeric(res$offset)
        nclip <- nclip + (res$n_clipped_total %||% 0L)
      }
      clip_pct <- if (length(X) > 0) 100 * nclip / length(X) else 0
      list(q = q, scale = scale, offset = offset,
           n_clipped_total = as.integer(nclip), clip_pct = clip_pct)
    }
  }

  quantize_once <- function(bits_use) {
    d <- dim(data)
    if (identical(scope, "global")) {
      x <- as.numeric(data)
      res <- if (method == "robust") robust_global(x, bits_use, center) else .quantize_global(x, bits_use, method, center)
      res$q <- array(as.integer(res$q), dim = d)
      res$scale <- as.numeric(res$scale)
      res$offset <- as.numeric(res$offset)
      return(res)
    }

    # scope == "voxel"
    if (length(d) == 4) {
      if (method == "robust") {
        return(robust_voxel_4d(data, bits_use, center))
      }
      # Collect clipping stats for reporting
      return(.quantize_voxel_block(data, bits_use, method, center))
    }

    if (length(d) == 2) {
      return(voxel_matrix(data, bits_use, method, center))
    }

    abort_lna(
      "scope='voxel' currently supports 2D matrices or 4D arrays",
      .subclass = "lna_error_validation",
      location = "experimental_workflow::quantize_data:scope"
    )
  }

  res <- quantize_once(bits)
  while (isTRUE(auto_adjust) && isTRUE(res$clip_pct > clip_threshold) && bits < 16L) {
    bits <- min(16L, bits + 2L)
    res <- quantize_once(bits)
  }

  original_size <- as.numeric(object.size(data))
  quantized_size <- length(data) * as.numeric(bits) / 8
  compression_ratio <- if (quantized_size > 0) original_size / quantized_size else Inf

  out <- list(
    q = res$q,
    quantized_data = res$q,
    scale = res$scale,
    offset = res$offset,
    bits = bits,
    method = method,
    scope = scope,
    center = center,
    params = list(
      bits = bits,
      method = method,
      scale_scope = scope,
      center = center,
      scale = res$scale,
      offset = res$offset
    ),
    compression_ratio = compression_ratio,
    clip_report = list(
      n_clipped_total = res$n_clipped_total %||% 0L,
      clip_pct = res$clip_pct %||% 0
    )
  )

  if (isTRUE(return_stats)) {
    xhat <- dequantize_data(out)
    mse <- mean((as.numeric(data) - as.numeric(xhat))^2)
    var_x <- stats::var(as.numeric(data))
    snr_db <- if (is.finite(mse) && mse > 0 && is.finite(var_x) && var_x > 0) {
      10 * log10(var_x / mse)
    } else if (isTRUE(mse == 0)) {
      Inf
    } else {
      NA_real_
    }
    out$mse <- mse
    out$snr_db <- snr_db
  }

  class(out) <- c("lna_quant_result", "list")
  out
}

#' Dequantize Data
#'
#' @description
#' Reconstruct original data from a quantized representation produced by
#' \code{\link{quantize_data}}. Reverses the quantization by applying stored
#' scale/offset parameters to map integer values back to floating-point.
#'
#' @param quant_result Result from quantize_data()
#' @return Numeric matrix of reconstructed (dequantized) data
#' @export
dequantize_data <- function(quant_result) {
  if (!is.list(quant_result)) {
    abort_lna(
      "quant_result must be a list (as returned by quantize_data())",
      .subclass = "lna_error_validation",
      location = "experimental_workflow::dequantize_data:quant_result"
    )
  }
  q <- quant_result$q %||% quant_result$quantized_data
  if (is.null(q)) {
    abort_lna(
      "quant_result is missing quantized data (field 'q' or 'quantized_data')",
      .subclass = "lna_error_validation",
      location = "experimental_workflow::dequantize_data:q"
    )
  }
  scale <- quant_result$scale %||% quant_result$params$scale
  offset <- quant_result$offset %||% quant_result$params$offset
  if (is.null(scale) || is.null(offset)) {
    abort_lna(
      "quant_result is missing reconstruction parameters (scale/offset)",
      .subclass = "lna_error_validation",
      location = "experimental_workflow::dequantize_data:params"
    )
  }

  # Common cases:
  # - global: scale/offset are scalars
  # - voxel, 2D: q is matrix, scale/offset are vectors (per-column)
  # - voxel, 4D: q is 4D, scale/offset are 3D (per-voxel, broadcast over time)
  dq <- dim(q)
  ds <- dim(scale)
  do <- dim(offset)

  if (is.null(ds) && is.null(do)) {
    return(q * scale + offset)
  }

  if (is.matrix(q) && is.null(ds) && is.null(do) &&
      length(scale) == ncol(q) && length(offset) == ncol(q)) {
    out <- sweep(q, 2, scale, FUN = "*")
    out <- sweep(out, 2, offset, FUN = "+")
    return(out)
  }

  if (!is.null(dq) && !is.null(ds) && !is.null(do) &&
      length(dq) == length(ds) + 1L &&
      identical(ds, do) &&
      all(ds == dq[seq_along(ds)])) {
    margins <- seq_along(ds)
    out <- sweep(q, margins, scale, FUN = "*")
    out <- sweep(out, margins, offset, FUN = "+")
    return(out)
  }

  # Fallback: let base R try; will error if non-conformable.
  q * scale + offset
}

#' Sparse Principal Component Analysis
#'
#' @description
#' Perform sparse PCA to find interpretable components with automatic
#' or manual sparsity control. Useful for dimensionality reduction
#' while maintaining interpretability.
#'
#' @param data Matrix with observations in rows, variables in columns
#' @param n_components Number of components to extract (default: 10)
#' @param sparsity_type Method for inducing sparsity: "lasso", "elastic_net", "hard_threshold"
#' @param sparsity_param Sparsity parameter (lambda for lasso, proportion for threshold)
#' @param max_iter Maximum iterations for convergence (default: 100)
#' @param tolerance Convergence tolerance (default: 1e-6)
#' @param standardize Standardize variables before analysis (default: TRUE)
#'
#' @return List with components, loadings, and diagnostic information
#'
#' @examples
#' \donttest{
#' # Create example data matrix (observations x variables)
#' set.seed(123)
#' my_data <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
#' # Add some structure
#' my_data[, 1:10] <- my_data[, 1:10] + rnorm(100, mean = 2)
#'
#' # Basic sparse PCA
#' spca_result <- sparse_pca(my_data, n_components = 5)
#'
#' # High sparsity for interpretability
#' sparse_result <- sparse_pca(
#'   my_data,
#'   sparsity_type = "hard_threshold",
#'   sparsity_param = 0.1  # Keep only 10% of loadings
#' )
#' }
#'
#' @export
sparse_pca <- function(data,
                      n_components = 10L,
                      sparsity_type = c("lasso", "elastic_net", "hard_threshold"),
                      sparsity_param = 0.1,
                      max_iter = 100L,
                      tolerance = 1e-6,
                      standardize = TRUE) {

  sparsity_type <- match.arg(sparsity_type)

  if (!(is.matrix(data) && is.numeric(data))) {
    abort_lna(
      "data must be a numeric matrix (observations x variables)",
      .subclass = "lna_error_validation",
      location = "experimental_workflow::sparse_pca:data"
    )
  }
  if (!(is.numeric(n_components) && length(n_components) == 1 && !is.na(n_components) &&
        n_components == as.integer(n_components) && n_components >= 0)) {
    abort_lna(
      "n_components must be a non-negative integer",
      .subclass = "lna_error_validation",
      location = "experimental_workflow::sparse_pca:n_components"
    )
  }
  if (!(is.logical(standardize) && length(standardize) == 1 && !is.na(standardize))) {
    abort_lna(
      "standardize must be a single logical value",
      .subclass = "lna_error_validation",
      location = "experimental_workflow::sparse_pca:standardize"
    )
  }
  if (!(is.numeric(sparsity_param) && length(sparsity_param) == 1 && !is.na(sparsity_param))) {
    abort_lna(
      "sparsity_param must be a single numeric value",
      .subclass = "lna_error_validation",
      location = "experimental_workflow::sparse_pca:sparsity_param"
    )
  }

  X <- data
  center_vec <- NULL
  scale_vec <- NULL
  if (isTRUE(standardize)) {
    center_vec <- colMeans(X)
    scale_vec <- apply(X, 2, stats::sd)
    scale_vec[scale_vec == 0] <- 1
    X <- scale(X, center = center_vec, scale = scale_vec)
  }

  k <- as.integer(n_components)
  if (k == 0L) {
    out <- list(
      loadings = matrix(numeric(0), nrow = 0, ncol = ncol(X)),
      scores = matrix(numeric(0), nrow = nrow(X), ncol = 0),
      backend = "k_zero_case",
      params = list(
        k = k,
        sparsity_type = sparsity_type,
        sparsity_param = sparsity_param,
        max_iter = as.integer(max_iter),
        tolerance = tolerance,
        standardize = standardize
      ),
      center = center_vec,
      scale = scale_vec
    )
    class(out) <- c("lna_sparsepca_result", "list")
    return(out)
  }

  k_eff <- min(k, min(dim(X)))
  backend <- "svd"
  loadings <- NULL
  scores <- NULL
  d <- NULL

  if (sparsity_type %in% c("lasso", "elastic_net") &&
      requireNamespace("sparsepca", quietly = TRUE)) {
    backend <- "sparsepca"
    fit <- sparsepca::spca(
      X,
      k = k_eff,
      alpha = sparsity_param,
      verbose = FALSE,
      center = FALSE,
      scale = FALSE
    )
    # sparsepca returns p x k loadings
    loadings <- t(fit$loadings)
    scores <- X %*% fit$loadings
    d <- fit$d %||% rep(NA_real_, ncol(scores))
  } else {
    sv <- svd(X, nu = k_eff, nv = k_eff)
    d <- sv$d[seq_len(k_eff)]
    loadings <- t(sv$v[, seq_len(k_eff), drop = FALSE])
    scores <- sv$u[, seq_len(k_eff), drop = FALSE] %*% diag(d, nrow = k_eff, ncol = k_eff)

    # Optional hard-threshold sparsification.
    if (sparsity_type == "hard_threshold") {
      if (!(sparsity_param > 0 && sparsity_param <= 1)) {
        abort_lna(
          "For sparsity_type='hard_threshold', sparsity_param must be in (0, 1]",
          .subclass = "lna_error_validation",
          location = "experimental_workflow::sparse_pca:hard_threshold"
        )
      }
      keep_n <- max(1L, ceiling(ncol(loadings) * sparsity_param))
      for (i in seq_len(nrow(loadings))) {
        ord <- order(abs(loadings[i, ]), decreasing = TRUE)
        drop_idx <- ord[(keep_n + 1L):length(ord)]
        if (length(drop_idx) > 0) loadings[i, drop_idx] <- 0
      }
    }
  }

  out <- list(
    loadings = loadings,
    scores = scores,
    singular_values = d,
    backend = backend,
    params = list(
      k = k,
      k_eff = k_eff,
      sparsity_type = sparsity_type,
      sparsity_param = sparsity_param,
      max_iter = as.integer(max_iter),
      tolerance = tolerance,
      standardize = standardize
    ),
    center = center_vec,
    scale = scale_vec
  )
  class(out) <- c("lna_sparsepca_result", "list")
  out
}

#' Complete fMRI preprocessing and compression workflow
#'
#' @description
#' High-level function that combines multiple neuroarchive transforms
#' for a complete fMRI preprocessing and compression workflow. This
#' function automatically suggests optimal parameters based on data
#' characteristics and study type.
#'
#' @param data Matrix with time in rows, voxels in columns
#' @param mask LogicalNeuroVol defining spatial structure (for spatial transforms)
#' @param TR Repetition time in seconds (for temporal transforms)
#' @param study_type Type of fMRI study: "resting", "task", "event"
#' @param workflow Type of workflow:
#'   \itemize{
#'     \item \code{"compress"}: Temporal + spatial compression
#'     \item \code{"denoise"}: Temporal denoising + compression
#'     \item \code{"full"}: Complete preprocessing pipeline
#'   }
#' @param compression_target Target compression ratio (default: 10)
#' @param denoise_level Denoising aggressiveness: "conservative", "moderate", "aggressive"
#' @param return_intermediate Return intermediate results for inspection (default: FALSE)
#' @param verbose Print progress and recommendations (default: TRUE)
#'
#' @return List containing:
#' \itemize{
#'   \item \code{compressed_data}: Final compressed representation
#'   \item \code{compression_ratio}: Achieved compression ratio
#'   \item \code{parameters}: Parameters used for each transform
#'   \item \code{reconstruction_function}: Function to reconstruct original data
#'   \item \code{intermediate_results}: Intermediate outputs (if requested)
#'   \item \code{recommendations}: Suggestions for parameter optimization
#' }
#'
#' @examples
#' \donttest{
#' # Create example fMRI data and mask
#' library(neuroim2)
#' set.seed(123)
#'
#' # Simulate 64x64x20 brain with 300 timepoints, TR=2s
#' dims <- c(64, 64, 20)
#' mask_array <- array(TRUE, dims)
#' # Create a brain-like mask (remove edges)
#' mask_array[1:5, , ] <- FALSE
#' mask_array[60:64, , ] <- FALSE
#' mask_array[, 1:5, ] <- FALSE
#' mask_array[, 60:64, ] <- FALSE
#'
#' mask <- LogicalNeuroVol(mask_array, NeuroSpace(dims, spacing = c(3, 3, 3)))
#' n_voxels <- sum(mask_array)
#'
#' # Simulate fMRI data with realistic characteristics
#' fmri_data <- matrix(rnorm(300 * n_voxels, mean = 1000, sd = 50),
#'                     nrow = 300, ncol = n_voxels)
#' # Add temporal drift
#' drift <- outer(1:300, rep(1, n_voxels)) * 0.5
#' fmri_data <- fmri_data + drift
#'
#' # Conservative compression workflow for resting-state
#' result_rest <- fmri_workflow(
#'   fmri_data, mask,
#'   TR = 2.0,
#'   study_type = "resting",
#'   workflow = "compress",
#'   compression_target = 20,
#'   denoise_level = "conservative"
#' )
#'
#' # Aggressive denoising workflow for task fMRI
#' result_task <- fmri_workflow(
#'   fmri_data, mask,
#'   TR = 2.0,
#'   study_type = "task",
#'   workflow = "denoise",
#'   denoise_level = "moderate",
#'   return_intermediate = TRUE
#' )
#'
#' # Full preprocessing pipeline
#' result_full <- fmri_workflow(
#'   fmri_data, mask,
#'   TR = 2.0,
#'   study_type = "event",
#'   workflow = "full",
#'   compression_target = 15
#' )
#'
#' # Reconstruct data to check quality
#' reconstructed <- result_rest$reconstruction_function()
#' mse <- mean((fmri_data - reconstructed)^2)
#' cat("Compression ratio:", result_rest$compression_ratio, "\n")
#' cat("Reconstruction MSE:", mse, "\n")
#' }
#'
#' @export
fmri_workflow <- function(data,
                         mask = NULL,
                         TR = 2.0,
                         study_type = c("resting", "task", "event"),
                         workflow = c("compress", "denoise", "full"),
                         compression_target = 10,
                         denoise_level = c("conservative", "moderate", "aggressive"),
                         return_intermediate = FALSE,
                         verbose = TRUE) {

  study_type <- match.arg(study_type)
  workflow <- match.arg(workflow)
  denoise_level <- match.arg(denoise_level)

  stopifnot(is.matrix(data), is.numeric(data))
  n_time <- nrow(data)
  n_voxels <- ncol(data)

  if (verbose) {
    cat("=== fMRI Workflow Configuration ===\n")
    cat(sprintf("Data: %d timepoints × %d voxels\n", n_time, n_voxels))
    cat(sprintf("TR: %.2f seconds\n", TR))
    cat(sprintf("Study type: %s\n", study_type))
    cat(sprintf("Workflow: %s\n", workflow))
    cat(sprintf("Target compression: %dx\n", compression_target))
  }

  # Container for all parameters and intermediate results
  workflow_params <- list()
  intermediate_results <- list()
  original_data <- data

  # Step 1: Temporal preprocessing and compression
  if (verbose) cat("\n--- Step 1: Temporal Transform ---\n")

  # Get optimal DPSS parameters for this study type
  dpss_params <- suggest_dpss_fmri(TR, n_time, study_type)
  if (verbose) {
    cat(sprintf("DPSS parameters: NW=%.2f, n_basis=%d\n",
                dpss_params$time_bandwidth_product, dpss_params$n_basis))
    cat(sprintf("Preserved frequency range: 0-%.3f Hz\n", dpss_params$preserved_freq))
  }

  # Configure thresholding based on workflow and denoise level
  if (workflow %in% c("denoise", "full")) {
    threshold_config <- switch(denoise_level,
      "conservative" = list(type = "energy", keep_energy = 0.98),
      "moderate" = list(type = "energy", keep_energy = 0.95),
      "aggressive" = list(type = "energy", keep_energy = 0.90)
    )
  } else {
    threshold_config <- list(type = "none")
  }

  # Apply temporal transform
  temporal_result <- temporal_project_denoise(
    data,
    basis_type = "dpss",
    n_components = dpss_params$n_basis,
    time_bandwidth_product = dpss_params$time_bandwidth_product,
    threshold_type = threshold_config$type,
    keep_energy = threshold_config$keep_energy,
    return_basis = TRUE
  )

  workflow_params$temporal <- list(
    basis_type = "dpss",
    time_bandwidth_product = dpss_params$time_bandwidth_product,
    n_basis = dpss_params$n_basis,
    threshold_type = threshold_config$type,
    keep_energy = threshold_config$keep_energy
  )

  if (return_intermediate) {
    intermediate_results$temporal <- temporal_result
  }

  current_data <- temporal_result$coefficients
  temporal_compression <- temporal_result$compression_ratio

  if (verbose) {
    cat(sprintf("Temporal compression: %.1fx\n", temporal_compression))
    if (!is.null(temporal_result$denoising_stats)) {
      cat(sprintf("Coefficient sparsity: %.1f%%\n",
                  temporal_result$denoising_stats$sparsity * 100))
    }
  }

  # Step 2: Spatial compression (if mask provided and workflow requires it)
  spatial_compression <- 1.0
  spatial_result <- NULL

  if (!is.null(mask) && workflow %in% c("compress", "full")) {
    if (verbose) cat("\n--- Step 2: Spatial Transform ---\n")

    # Determine spatial compression target based on overall compression goal
    remaining_compression_needed <- compression_target / temporal_compression

    if (remaining_compression_needed > 1.5) {
      spatial_result <- basis_space_embed(
        t(current_data),  # Transpose for spatial embedding
        mask,
        compression_ratio = 1 / remaining_compression_needed,
        return_basis = return_intermediate
      )

      current_data <- t(spatial_result$coefficients)  # Transpose back
      spatial_compression <- spatial_result$compression_ratio

      workflow_params$spatial <- spatial_result$basis_params

      if (return_intermediate) {
        intermediate_results$spatial <- spatial_result
      }

      if (verbose) {
        cat(sprintf("Spatial compression: %.1fx\n", spatial_compression))
      }
    } else {
      if (verbose) cat("Spatial compression not needed (temporal compression sufficient)\n")
    }
  }

  # Step 3: Quantization (for full workflow)
  quantization_compression <- 1.0
  quantization_result <- NULL

  if (workflow == "full") {
    if (verbose) cat("\n--- Step 3: Quantization ---\n")

    # Use adaptive quantization
    quantization_result <- quantize_data(
      current_data,
      bits = 8,
      scope = "adaptive",
      method = "robust",
      auto_adjust = TRUE,
      return_stats = TRUE
    )

    current_data <- quantization_result$quantized_data
    quantization_compression <- quantization_result$compression_ratio

    workflow_params$quantization <- quantization_result$params

    if (return_intermediate) {
      intermediate_results$quantization <- quantization_result
    }

    if (verbose) {
      cat(sprintf("Quantization compression: %.1fx\n", quantization_compression))
      if (!is.null(quantization_result$mse)) {
        cat(sprintf("Quantization SNR: %.1f dB\n", quantization_result$snr_db))
      }
    }
  }

  # Calculate total compression
  total_compression <- temporal_compression * spatial_compression * quantization_compression

  if (verbose) {
    cat("\n=== Final Results ===\n")
    cat(sprintf("Total compression ratio: %.1fx\n", total_compression))
    cat(sprintf("Target compression: %.1fx\n", compression_target))

    if (total_compression >= compression_target * 0.8) {
      cat("✓ Compression target achieved\n")
    } else {
      cat("⚠ Compression target not fully achieved\n")
    }
  }

  # Create reconstruction function
  reconstruction_function <- function(output_file = NULL) {
    # Reverse the transform pipeline to reconstruct original data
    reconstructed <- current_data

    # Step 1: Dequantize if quantization was applied
    if (!is.null(quantization_result)) {
      reconstructed <- dequantize_data(quantization_result)
    }

    # Step 2: Reconstruct from spatial basis if spatial transform was applied
    if (!is.null(spatial_result)) {
      # Transpose to get proper orientation for reconstruction
      spatial_coeffs <- list(
        coefficients = t(reconstructed),
        params = spatial_result$params,
        basis = if (!is.null(spatial_result$basis)) spatial_result$basis else NULL
      )
      class(spatial_coeffs) <- "lna_basis_embed_result"
      reconstructed <- t(basis_space_reconstruct(spatial_coeffs, mask))
    }

    # Step 3: Reconstruct from temporal basis
    if (!is.null(temporal_result$basis)) {
      reconstructed <- temporal_result$basis %*% reconstructed
    }

    return(reconstructed)
  }

  # Generate recommendations
  recommendations <- list()

  if (total_compression < compression_target * 0.8) {
    recommendations <- append(recommendations,
      "Consider more aggressive temporal thresholding or higher spatial compression")
  }

  if (workflow != "full" && compression_target > 15) {
    recommendations <- append(recommendations,
      "Consider 'full' workflow with quantization for higher compression ratios")
  }

  if (study_type == "resting" && denoise_level == "conservative") {
    recommendations <- append(recommendations,
      "Resting-state data may benefit from moderate denoising")
  }

  # Return comprehensive results
  list(
    compressed_data = current_data,
    compression_ratio = total_compression,
    parameters = workflow_params,
    reconstruction_function = reconstruction_function,
    intermediate_results = if (return_intermediate) intermediate_results else NULL,
    recommendations = recommendations,
    workflow_summary = list(
      study_type = study_type,
      workflow = workflow,
      denoise_level = denoise_level,
      temporal_compression = temporal_compression,
      spatial_compression = spatial_compression,
      quantization_compression = quantization_compression
    )
  )
}
