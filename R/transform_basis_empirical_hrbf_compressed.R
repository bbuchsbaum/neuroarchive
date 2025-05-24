#' Inverse step for 'basis.empirical_hrbf_compressed'
#'
#' Reconstructs an empirical basis matrix from quantized HRBF
#' dictionary codes and the stored SVD \eqn{V^T} matrix. The HRBF
#' dictionary is regenerated from the descriptor referenced by
#' `hrbf_dictionary_descriptor_path`.
#' @keywords internal
invert_step.basis.empirical_hrbf_compressed <- function(type, desc, handle) {
  p <- desc$params %||% list()
  dict_path <- p$hrbf_dictionary_descriptor_path
  if (is.null(dict_path)) {
    abort_lna("hrbf_dictionary_descriptor_path missing in descriptor",
              .subclass = "lna_error_descriptor",
              location = "invert_step.basis.empirical_hrbf_compressed")
  }

  vt_path <- NULL
  if (!is.null(desc$datasets)) {
    roles <- vapply(desc$datasets, function(d) d$role, character(1))
    idx <- which(roles == "svd_vt")
    if (length(idx) > 0) vt_path <- desc$datasets[[idx[1]]]$path
  }
  if (is.null(vt_path)) {
    abort_lna("svd_vt path not found in descriptor datasets",
              .subclass = "lna_error_descriptor",
              location = "invert_step.basis.empirical_hrbf_compressed:vt_path")
  }

  codes_key <- desc$outputs[[1]] %||% "hrbf_codes"
  input_key <- desc$inputs[[1]] %||% "basis_matrix"
  if (!handle$has_key(codes_key)) {
    return(handle)
  }
  codes <- handle$get_inputs(codes_key)[[codes_key]]

  # Handle structured codes format from forward step
  if (is.list(codes) && length(codes) > 0 && is.list(codes[[1]])) {
    # Reconstruct sparse codes matrix directly - no dense intermediate
    n_components <- length(codes)
    
    # Load HRBF dictionary to get dimensions
    root <- handle$h5[["/"]]
    path_parts <- strsplit(dict_path, "/")[[1]]
    dname <- tail(path_parts, 1)
    gpath <- paste(head(path_parts, -1), collapse = "/")
    if (gpath == "") gpath <- "/"
    tf_group <- root[[gpath]]
    dict_desc <- read_json_descriptor(tf_group, dname)
    if (inherits(tf_group, "H5Group")) tf_group$close()
    
    mask_neurovol <- handle$mask_info$mask
    if (is.null(mask_neurovol)) {
      abort_lna("mask_info$mask missing",
                .subclass = "lna_error_validation",
                location = "invert_step.basis.empirical_hrbf_compressed:mask")
    }
    B_dict <- hrbf_basis_from_params(dict_desc$params, mask_neurovol)
    n_atoms <- nrow(B_dict)
    
    # Pre-allocate vectors to avoid quadratic concatenation
    lens <- vapply(codes, function(z) length(z$indices), integer(1))
    nz_tot <- sum(lens)
    i_idx <- integer(nz_tot)
    j_idx <- integer(nz_tot) 
    x_val <- numeric(nz_tot)
    pos <- 1L
    
    for (j in seq_len(n_components)) {
      m <- lens[j]
      if (m > 0) {
        # Dequantize weights
        weights <- codes[[j]]$q * codes[[j]]$scale
        rng <- pos:(pos + m - 1L)
        i_idx[rng] <- j
        j_idx[rng] <- codes[[j]]$indices
        x_val[rng] <- weights
        pos <- pos + m
      }
    }
    
    codes <- Matrix::sparseMatrix(i = i_idx, j = j_idx, x = x_val,
                                  dims = c(n_components, n_atoms))
  }

  root <- handle$h5[["/"]]
  Vt <- h5_read(root, vt_path)

  # B_dict should already be loaded from the structured-code branch above
  # If not, we have a logic error since codes should always be structured format
  if (!exists("B_dict")) {
    abort_lna("B_dict not found - codes format not recognized",
              .subclass = "lna_error_validation", 
              location = "invert_step.basis.empirical_hrbf_compressed:missing_dict")
  }

  bits <- p$omp_quant_bits %||% 5
  if (inherits(codes, "integer")) {
    codes_num <- as.numeric(codes)
  } else {
    codes_num <- codes
  }
  codes_num <- codes_num / (2^bits - 1)

  U_sigma <- codes_num %*% B_dict
  basis_reco <- t(Vt) %*% U_sigma
  
  # Convert to regular matrix if sparse
  if (inherits(basis_reco, "sparseMatrix")) {
    basis_reco <- as.matrix(basis_reco)
  }

  handle$update_stash(keys = codes_key,
                      new_values = setNames(list(basis_reco), input_key))
}

#' Default parameters for 'basis.empirical_hrbf_compressed'
#' @export
#' @keywords internal
lna_default.basis.empirical_hrbf_compressed <- function() {
  default_params("basis.empirical_hrbf_compressed")
}


#' Empirical HRBF Compressed Basis - Forward Step
#'
#' Compresses a dense empirical basis matrix via SVD and sparse
#' HRBF re-expansion. This is a minimal implementation following the
#' HRBF proposal. The dictionary regeneration and OMP coding are
#' simplified and may be extended.
#' @keywords internal
forward_step.basis.empirical_hrbf_compressed <- function(type, desc, handle) {
  p <- desc$params %||% list()
  svd_rank <- p$svd_rank %||% 120L
  omp_tol <- p$omp_tol %||% 0.01
  omp_sparsity_limit <- p$omp_sparsity_limit %||% 32L
  omp_quant_bits <- p$omp_quant_bits %||% 5L
  dict_path <- p$hrbf_dictionary_descriptor_path
  if (is.null(dict_path)) {
    abort_lna("hrbf_dictionary_descriptor_path required",
              .subclass = "lna_error_validation",
              location = "forward_step.basis.empirical_hrbf_compressed:param")
  }

  inp <- handle$pull_first(c("dense_basis_matrix", "basis_matrix", "input"))
  input_key <- inp$key
  B <- as_dense_mat(inp$value)

  r <- min(as.integer(svd_rank), min(dim(B)))
  # Use irlba only if the rank is small enough relative to matrix dimensions
  use_irlba <- requireNamespace("irlba", quietly = TRUE) && r < min(dim(B))
  if (use_irlba) {
    sv <- irlba::irlba(B, nv = r, nu = r)
    U <- sv$u
    V <- sv$v
    d <- sv$d
  } else {
    sv <- svd(B, nu = r, nv = r)
    U <- sv$u
    V <- sv$v
    d <- sv$d
  }
  U_sigma <- U[, seq_len(r), drop = FALSE] %*%
    diag(d[seq_len(r)], nrow = r)
  Vt <- t(V[, seq_len(r), drop = FALSE])

  root <- handle$h5[["/"]]
  dict_desc <- read_json_descriptor(root, dict_path)
  B_dict <- .regenerate_hrbf_basis(dict_desc$params, handle$mask_info$mask)
  D <- Matrix::t(B_dict) # voxels x atoms

  codes <- vector("list", ncol(U_sigma))
  for (j in seq_len(ncol(U_sigma))) {
    y <- U_sigma[, j]
    enc <- .omp_encode(y, D, tol = omp_tol,
                       max_nonzero = omp_sparsity_limit)
    q <- .quantize_weights(enc$coeff, bits = omp_quant_bits)
    codes[[j]] <- list(indices = as.integer(enc$idx),
                       q = q$q,
                       scale = q$scale,
                       max_iter_reached = isTRUE(enc$max_iter_reached))
  }

  plan <- handle$plan
  fname <- plan$get_next_filename(type)
  base <- tools::file_path_sans_ext(fname)
  vt_path <- paste0("/basis/", base, "/vt_matrix")
  codes_path <- paste0("/basis/", base, "/hrbf_codes")
  params_json <- as.character(jsonlite::toJSON(p, auto_unbox = TRUE))

  desc$version <- "1.0"
  desc$inputs <- c(input_key)
  desc$outputs <- c("hrbf_codes")
  desc$datasets <- list(list(path = vt_path, role = "svd_vt"),
                        list(path = codes_path, role = "hrbf_codes"))
  desc$params <- p

  plan$add_descriptor(fname, desc)
  plan$add_payload(vt_path, Vt)
  plan$add_dataset_def(vt_path, "svd_vt", as.character(type),
                       plan$origin_label, as.integer(plan$next_index - 1L),
                       params_json, vt_path, "eager", dtype = NA_character_)
  plan$add_payload(codes_path, codes)
  plan$add_dataset_def(codes_path, "hrbf_codes", as.character(type),
                       plan$origin_label, as.integer(plan$next_index - 1L),
                       params_json, codes_path, "eager", dtype = NA_character_)

  handle$plan <- plan
  handle$update_stash(keys = character(),
                      new_values = list(hrbf_vt = Vt,
                                         hrbf_codes = codes))
}

# ----------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------

.regenerate_hrbf_basis <- function(p, mask_neurovol) {
  sigma0 <- p$sigma0 %||% 6
  levels <- p$levels %||% 3L
  radius_factor <- p$radius_factor %||% 2.5

  num_extra_fine_levels <- p$num_extra_fine_levels %||% 0L
  kernel_type <- p$kernel_type %||% "gaussian"

  seed <- p$seed %||% 1L

  voxel_to_world <- function(vox_mat) {
    spc <- tryCatch(space(mask_neurovol), error = function(e) NULL)
    spacing_vec <- tryCatch(spacing(spc), error = function(e) c(1,1,1))
    origin_vec <- tryCatch(origin(spc), error = function(e) c(0,0,0))
    sweep(vox_mat - 1, 2, spacing_vec, `*`) +
      matrix(origin_vec, nrow(vox_mat), 3, byrow = TRUE)
  }

  centres_list <- list(); sigs <- numeric(); level_vec <- integer()
  for (j in seq_len(levels + 1L) - 1L) {
    sigma_j <- sigma0 / (2^j)
    r_j <- radius_factor * sigma_j
    vox_centres <- poisson_disk_sample_neuroim2(mask_neurovol, r_j, seed + j)
    if (nrow(vox_centres) > 0) {
      centres_list[[length(centres_list) + 1L]] <- voxel_to_world(vox_centres)
      n_new <- nrow(vox_centres)
      sigs <- c(sigs, rep(sigma_j, n_new))
      level_vec <- c(level_vec, rep(j, n_new))
    }
  }
  if (num_extra_fine_levels > 0L) {
    for (j_extra in seq_len(num_extra_fine_levels)) {
      j_level <- levels + j_extra
      sigma_new <- sigma0 / (2^(levels + j_extra))
      r_new <- radius_factor * sigma_new
      vox_centres <- poisson_disk_sample_neuroim2(
        mask_neurovol, r_new, seed + levels + j_extra
      )
      if (nrow(vox_centres) > 0) {
        centres_list[[length(centres_list) + 1L]] <- voxel_to_world(vox_centres)
        n_new <- nrow(vox_centres)
        sigs <- c(sigs, rep(sigma_new, n_new))
        level_vec <- c(level_vec, rep(j_level, n_new))
      }
    }
  }
  C_total <- if (length(centres_list) > 0) do.call(rbind, centres_list)
             else matrix(numeric(0), ncol = 3)
  sigma_vec <- sigs
  level_vec <- level_vec

  mask_arr <- as.array(mask_neurovol)
  mask_linear_indices <- which(mask_arr)
  vox_coords <- which(mask_arr, arr.ind = TRUE)
  mask_coords_world <- voxel_to_world(vox_coords)
  n_total_vox <- length(mask_arr)
  k_actual <- nrow(C_total)

  if (k_actual > 0) {
    # Pre-allocate lists to avoid quadratic concatenation
    trip_i <- vector("list", k_actual)
    trip_j <- vector("list", k_actual)
    trip_x <- vector("list", k_actual)
    
    for (kk in seq_len(k_actual)) {
      atom <- generate_hrbf_atom(mask_coords_world, mask_linear_indices,
                                 C_total[kk,], sigma_vec[kk],
                                 level_vec[kk], levels, p)
      trip_i[[kk]] <- rep.int(kk, length(atom$indices))
      trip_j[[kk]] <- atom$indices
      trip_x[[kk]] <- atom$values
    }
    
    # Flatten lists efficiently
    i_idx <- unlist(trip_i, use.names = FALSE)
    j_idx <- unlist(trip_j, use.names = FALSE)
    x_val <- unlist(trip_x, use.names = FALSE)
    
    Matrix::sparseMatrix(i = i_idx, j = j_idx, x = x_val,
                         dims = c(k_actual, n_total_vox))
  } else {
    Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(),
                         dims = c(0, n_total_vox))
  }
}

.omp_encode_R <- function(y, D, tol = 1e-3, max_nonzero = 32L) {
  r <- y
  idx <- integer(); coeff <- numeric()
  while (sum(r^2) > tol^2 && length(idx) < max_nonzero) {
    corr <- as.numeric(Matrix::crossprod(D, r))
    j <- which.max(abs(corr))
    idx <- unique(c(idx, j))
    D_sub <- D[, idx, drop = FALSE]
    coeff <- as.numeric(Matrix::solve(Matrix::crossprod(D_sub), Matrix::crossprod(D_sub, y)))
    r <- y - D_sub %*% coeff
  }
  list(idx = idx, coeff = coeff, max_iter_reached = length(idx) >= max_nonzero && sum(r^2) > tol^2)
}

.omp_encode <- function(y, D, tol = 1e-3, max_nonzero = 32L) {
  use_rcpp <- isTRUE(getOption("lna.hrbf.use_rcpp_helpers", TRUE)) &&
    exists("omp_encode_rcpp")

  if (use_rcpp) {
    res <- tryCatch(
      omp_encode_rcpp(as.numeric(y), D, tol^2, as.integer(max_nonzero)),
      error = function(e) NULL
    )
    if (!is.null(res)) {
      if (isTRUE(res$max_iter_reached)) {
        warning("OMP residual tolerance not reached before hitting sparsity limit", call. = FALSE)
      }
      return(list(idx = as.integer(res$indices_0based) + 1L,
                  coeff = as.numeric(res$coefficients),
                  max_iter_reached = isTRUE(res$max_iter_reached)))
    }
  }

  if (!exists(".omp_encode_R", inherits = FALSE)) {
    stop("omp_encode_rcpp unavailable and R fallback missing", call. = FALSE)
  }

  .omp_encode_R(y, D, tol = tol, max_nonzero = max_nonzero)
}

.quantize_weights <- function(w, bits = 5L) {
  qmax <- 2^(bits - 1) - 1
  max_val <- max(abs(w))
  scale <- if (max_val > 0) max_val / qmax else 1
  q <- as.integer(round(w / scale))
  q[q > qmax] <- qmax
  q[q < -qmax] <- -qmax
  list(q = q, scale = scale)
}

#' Default parameters for empirical HRBF compression
#' @export
#' @keywords internal
lna_default.basis.empirical_hrbf_compressed <- function() {
  list(svd_rank = 120L, omp_tol = 0.01, omp_sparsity_limit = 32L,
       omp_quant_bits = 5L)
}

