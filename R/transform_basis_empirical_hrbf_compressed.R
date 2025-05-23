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
  if (requireNamespace("irlba", quietly = TRUE)) {
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
                       scale = q$scale)
  }

  plan <- handle$plan
  fname <- plan$get_next_filename(type)
  base <- tools::file_path_sans_ext(fname)
  vt_path <- paste0("/basis/", base, "/vt_matrix")
  codes_path <- paste0("/basis/", base, "/hrbf_codes")
  params_json <- as.character(jsonlite::toJSON(p, auto_unbox = TRUE))

  desc$version <- "1.0"
  desc$inputs <- c(input_key)
  desc$outputs <- character()
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
  kernel_type <- p$kernel_type %||% "gaussian"
  seed <- p$seed %||% 1L

  voxel_to_world <- function(vox_mat) {
    spc <- tryCatch(space(mask_neurovol), error = function(e) NULL)
    spacing_vec <- tryCatch(spacing(spc), error = function(e) c(1,1,1))
    origin_vec <- tryCatch(origin(spc), error = function(e) c(0,0,0))
    sweep(vox_mat - 1, 2, spacing_vec, `*`) +
      matrix(origin_vec, nrow(vox_mat), 3, byrow = TRUE)
  }

  centres_list <- list(); sigs <- numeric()
  for (j in seq_len(levels + 1L) - 1L) {
    sigma_j <- sigma0 / (2^j)
    r_j <- radius_factor * sigma_j
    vox_centres <- poisson_disk_sample_neuroim2(mask_neurovol, r_j, seed + j)
    if (nrow(vox_centres) > 0) {
      centres_list[[length(centres_list) + 1L]] <- voxel_to_world(vox_centres)
      sigs <- c(sigs, rep(sigma_j, nrow(vox_centres)))
    }
  }
  C_total <- if (length(centres_list) > 0) do.call(rbind, centres_list)
             else matrix(numeric(0), ncol = 3)
  sigma_vec <- sigs

  mask_arr <- as.array(mask_neurovol)
  mask_linear_indices <- which(mask_arr)
  vox_coords <- which(mask_arr, arr.ind = TRUE)
  mask_coords_world <- voxel_to_world(vox_coords)
  n_total_vox <- length(mask_arr)
  k_actual <- nrow(C_total)

  if (k_actual > 0) {
    i_idx <- integer(); j_idx <- integer(); x_val <- numeric()
    for (kk in seq_len(k_actual)) {
      atom <- generate_hrbf_atom(mask_coords_world, mask_linear_indices,
                                 C_total[kk,], sigma_vec[kk], kernel_type)
      i_idx <- c(i_idx, rep.int(kk, length(atom$indices)))
      j_idx <- c(j_idx, atom$indices)
      x_val <- c(x_val, atom$values)
    }
    Matrix::sparseMatrix(i = i_idx, j = j_idx, x = x_val,
                         dims = c(k_actual, n_total_vox))
  } else {
    Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(),
                         dims = c(0, n_total_vox))
  }
}

.omp_encode <- function(y, D, tol = 1e-3, max_nonzero = 32L) {
  r <- y
  idx <- integer(); coeff <- numeric()
  while (sum(r^2) > tol^2 && length(idx) < max_nonzero) {
    corr <- as.numeric(crossprod(D, r))
    j <- which.max(abs(corr))
    idx <- unique(c(idx, j))
    D_sub <- D[, idx, drop = FALSE]
    coeff <- as.numeric(qr.solve(D_sub, y))
    r <- y - D_sub %*% coeff
  }
  list(idx = idx, coeff = coeff)
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
