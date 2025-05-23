#' Label connected components in a logical 3D array
#'
#' Provides a minimal 6-neighbourhood connected component labeller used by the
#' HRBF helpers. When available and enabled, a faster Rcpp implementation is
#' utilised. Otherwise a pure-R fallback is employed.
#'
#' @param mask_arr_3d Logical array representing the mask.
#' @return A list with elements `count` (number of components) and `labels`
#'   (integer array of the same dimensions as `mask_arr_3d`).
#' @keywords internal
label_components <- function(mask_arr_3d) {
  dims <- dim(mask_arr_3d)
  use_rcpp <- isTRUE(getOption("lna.hrbf.use_rcpp_helpers", TRUE)) &&
    exists("label_components_6N_rcpp")

  if (use_rcpp) {
    labels_flat <- tryCatch(
      label_components_6N_rcpp(as.logical(mask_arr_3d), dims),
      error = function(e) NULL
    )
    if (!is.null(labels_flat)) {
      label_arr <- array(labels_flat, dim = dims)
      return(list(count = max(labels_flat), labels = label_arr))
    }
  }

  visited <- array(FALSE, dim = dims)
  labels <- array(0L, dim = dims)
  comp_id <- 0L
  neighbours <- matrix(c(1,0,0,-1,0,0,0,1,0,0,-1,0,0,0,1,0,0,-1),
                       ncol = 3, byrow = TRUE)

  for (i in seq_len(dims[1])) {
    for (j in seq_len(dims[2])) {
      for (k in seq_len(dims[3])) {
        if (mask_arr_3d[i, j, k] && !visited[i, j, k]) {
          comp_id <- comp_id + 1L
          q <- list(c(i, j, k))
          while (length(q) > 0) {
            pt <- q[[1]]
            q <- q[-1]
            ii <- pt[1]; jj <- pt[2]; kk <- pt[3]
            if (visited[ii, jj, kk]) next
            visited[ii, jj, kk] <- TRUE
            labels[ii, jj, kk] <- comp_id
            for (n in seq_len(nrow(neighbours))) {
              nn <- pt + neighbours[n, ]
              ni <- nn[1]; nj <- nn[2]; nk <- nn[3]
              if (ni >= 1 && ni <= dims[1] &&
                  nj >= 1 && nj <= dims[2] &&
                  nk >= 1 && nk <= dims[3] &&
                  mask_arr_3d[ni, nj, nk] && !visited[ni, nj, nk]) {
                q[[length(q) + 1L]] <- c(ni, nj, nk)
              }
            }
          }
        }
      }
    }
  }
  list(count = comp_id, labels = labels)
}

#' Poisson-disk sampling for LogicalNeuroVol masks
#'
#' @description Internal helper implementing a basic Poisson-disk sampler for
#' `neuroim2::LogicalNeuroVol` objects. The algorithm works in voxel space and
#' uses a simple sequential rejection scheme. When called on a mask with
#' multiple disconnected components, sampling is performed independently per
#' component. Components are identified using a minimal connectivity routine and
#' reseeded with `seed + component_id` so results remain deterministic.
#'
#' @param mask_neurovol A `LogicalNeuroVol` mask.
#' @param radius_mm Sampling radius in millimetres.
#' @param seed Integer RNG seed.
#' @param component_id_for_seed_offset Integer offset added to the seed when the
#'   mask represents a single component. Users typically call this function on
#'   the full mask with the default `0`.
#'
#' @return Integer matrix with columns `i`, `j`, `k` containing voxel
#'   coordinates of sampled centres.
#' @keywords internal
poisson_disk_sample_neuroim2 <- function(mask_neurovol, radius_mm, seed,
                                         component_id_for_seed_offset = 0) {
  if (!inherits(mask_neurovol, "LogicalNeuroVol")) {
    abort_lna("mask_neurovol must be LogicalNeuroVol",
              .subclass = "lna_error_validation",
              location = "poisson_disk_sample_neuroim2")
  }

  mask_arr <- as.array(mask_neurovol)
  spc <- tryCatch(space(mask_neurovol), error = function(e) NULL)
  spacing_vec <- tryCatch(spacing(spc), error = function(e) c(1, 1, 1))

  comp_info <- label_components(mask_arr)
  use_rcpp <- isTRUE(getOption("lna.hrbf.use_rcpp_helpers", TRUE)) &&
    exists("poisson_disk_sample_component_rcpp")

  radius_vox <- radius_mm / mean(spacing_vec)
  r2 <- radius_vox^2

  sample_component_R <- function(coords, base_seed) {
    set.seed(as.integer(base_seed))
    remaining <- coords[sample(nrow(coords)), , drop = FALSE]
    selected <- matrix(numeric(0), ncol = 3)
    while (nrow(remaining) > 0) {
      cand <- remaining[1, , drop = FALSE]
      remaining <- remaining[-1, , drop = FALSE]
      if (nrow(selected) == 0) {
        selected <- rbind(selected, cand)
      } else {
        d2 <- rowSums((selected - matrix(cand, nrow = nrow(selected), ncol = 3,
                                          byrow = TRUE))^2)
        if (all(d2 >= r2)) {
          selected <- rbind(selected, cand)
        }
      }
    }
    selected
  }

  gather <- function(coords, seed_val) {
    if (use_rcpp) {
      out <- poisson_disk_sample_component_rcpp(coords - 1L, r2, seed_val) + 1L
    } else {
      out <- sample_component_R(coords, seed_val)
    }
    if (nrow(out) == 0 && nrow(coords) < 150) {
      centroid_pt <- round(colMeans(coords))
      out <- matrix(centroid_pt, nrow = 1)
    }
    out
  }

  if (comp_info$count > 1L && component_id_for_seed_offset == 0) {
    centres <- lapply(seq_len(comp_info$count), function(id) {
      coords <- which(comp_info$labels == id, arr.ind = TRUE)
      if (nrow(coords) == 0) return(matrix(numeric(0), ncol = 3))
      gather(coords, as.integer(seed) + id)
    })
    res <- do.call(rbind, centres)
    colnames(res) <- c("i", "j", "k")
    return(res)
  }

  vox_coords <- which(mask_arr, arr.ind = TRUE)
  if (nrow(vox_coords) == 0) {
    return(matrix(integer(0), ncol = 3, dimnames = list(NULL, c("i","j","k"))))
  }

  selected <- gather(vox_coords, as.integer(seed) + as.integer(component_id_for_seed_offset))
  colnames(selected) <- c("i", "j", "k")
  selected
}

#' Generate an analytic HRBF atom over a mask
#'
#' @description Internal helper that evaluates a radial basis function centred at
#' `centre_coord_world` on all voxels of a mask. The mask is provided via its
#' voxel world coordinates and corresponding linear indices. Optionally the atom
#' is \eqn{L_2}-normalised over the mask voxels.
#'
#' @param mask_coords_world Numeric matrix of world coordinates for the mask
#'   voxels (\eqn{N_{maskvox} \times 3}).
#' @param mask_linear_indices Integer vector of the same length giving the voxel
#'   linear indices within the full volume.
#' @param centre_coord_world Numeric vector of length 3 giving the RBF centre in
#'   world coordinates.
#' @param sigma_mm Numeric width parameter in millimetres.
#' @param kernel_type Either \code{"gaussian"} or \code{"wendland_c4"}.
#' @param normalize_over_mask Logical; if \code{TRUE} the returned values are
#'   normalised to unit \eqn{L_2} norm over the mask.
#'
#' @return A list with elements \code{values} (numeric vector of length
#'   \code{nrow(mask_coords_world)}) and \code{indices}
#'   (\code{mask_linear_indices}).
#' @keywords internal
generate_hrbf_atom <- function(mask_coords_world, mask_linear_indices,
                               centre_coord_world, sigma_mm,
                               kernel_type = c("gaussian", "wendland_c4"),
                               normalize_over_mask = TRUE) {
  mask_coords_world <- as.matrix(mask_coords_world)
  if (ncol(mask_coords_world) != 3) {
    abort_lna("mask_coords_world must have 3 columns",
              .subclass = "lna_error_validation",
              location = "generate_hrbf_atom")
  }
  if (length(mask_linear_indices) != nrow(mask_coords_world)) {
    abort_lna("mask_linear_indices length mismatch",
              .subclass = "lna_error_validation",
              location = "generate_hrbf_atom")
  }
  centre_coord_world <- as.numeric(centre_coord_world)
  if (length(centre_coord_world) != 3) {
    abort_lna("centre_coord_world must be length 3",
              .subclass = "lna_error_validation",
              location = "generate_hrbf_atom")
  }
  kernel_type <- match.arg(kernel_type)

  diffs <- sweep(mask_coords_world, 2, centre_coord_world, FUN = "-")
  dist_mm <- sqrt(rowSums(diffs^2))

  if (kernel_type == "gaussian") {
    phi <- exp(-(dist_mm^2) / (2 * sigma_mm^2))
  } else { # wendland_c4
    r <- dist_mm / sigma_mm
    base <- pmax(0, 1 - r)
    phi <- base^8 * (32 * r^3 + 25 * r^2 + 8 * r + 1)
    phi[r >= 1] <- 0
  }

  if (normalize_over_mask) {
    norm_val <- sqrt(sum(phi^2))
    if (norm_val > 0) {
      phi <- phi / norm_val
    }
  }

  list(values = phi, indices = mask_linear_indices)
}

#' Regenerate an analytic HRBF basis matrix from descriptor parameters
#'
#' @param params List of HRBF parameters (as from a descriptor).
#' @param mask_neurovol `LogicalNeuroVol` mask defining the voxel grid.
#' @param h5_root Optional H5 group if centres are stored in the file.
#' @return Sparse matrix with one row per HRBF atom and columns matching
#'   mask voxels.
#' @keywords internal
hrbf_basis_from_params <- function(params, mask_neurovol, h5_root = NULL) {
  sigma0 <- params$sigma0 %||% 6
  levels <- params$levels %||% 3L
  radius_factor <- params$radius_factor %||% 2.5
  kernel_type <- params$kernel_type %||% "gaussian"
  num_extra_fine_levels <- params$num_extra_fine_levels %||% 0L
  seed <- params$seed
  centres_path <- params$centres_path
  sigma_vec_path <- params$sigma_vec_path
  centres_stored <- isTRUE(params$centres_stored)

  voxel_to_world <- function(vox_mat) {
    spc <- tryCatch(space(mask_neurovol), error = function(e) NULL)
    spacing_vec <- tryCatch(spacing(spc), error = function(e) c(1,1,1))
    origin_vec <- tryCatch(origin(spc), error = function(e) c(0,0,0))
    sweep(vox_mat - 1, 2, spacing_vec, `*`) +
      matrix(origin_vec, nrow(vox_mat), 3, byrow = TRUE)
  }

  if (centres_stored && !is.null(centres_path) && !is.null(sigma_vec_path) &&
      !is.null(h5_root)) {
    C_total <- h5_read(h5_root, centres_path)
    sigma_vec <- as.numeric(h5_read(h5_root, sigma_vec_path))
  } else if (!is.null(seed)) {
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
    if (num_extra_fine_levels > 0L) {
      for (j_extra in seq_len(num_extra_fine_levels)) {
        sigma_new <- sigma0 / (2^(levels + j_extra))
        r_new <- radius_factor * sigma_new
        vox_centres <- poisson_disk_sample_neuroim2(
          mask_neurovol, r_new, seed + levels + j_extra
        )
        if (nrow(vox_centres) > 0) {
          centres_list[[length(centres_list) + 1L]] <- voxel_to_world(vox_centres)
          sigs <- c(sigs, rep(sigma_new, nrow(vox_centres)))
        }
      }
    }
    C_total <- if (length(centres_list) > 0) do.call(rbind, centres_list)
               else matrix(numeric(0), ncol = 3)
    sigma_vec <- sigs
  } else {
    abort_lna("Insufficient parameters to regenerate HRBF basis",
              .subclass = "lna_error_descriptor",
              location = "hrbf_basis_from_params")
  }

  mask_arr <- as.array(mask_neurovol)
  mask_coords_vox <- which(mask_arr, arr.ind = TRUE)
  mask_coords_world <- voxel_to_world(mask_coords_vox)
  mask_linear_indices <- as.integer(which(mask_arr))
  n_total_vox <- length(mask_arr)
  k_actual <- nrow(C_total)

  if (k_actual > 0) {
    triplet_i_list <- vector("list", k_actual)
    triplet_j_list <- vector("list", k_actual)
    triplet_x_list <- vector("list", k_actual)
    for (kk in seq_len(k_actual)) {
      atom <- generate_hrbf_atom(mask_coords_world, mask_linear_indices,
                                 C_total[kk, ], sigma_vec[kk], kernel_type)
      triplet_i_list[[kk]] <- rep.int(kk, length(atom$indices))
      triplet_j_list[[kk]] <- atom$indices
      triplet_x_list[[kk]] <- atom$values
    }
    i_idx <- unlist(triplet_i_list, use.names = FALSE)
    j_idx <- unlist(triplet_j_list, use.names = FALSE)
    x_val <- unlist(triplet_x_list, use.names = FALSE)
    Matrix::sparseMatrix(i = i_idx, j = j_idx, x = x_val,
                         dims = c(k_actual, n_total_vox))
  } else {
    Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(),
                         dims = c(0, n_total_vox))
  }
}

