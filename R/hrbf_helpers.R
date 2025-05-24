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
#' @param edge_binary_map Optional logical array the same dimensions as the mask
#'   indicating edge voxels for adaptive sampling.
#' @param density_factor Numeric factor (>1) to reduce radius when sampling
#'   candidate centres in edge regions.
#'
#' @return Integer matrix with columns `i`, `j`, `k` containing voxel
#'   coordinates of sampled centres.
#' @keywords internal
poisson_disk_sample_neuroim2 <- function(mask_neurovol, radius_mm, seed,
                                         component_id_for_seed_offset = 0,
                                         edge_binary_map = NULL,
                                         density_factor = 1.5) {
  if (!inherits(mask_neurovol, "LogicalNeuroVol")) {
    abort_lna("mask_neurovol must be LogicalNeuroVol",
              .subclass = "lna_error_validation",
              location = "poisson_disk_sample_neuroim2")
  }

  mask_arr <- as.array(mask_neurovol)
  spc <- tryCatch(space(mask_neurovol), error = function(e) NULL)
  spacing_vec <- tryCatch(spacing(spc), error = function(e) c(1, 1, 1))

  if (!is.null(edge_binary_map)) {
    if (!is.logical(edge_binary_map) || !all(dim(edge_binary_map) == dim(mask_arr))) {
      abort_lna("edge_binary_map dims mismatch mask",
                .subclass = "lna_error_validation",
                location = "poisson_disk_sample_neuroim2")
    }
  }

  comp_info <- label_components(mask_arr)
  use_rcpp <- isTRUE(getOption("lna.hrbf.use_rcpp_helpers", TRUE)) &&
    exists("poisson_disk_sample_component_rcpp") && is.null(edge_binary_map)

  radius_vox <- radius_mm / mean(spacing_vec)
  r2 <- radius_vox^2
  r2_edge <- (radius_vox / density_factor)^2

  sample_component_R <- function(coords, base_seed) {
    set.seed(as.integer(base_seed))
    remaining <- coords[sample(nrow(coords)), , drop = FALSE]
    selected <- matrix(numeric(0), ncol = 3)
    while (nrow(remaining) > 0) {
      cand <- remaining[1, , drop = FALSE]
      remaining <- remaining[-1, , drop = FALSE]
      cand_edge <- !is.null(edge_binary_map) && edge_binary_map[cand]
      r2_eff <- if (cand_edge) r2_edge else r2
      if (nrow(selected) == 0) {
        selected <- rbind(selected, cand)
      } else {
        d2 <- rowSums((selected - matrix(cand, nrow = nrow(selected), ncol = 3,
                                          byrow = TRUE))^2)
        if (all(d2 >= r2_eff)) {
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
#' @param current_level_j Integer level index of this atom (0-indexed).
#' @param total_levels Total number of levels in the pyramid (\eqn{J}).
#' @param params Parameter list defining \code{kernel_type},
#'   \code{kernel_type_fine_levels} and \code{num_fine_levels_alt_kernel}.
#' @param normalize_over_mask Logical; if \code{TRUE} the returned values are
#'   normalised to unit \eqn{L_2} norm over the mask.
#'
#' @return A list with elements \code{values} (numeric vector of length
#'   \code{nrow(mask_coords_world)}) and \code{indices}
#'   (\code{mask_linear_indices}).
#' @keywords internal
generate_hrbf_atom <- function(mask_coords_world, mask_linear_indices,
                               centre_coord_world, sigma_mm,
                               current_level_j, total_levels, params,
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
  p_kernel_type <- params$kernel_type %||% "gaussian"
  p_kernel_type_fine <- params$kernel_type_fine_levels %||% "wendland_c4"
  num_alt <- params$num_fine_levels_alt_kernel %||% 0L
  use_alt_kernel <- current_level_j > (total_levels - num_alt)
  eff_kernel <- if (use_alt_kernel) p_kernel_type_fine else p_kernel_type

  diffs <- sweep(mask_coords_world, 2, centre_coord_world, FUN = "-")
  dist_mm <- sqrt(rowSums(diffs^2))

  if (eff_kernel == "gaussian") {
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
#'
#' The parameter list may include `num_extra_fine_levels` to generate
#' additional dyadic levels beyond `levels`. Each extra level halves
#' `sigma` and samples a denser set of centres.
#' @keywords internal
hrbf_basis_from_params <- function(params, mask_neurovol, h5_root = NULL,
                                   mask_world_coords = NULL, mask_arr = NULL,
                                   mask_linear_indices = NULL) {
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
    level_vec <- as.integer(round(log2(sigma0 / sigma_vec)))
  } else if (!is.null(seed)) {
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
  } else {
    abort_lna("Insufficient parameters to regenerate HRBF basis",
              .subclass = "lna_error_descriptor",
              location = "hrbf_basis_from_params")
  }

  if (is.null(mask_arr)) {
    mask_arr <- as.array(mask_neurovol)
  }
  if (is.null(mask_world_coords)) {
    mask_coords_vox <- which(mask_arr, arr.ind = TRUE)
    mask_coords_world <- voxel_to_world(mask_coords_vox)
  } else {
    mask_coords_world <- mask_world_coords
  }
  if (is.null(mask_linear_indices)) {
    mask_linear_indices <- as.integer(which(mask_arr))
  }
  n_total_vox <- length(mask_arr)
  k_actual <- nrow(C_total)

  if (k_actual > 0) {
    triplet_i_list <- vector("list", k_actual)
    triplet_j_list <- vector("list", k_actual)
    triplet_x_list <- vector("list", k_actual)
    for (kk in seq_len(k_actual)) {
      atom <- generate_hrbf_atom(mask_coords_world, mask_linear_indices,
                                 C_total[kk, ], sigma_vec[kk],
                                 level_vec[kk], levels, params)
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


#' Compute edge binary map for neuroim2 volumes
#'
#' Helper used for edge-adaptive HRBF sampling. Depending on
#' `source_spec`, either computes a Sobel gradient magnitude from the
#' current data or loads a pre-computed gradient map.
#'
#' @details 
#' **Performance Warning:** The `self_mean` mode computes 3D Sobel gradients 
#' using pure R, which is extremely slow on large volumes (≈10 min for 64³). 
#' For performance-critical applications, consider:
#' \itemize{
#'   \item Pre-computing structural gradients and using `structural_path` mode
#'   \item Installing the Rcpp acceleration (when available): 30-100× faster
#'   \item Disabling edge-adaptive sampling via parameters
#' }
#'
#' If `structural_to_epi_affine_path` is supplied and the structural
#' gradient map's dimensions do not match the mask, the map is
#' resampled using `neuroim2` with that affine. When dimensions differ
#' and no affine is provided, a `lna_error_validation` is thrown.
#'
#' @param source_spec Character string specifying the source of the
#'   gradient map ("self_mean" or "structural_path").
#' @param data_handle A `DataHandle` providing access to the input data
#'   and HDF5 file if required.
#' @param params_edge_adaptive List of edge adaptive parameters. Fields
#'   `structural_path`, `structural_to_epi_affine_path`, and
#'   `edge_thresh_k` are used.
#' @return Logical 3D array indicating edge voxels.
#' @keywords internal
compute_edge_map_neuroim2 <- function(source_spec, data_handle,
                                      params_edge_adaptive) {
  mask_nv <- data_handle$mask_info$mask
  if (is.null(mask_nv)) {
    abort_lna("mask_info$mask missing", .subclass = "lna_error_validation",
              location = "compute_edge_map_neuroim2")
  }
  mask_arr <- as.array(mask_nv)
  dims <- dim(mask_arr)
  p <- params_edge_adaptive %||% list()
  thresh_k <- p$edge_thresh_k %||% 3.0

  sobel_mag <- function(vol) {
    # Try fast Rcpp implementation first
    use_rcpp <- isTRUE(getOption("lna.edge_adaptive.use_rcpp", TRUE))
    if (use_rcpp) {
      # Check if function is available in the package DLL
      dll_routines <- tryCatch(
        getDLLRegisteredRoutines("neuroarchive"),
        error = function(e) NULL
      )
      use_rcpp <- !is.null(dll_routines) && 
        "sobel3d_magnitude_rcpp" %in% names(dll_routines$`.Call`)
    }
    
    if (use_rcpp) {
      result <- tryCatch(
        sobel3d_magnitude_rcpp(vol),
        error = function(e) NULL
      )
      if (!is.null(result)) {
        return(result)
      }
    }
    
    # Performance warning for large volumes
    nvox <- prod(dim(vol))
    if (nvox > 50000) {  # Warn for volumes larger than ~37³
      warn_lna(
        sprintf("Computing 3D Sobel gradients on %dx%dx%d volume in pure R. This may take several minutes. Consider pre-computing gradients or installing Rcpp acceleration.",
                dim(vol)[1], dim(vol)[2], dim(vol)[3]),
        .subclass = "lna_warning_performance",
        location = "compute_edge_map_neuroim2:sobel_mag"
      )
    }
    
    # Pure R fallback implementation
    w <- matrix(c(1,2,1,2,4,2,1,2,1), nrow = 3, byrow = TRUE)
    kx <- array(0, c(3,3,3)); ky <- array(0, c(3,3,3)); kz <- array(0, c(3,3,3))
    for (i in 1:3) for (j in 1:3) {
      kx[1,i,j] <- -w[i,j]; kx[3,i,j] <- w[i,j]
      ky[i,1,j] <- -w[i,j]; ky[i,3,j] <- w[i,j]
      kz[i,j,1] <- -w[i,j]; kz[i,j,3] <- w[i,j]
    }
    conv3d <- function(arr, ker) {
      d <- dim(arr); out <- array(0, d)
      for (x in 2:(d[1]-1)) for (y in 2:(d[2]-1)) for (z in 2:(d[3]-1)) {
        sub <- arr[(x-1):(x+1), (y-1):(y+1), (z-1):(z+1)]
        out[x,y,z] <- sum(sub * ker)
      }
      out
    }
    gx <- conv3d(vol, kx); gy <- conv3d(vol, ky); gz <- conv3d(vol, kz)
    sqrt(gx^2 + gy^2 + gz^2)
  }

  if (identical(source_spec, "self_mean")) {
    inp <- data_handle$pull_first(c("input_dense_mat", "dense_mat", "input"))
    X <- as_dense_mat(inp$value)
    mean_vec <- colMeans(X)
    vol <- array(mean_vec, dim = dims)
    grad_map <- sobel_mag(vol)
  } else if (identical(source_spec, "structural_path")) {
    if (is.null(data_handle$h5)) {
      abort_lna("H5 handle required for structural_path",
                .subclass = "lna_error_validation",
                location = "compute_edge_map_neuroim2")
    }
    root <- data_handle$h5[["/"]]
    grad_map <- h5_read(root, p$structural_path)
    if (!all(dim(grad_map) == dims)) {
      if (!is.null(p$structural_to_epi_affine_path)) {
        affine <- h5_read(root, p$structural_to_epi_affine_path)
        if (!is.matrix(affine) || !all(dim(affine) == c(4L, 4L))) {
          abort_lna("Invalid structural_to_epi affine",
                    .subclass = "lna_error_validation",
                    location = "compute_edge_map_neuroim2")
        }

        safe_call <- function(fn_name, ...) {
          if (exists(fn_name, envir = .GlobalEnv, mode = "function")) {
            get(fn_name, envir = .GlobalEnv)(...)
          } else if (exists(fn_name,
                            envir = asNamespace("neuroim2"), mode = "function")) {
            get(fn_name, envir = asNamespace("neuroim2"))(...)
          } else {
            stop(sprintf("Function %s not found", fn_name), call. = FALSE)
          }
        }

        resampled <- safe_call("resample", grad_map, mask_nv, affine)
        grad_map <- if (inherits(resampled, "NeuroObj")) {
          safe_call("as.array", resampled)
        } else {
          resampled
        }
        if (!all(dim(grad_map) == dims)) {
          abort_lna("Resampled gradient map dims mismatch mask",
                    .subclass = "lna_error_validation",
                    location = "compute_edge_map_neuroim2")
        }
      } else {
        abort_lna("Gradient map dims mismatch mask",
                  .subclass = "lna_error_validation",
                  location = "compute_edge_map_neuroim2")
      }
    }
  } else {
    abort_lna("Unknown edge map source", .subclass = "lna_error_validation",
              location = "compute_edge_map_neuroim2")
  }

  med_val <- median(abs(as.numeric(grad_map[mask_arr])), na.rm = TRUE)
  edge_binary <- grad_map > (thresh_k * med_val)
  array(as.logical(edge_binary), dim = dims)
}

#' Check edge-adaptive HRBF sampling availability and performance
#'
#' Helper to determine if edge-adaptive sampling should be used based on
#' volume size and available acceleration.
#'
#' @param dims Vector of 3D dimensions
#' @param warn_large Logical, whether to warn about large volumes
#' @return List with recommendations and performance info
#' @export
check_edge_adaptive_performance <- function(dims, warn_large = TRUE) {
  nvox <- prod(dims)
  # Check if Rcpp implementation is available
  dll_routines <- tryCatch(
    getDLLRegisteredRoutines("neuroarchive"),
    error = function(e) NULL
  )
  has_rcpp <- !is.null(dll_routines) && 
    "sobel3d_magnitude_rcpp" %in% names(dll_routines$`.Call`)
  
  # Performance thresholds
  small_vol <- nvox <= 10000    # ~21³, fast even in R
  medium_vol <- nvox <= 100000  # ~46³, manageable in R
  large_vol <- nvox > 100000    # >46³, very slow in R
  
  recommendation <- if (has_rcpp) {
    "fast" # Can handle any size with Rcpp
  } else if (small_vol) {
    "acceptable" # Small enough for R
  } else if (medium_vol) {
    "slow" # Will be slow but doable
  } else {
    "disable" # Too slow, recommend disabling
  }
  
  estimated_time_r <- if (nvox <= 1000) {
    "< 1 sec"
  } else if (nvox <= 10000) {
    "< 10 sec" 
  } else if (nvox <= 50000) {
    "< 1 min"
  } else if (nvox <= 100000) {
    "1-5 min"
  } else {
    "> 5 min"
  }
  
  if (warn_large && recommendation %in% c("slow", "disable")) {
    if (has_rcpp) {
      message(sprintf("Volume %dx%dx%d is large but Rcpp acceleration available", 
                     dims[1], dims[2], dims[3]))
    } else {
      msg <- sprintf(paste0(
        "Volume %dx%dx%d will be very slow for edge-adaptive sampling (est. %s). ",
        "Consider: (1) Installing Rcpp acceleration, (2) Pre-computing structural gradients, ",
        "or (3) Disabling edge-adaptive sampling"),
        dims[1], dims[2], dims[3], estimated_time_r)
      if (recommendation == "disable") {
        warning(msg, call. = FALSE)
      } else {
        message(msg)
      }
    }
  }
  
  list(
    nvox = nvox,
    has_rcpp = has_rcpp,
    recommendation = recommendation,
    estimated_time_r = estimated_time_r,
    dims_str = paste(dims, collapse = "×")
  )
}

