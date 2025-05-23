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

  # Helper to label connected components (6-neighbourhood)
  label_components <- function(arr) {
    dims <- dim(arr)
    visited <- array(FALSE, dim = dims)
    labels <- array(0L, dim = dims)
    comp_id <- 0L
    neighbours <- matrix(c(1,0,0,-1,0,0,0,1,0,0,-1,0,0,0,1,0,0,-1),
                         ncol = 3, byrow = TRUE)

    for (i in seq_len(dims[1])) {
      for (j in seq_len(dims[2])) {
        for (k in seq_len(dims[3])) {
          if (arr[i, j, k] && !visited[i, j, k]) {
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
                    arr[ni, nj, nk] && !visited[ni, nj, nk]) {
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

  comp_info <- label_components(mask_arr)

  if (comp_info$count > 1L && component_id_for_seed_offset == 0) {
    centres <- lapply(seq_len(comp_info$count), function(id) {
      sub_mask <- comp_info$labels == id
      sub_vol <- structure(list(arr = sub_mask), class = "LogicalNeuroVol")
      if (!is.null(spc)) attr(sub_vol, "space") <- spc
      poisson_disk_sample_neuroim2(sub_vol, radius_mm, seed, id)
    })
    res <- do.call(rbind, centres)
    colnames(res) <- c("i", "j", "k")
    return(res)
  }

  set.seed(as.integer(seed) + as.integer(component_id_for_seed_offset))

  radius_vox <- radius_mm / mean(spacing_vec)
  r2 <- radius_vox^2

  vox_coords <- which(mask_arr, arr.ind = TRUE)
  if (nrow(vox_coords) == 0) {
    return(matrix(integer(0), ncol = 3, dimnames = list(NULL, c("i","j","k"))))
  }
  remaining <- vox_coords[sample(nrow(vox_coords)), , drop = FALSE]
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

  if (nrow(selected) == 0 && sum(mask_arr) < 150) {
    centroid_pt <- round(colMeans(vox_coords))
    selected <- matrix(centroid_pt, nrow = 1)
  }

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
