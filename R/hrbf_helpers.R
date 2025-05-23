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
