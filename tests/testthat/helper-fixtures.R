# Utility functions for Haar octwave tests

make_ball_mask <- function(radius = 2L) {
  dims <- rep(radius*2L + 1L, 3)
  centre <- radius + 1L
  mask <- array(FALSE, dim = dims)
  for (i in seq_len(dims[1])) {
    for (j in seq_len(dims[2])) {
      for (k in seq_len(dims[3])) {
        if (sqrt((i-centre)^2 + (j-centre)^2 + (k-centre)^2) <= radius) {
          mask[i,j,k] <- TRUE
        }
      }
    }
  }
  mask
}

make_disjoint_boxes_mask <- function() {
  mask <- array(FALSE, dim = c(5,5,5))
  mask[1:2,1:2,1:2] <- TRUE
  mask[4:5,4:5,4:5] <- TRUE
  mask
}

make_synthetic_fmri <- function(time = 10L, mask) {
  nvox <- sum(mask)
  array(rnorm(time * nvox), dim = c(dim(mask), time))
}
