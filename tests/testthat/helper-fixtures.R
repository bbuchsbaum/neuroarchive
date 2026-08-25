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

create_valid_lna <- function(path, checksum = TRUE) {
  h5 <- neuroarchive:::open_h5(path, mode = "w")
  plan <- Plan$new()
  plan$add_descriptor("00_dummy.json", list(type = "dummy"))
  plan$add_payload("payload", matrix(1:4, nrow = 2))
  plan$add_dataset_def(
    "/scans/run-01/data",
    "data",
    "dummy",
    "run-01",
    0L,
    "{}",
    "payload",
    "eager",
    dtype = NA_character_
  )
  materialise_plan(h5, plan, checksum = if (checksum) "sha256" else "none")
}

create_schema_lna <- function(path) {
  h5 <- neuroarchive:::open_h5(path, mode = "w")
  plan <- Plan$new()
  plan$add_descriptor("00_basis.json", list(type = "basis", basis_path = "foo"))
  plan$add_descriptor("01_embed.json", list(type = "embed", basis_path = "foo"))
  materialise_plan(h5, plan, checksum = "none")
}

mutate_lna_file <- function(path, mutate) {
  h5 <- neuroarchive:::open_h5(path, mode = "r+")
  on.exit(neuroarchive:::close_h5_safely(h5), add = TRUE)
  mutate(h5)
  invisible(path)
}
test_inst_file <- function(...) {
  installed <- system.file(..., package = "neuroarchive")
  if (nzchar(installed)) return(installed)
  testthat::test_path("..", "..", "inst", ...)
}
