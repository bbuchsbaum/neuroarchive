identity_affine <- function() {
  matrix(
    c(1.75, 0.25, 0, 0,
      -0.5, 2.25, 0.125, 0,
      0, 0.375, 3.5, 0,
      -12, 8, 24, 1),
    nrow = 4L,
    ncol = 4L
  )
}

test_that("spatial identity golden digests are stable", {
  affine <- identity_affine()
  grid_id <- .lna_grid_id(
    c(5, 6, 7), affine, "mm", "first-axis-fastest"
  )
  support <- rep(c(TRUE, FALSE, TRUE), length.out = 5 * 6 * 7)

  expect_identical(
    .lna_affine_id(affine),
    "sha256:1a3b185fed618effe84064eb05777a86caca2451206cc0e61dc7b361960ecf42"
  )
  expect_identical(
    grid_id,
    "sha256:6e996d3cb42b9f82cce8a190e945d2092e29572172521a20d0791b72007b5c72"
  )
  expect_identical(
    .lna_support_id(grid_id, support),
    "sha256:bd49607ce863da8357da7518babdd86218de5ea64a0d8f2c06d8a9ab4f539fc2"
  )
})

test_that("public spatial identity constructors match canonical internals", {
  affine <- identity_affine()
  grid_id <- lna_grid_identity(c(5, 6, 7), affine)
  support <- rep(c(TRUE, FALSE), length.out = 5 * 6 * 7)

  expect_identical(lna_affine_identity(affine), .lna_affine_id(affine))
  expect_identical(
    grid_id,
    .lna_grid_id(c(5, 6, 7), affine, "mm", "first-axis-fastest")
  )
  expect_identical(
    lna_support_identity(grid_id, support),
    .lna_support_id(grid_id, support)
  )
})

test_that("equivalent R containers produce identical spatial identities", {
  affine <- identity_affine()
  affine_frame <- as.data.frame(affine)
  dims <- c(x = 5L, y = 6L, z = 7L)
  grid_matrix <- .lna_grid_id(dims, affine, "mm", "first-axis-fastest")
  grid_frame <- .lna_grid_id(
    unname(as.double(dims)), affine_frame, "mm", "first-axis-fastest"
  )
  support_array <- array(rep(c(1L, 0L), 105), dim = c(5, 6, 7))

  expect_identical(.lna_affine_id(affine), .lna_affine_id(as.vector(affine)))
  expect_identical(.lna_affine_id(affine), .lna_affine_id(affine_frame))
  expect_identical(grid_matrix, grid_frame)
  expect_identical(
    .lna_support_id(grid_matrix, support_array),
    .lna_support_id(grid_matrix, as.logical(as.vector(support_array)))
  )
})

test_that("affine identity changes for rotation and shear", {
  base <- diag(4)
  rotated <- base
  theta <- pi / 8
  rotated[1:2, 1:2] <- matrix(
    c(cos(theta), -sin(theta), sin(theta), cos(theta)),
    nrow = 2L,
    byrow = TRUE
  )
  sheared <- base
  sheared[1, 2] <- 0.2

  expect_false(identical(.lna_affine_id(base), .lna_affine_id(rotated)))
  expect_false(identical(.lna_affine_id(base), .lna_affine_id(sheared)))
  expect_false(identical(.lna_affine_id(rotated), .lna_affine_id(sheared)))
})

test_that("each grid field contributes to grid identity", {
  affine <- identity_affine()
  baseline <- .lna_grid_id(c(5, 6, 7), affine, "mm", "first-axis-fastest")
  changed_affine <- affine
  changed_affine[1, 4] <- changed_affine[1, 4] + 1

  alternatives <- c(
    .lna_grid_id(c(5, 6, 8), affine, "mm", "first-axis-fastest"),
    .lna_grid_id(c(5, 6, 7), changed_affine, "mm", "first-axis-fastest"),
    .lna_grid_id(c(5, 6, 7), affine, "m", "first-axis-fastest"),
    .lna_grid_id(c(5, 6, 7), affine, "mm", "last-axis-fastest")
  )

  expect_true(all(alternatives != baseline))
  expect_length(unique(c(baseline, alternatives)), 5L)
})

test_that("grid and support changes contribute to support identity", {
  affine <- identity_affine()
  grid <- .lna_grid_id(c(2, 2, 2), affine, "mm", "first-axis-fastest")
  other_grid <- .lna_grid_id(c(2, 2, 2), affine, "cm", "first-axis-fastest")
  empty <- rep(FALSE, 8)
  full <- rep(TRUE, 8)
  masked <- c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE)

  ids <- c(
    .lna_support_id(grid, empty),
    .lna_support_id(grid, full),
    .lna_support_id(grid, masked),
    .lna_support_id(other_grid, masked)
  )
  expect_length(unique(ids), 4L)
})

test_that("legacy parcel affine digest is preserved by the canonical API", {
  affine <- identity_affine()
  expect_identical(
    .lna_affine_id(affine),
    .parcel_mean_numeric_checksum(affine, "lna_affine_v1")
  )
})

test_that("spatial identity validation rejects ambiguous inputs", {
  expect_error(.lna_affine_id(1:15), class = "lna_error_spatial_identity")
  expect_error(
    .lna_grid_id(c(2, 2.5, 2), diag(4), "mm", "first-axis-fastest"),
    class = "lna_error_spatial_identity"
  )
  expect_error(
    .lna_support_id("not-a-grid-id", c(TRUE, FALSE)),
    class = "lna_error_spatial_identity"
  )
  expect_error(
    .lna_support_id(paste0("sha256:", strrep("a", 64)), c(0, 2)),
    class = "lna_error_spatial_identity"
  )
})
