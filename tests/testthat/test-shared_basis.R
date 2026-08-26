# Tests for shared basis infrastructure: write, read, resolve, registry, dataset

shared_basis_test_grid <- function(mask, affine = diag(4)) {
  list(
    dims = dim(mask),
    transform = affine,
    spatial_units = "mm",
    voxel_ordering = "first-axis-fastest"
  )
}

shared_basis_test_identity <- function(mask, affine = diag(4),
                                       support_type = "masked") {
  grid <- shared_basis_test_grid(mask, affine)
  grid_id <- .lna_grid_id(
    grid$dims, grid$transform, grid$spatial_units, grid$voxel_ordering
  )
  list(
    basis_id = "fixture-basis",
    basis_ids = "fixture-basis",
    checksum = paste0("sha256:", strrep("a", 64)),
    grid_id = grid_id,
    support_id = .lna_support_id(grid_id, mask),
    support_type = support_type,
    spatial_units = grid$spatial_units,
    voxel_ordering = grid$voxel_ordering,
    grid_dims = grid$dims
  )
}

shared_basis_test_handle <- function(mask = NULL, affine = diag(4),
                                     dims = NULL) {
  if (is.null(dims)) dims <- if (is.null(mask)) c(2L, 2L, 1L) else dim(mask)
  DataHandle$new(
    initial_meta = list(header = list(
      dims = dims,
      transform = affine,
      spatial_units = "mm",
      voxel_ordering = "first-axis-fastest"
    )),
    mask_info = if (is.null(mask)) NULL else list(mask = mask)
  )
}

test_that("write_shared_basis creates artifact and registry", {
  dir <- withr::local_tempdir()
  bases_dir <- file.path(dir, "bases")

  # Create a fake template as a raw list
  set.seed(42)
  n_vox <- 100
  n_atoms <- 5
  loadings <- matrix(rnorm(n_vox * n_atoms), nrow = n_vox, ncol = n_atoms)
  mask_arr <- array(TRUE, dim = c(10, 10, 1))

  template <- list(loadings = loadings, mask = mask_arr, kind = "test_basis")

  result <- write_shared_basis(template, bases_dir, basis_id = "test-k5")

  # Check return value

  expect_true(file.exists(result$file))
  expect_match(result$checksum, "^sha256:")
  expect_equal(result$basis_id, "test-k5")
  expect_match(result$filename, "^sha256-.*\\.lna_basis\\.h5$")

  # Check registry was created
  registry_path <- file.path(bases_dir, "registry.json")
  expect_true(file.exists(registry_path))

  registry <- jsonlite::fromJSON(registry_path, simplifyVector = FALSE)
  expect_equal(registry$version, "1.0")
  expect_true("test-k5" %in% names(registry$bases))
  expect_equal(registry$bases[["test-k5"]]$checksum, result$checksum)
  expect_equal(registry$bases[["test-k5"]]$kind, "test_basis")
  expect_equal(registry$bases[["test-k5"]]$n_atoms, n_atoms)
  expect_equal(registry$bases[["test-k5"]]$n_voxels, n_vox)
})

test_that("read_shared_basis round-trips correctly", {
  dir <- withr::local_tempdir()
  bases_dir <- file.path(dir, "bases")

  set.seed(42)
  n_vox <- 50
  n_atoms <- 3
  loadings <- matrix(rnorm(n_vox * n_atoms), nrow = n_vox, ncol = n_atoms)
  mask_arr <- array(c(rep(TRUE, n_vox), rep(FALSE, 50)), dim = c(10, 10, 1))
  params <- list(k = 3, method = "test")

  template <- list(loadings = loadings, mask = mask_arr, kind = "pca",
                   params = params)
  written <- write_shared_basis(template, bases_dir, basis_id = "pca-k3")

  # Read back
  result <- read_shared_basis(written$file)

  expect_equal(result$basis_id, "pca-k3")
  expect_equal(result$basis_ids, "pca-k3")
  expect_equal(result$basis_kind, "pca")
  expect_equal(result$checksum, written$checksum)
  expect_equal(dim(result$loadings), c(n_vox, n_atoms))
  expect_equal(result$loadings, loadings, tolerance = 1e-10)
  expect_equal(dim(result$mask), c(10, 10, 1))
  expect_equal(sum(result$mask), n_vox)
  expect_equal(result$params$k, 3)
})

test_that("shared basis persists grid and support identity in immutable content", {
  dir <- withr::local_tempdir()
  mask <- array(c(TRUE, FALSE, TRUE, TRUE), dim = c(2, 2, 1))
  loadings <- diag(3)
  first <- write_shared_basis(
    list(
      loadings = loadings,
      mask = mask,
      kind = "identity-test",
      grid = shared_basis_test_grid(mask),
      support_type = "masked"
    ),
    file.path(dir, "bases"),
    "identity-one"
  )
  shifted <- diag(4)
  shifted[1, 4] <- 4
  second <- write_shared_basis(
    list(
      loadings = loadings,
      mask = mask,
      kind = "identity-test",
      grid = shared_basis_test_grid(mask, shifted),
      support_type = "masked"
    ),
    file.path(dir, "bases"),
    "identity-two"
  )

  restored <- read_shared_basis(first$file)
  expect_identical(restored$grid_id, first$grid_id)
  expect_identical(restored$support_id, first$support_id)
  expect_identical(restored$support_type, "masked")
  expect_false(identical(first$checksum, second$checksum))
})

test_that("shared basis writer rejects false full-grid declarations", {
  dir <- withr::local_tempdir()
  mask <- array(c(TRUE, FALSE, TRUE, TRUE), dim = c(2, 2, 1))
  expect_error(
    write_shared_basis(
      list(
        loadings = diag(3), mask = mask, kind = "invalid-full-grid",
        grid = shared_basis_test_grid(mask), support_type = "full-grid"
      ),
      file.path(dir, "bases"),
      "invalid-full-grid"
    ),
    "support at every grid voxel",
    class = "lna_error_spatial_identity"
  )
})

test_that("read_shared_basis detects checksum corruption", {
  dir <- withr::local_tempdir()
  bases_dir <- file.path(dir, "bases")

  loadings <- matrix(1:20, nrow = 10, ncol = 2)
  mask_arr <- array(TRUE, dim = c(10, 1, 1))
  template <- list(loadings = loadings, mask = mask_arr, kind = "test")

  written <- write_shared_basis(template, bases_dir, basis_id = "corrupt-test")

  # Tamper with the loadings in the HDF5 file
  h5 <- hdf5r::H5File$new(written$file, mode = "r+")
  ds <- h5[["loadings"]]
  ds[1:10, 1:2] <- matrix(999, nrow = 10, ncol = 2)
  h5$close_all()

  expect_error(
    read_shared_basis(written$file, verify_checksum = TRUE),
    class = "lna_error_checksum"
  )

  # Should succeed without verification
  result <- read_shared_basis(written$file, verify_checksum = FALSE)
  expect_equal(result$basis_id, "corrupt-test")
})

test_that("resolve_basis finds by locator", {
  dir <- withr::local_tempdir()
  bases_dir <- file.path(dir, "bases")
  sub_dir <- file.path(dir, "sub-01", "func")
  dir.create(sub_dir, recursive = TRUE)

  loadings <- matrix(rnorm(30), nrow = 10, ncol = 3)
  mask_arr <- array(TRUE, dim = c(10, 1, 1))
  template <- list(loadings = loadings, mask = mask_arr, kind = "test")

  written <- write_shared_basis(template, bases_dir, basis_id = "loc-test")

  # Simulate an archive in sub-01/func/
  fake_archive <- file.path(sub_dir, "sub-01_bold.lna.h5")
  file.create(fake_archive)

  locator <- file.path("../../bases", written$filename)

  resolved <- resolve_basis(
    basis_id = "loc-test",
    checksum = written$checksum,
    locator = locator,
    archive_path = fake_archive
  )

  expect_equal(
    normalizePath(resolved, winslash = "/"),
    normalizePath(written$file, winslash = "/")
  )
})

test_that("resolve_basis finds via registry", {
  dir <- withr::local_tempdir()
  bases_dir <- file.path(dir, "bases")
  sub_dir <- file.path(dir, "sub-01", "func")
  dir.create(sub_dir, recursive = TRUE)

  loadings <- matrix(rnorm(30), nrow = 10, ncol = 3)
  mask_arr <- array(TRUE, dim = c(10, 1, 1))
  template <- list(loadings = loadings, mask = mask_arr, kind = "test")

  written <- write_shared_basis(template, bases_dir, basis_id = "reg-test")

  fake_archive <- file.path(sub_dir, "sub-01_bold.lna.h5")
  file.create(fake_archive)

  # Resolve without locator — should find via registry walk-up
  resolved <- resolve_basis(
    basis_id = "reg-test",
    checksum = written$checksum,
    locator = NULL,
    archive_path = fake_archive
  )

  expect_equal(
    normalizePath(resolved, winslash = "/"),
    normalizePath(written$file, winslash = "/")
  )
})

test_that("resolve_basis fails with actionable error for missing basis", {
  expect_error(
    resolve_basis(
      basis_id = "nonexistent",
      checksum = "sha256:0000",
      archive_path = tempfile()
    ),
    class = "lna_error_basis_not_found"
  )
})

test_that("list_shared_bases returns correct data frame", {
  dir <- withr::local_tempdir()
  bases_dir <- file.path(dir, "bases")

  # Empty registry
  result <- list_shared_bases(bases_dir)
  expect_equal(nrow(result), 0)
  expect_true(all(c("basis_id", "checksum", "filename", "kind") %in% names(result)))

  # Add two bases
  m1 <- matrix(rnorm(20), 10, 2)
  m2 <- matrix(rnorm(30), 10, 3)
  mask <- array(TRUE, dim = c(10, 1, 1))
  write_shared_basis(list(loadings = m1, mask = mask, kind = "a"), bases_dir, "b1")
  write_shared_basis(list(loadings = m2, mask = mask, kind = "b"), bases_dir, "b2")

  result <- list_shared_bases(bases_dir)
  expect_equal(nrow(result), 2)
  expect_setequal(result$basis_id, c("b1", "b2"))
  expect_equal(result$n_atoms[result$basis_id == "b1"], 2L)
  expect_equal(result$n_atoms[result$basis_id == "b2"], 3L)
})

test_that("registry updates preserve existing entries", {
  dir <- withr::local_tempdir()
  bases_dir <- file.path(dir, "bases")
  mask <- array(TRUE, dim = c(5, 1, 1))

  write_shared_basis(list(loadings = matrix(1:10, 5, 2), mask = mask, kind = "a"),
                     bases_dir, "first")
  write_shared_basis(list(loadings = matrix(1:15, 5, 3), mask = mask, kind = "b"),
                     bases_dir, "second")

  registry <- jsonlite::fromJSON(
    file.path(bases_dir, "registry.json"), simplifyVector = FALSE
  )
  expect_true("first" %in% names(registry$bases))
  expect_true("second" %in% names(registry$bases))
})

test_that("same immutable artifact can have multiple aliases without mutation", {
  dir <- withr::local_tempdir()
  bases_dir <- file.path(dir, "bases")
  mask <- array(TRUE, dim = c(5, 1, 1))
  loadings <- matrix(rnorm(10), 5, 2)

  first <- write_shared_basis(
    list(loadings = loadings, mask = mask, kind = "alias-test"),
    bases_dir,
    "alias-one"
  )
  second <- write_shared_basis(
    list(loadings = loadings, mask = mask, kind = "alias-test"),
    bases_dir,
    "alias-two"
  )

  expect_equal(first$filename, second$filename)

  registry <- jsonlite::fromJSON(file.path(bases_dir, "registry.json"), simplifyVector = FALSE)
  expect_equal(registry$bases$`alias-one`$filename, first$filename)
  expect_equal(registry$bases$`alias-two`$filename, first$filename)

  artifact <- read_shared_basis(first$file)
  expect_null(artifact$basis_id)
  expect_setequal(artifact$basis_ids, c("alias-one", "alias-two"))
})

test_that("different masks produce different immutable artifacts", {
  dir <- withr::local_tempdir()
  bases_dir <- file.path(dir, "bases")
  loadings <- matrix(rnorm(12), 6, 2)
  mask1 <- array(c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE), dim = c(3, 2, 1))
  mask2 <- array(c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE), dim = c(3, 2, 1))

  first <- write_shared_basis(list(loadings = loadings, mask = mask1, kind = "mask-a"),
                              bases_dir, "mask-a")
  second <- write_shared_basis(list(loadings = loadings, mask = mask2, kind = "mask-b"),
                               bases_dir, "mask-b")

  expect_false(identical(first$filename, second$filename))
})

test_that("shared basis artifacts preserve sparse loadings", {
  skip_if_not_installed("Matrix")
  dir <- withr::local_tempdir()
  bases_dir <- file.path(dir, "bases")
  mask <- array(TRUE, dim = c(6, 1, 1))
  dense <- matrix(c(
    1, 0,
    0, 0,
    0, 2,
    3, 0,
    0, 0,
    0, 4
  ), nrow = 6, byrow = TRUE)
  loadings <- Matrix::Matrix(dense, sparse = TRUE)

  written <- write_shared_basis(
    list(loadings = loadings, mask = mask, kind = "sparse"),
    bases_dir,
    "sparse-basis"
  )
  result <- read_shared_basis(written$file)

  expect_s4_class(result$loadings, "dgCMatrix")
  expect_equal(as.matrix(result$loadings), dense)
})

test_that("embed.shared_basis round-trips through a nearby basis registry", {
  dir <- withr::local_tempdir()

  mask <- array(c(TRUE, TRUE, TRUE, TRUE), dim = c(2, 2, 1))
  loadings <- diag(4)
  written_basis <- write_shared_basis(
    list(
      loadings = loadings,
      mask = mask,
      kind = "test_embed",
      grid = shared_basis_test_grid(mask),
      support_type = "masked"
    ),
    file.path(dir, "bases"),
    "embed-test"
  )

  arr <- array(c(
    1, 2, 3, 4,
    2, 3, 4, 5,
    3, 4, 5, 6
  ), dim = c(2, 2, 1, 3))

  archive_dir <- file.path(dir, "sub-01", "func")
  dir.create(archive_dir, recursive = TRUE)
  archive <- file.path(archive_dir, "sub-01_task-rest_bold.lna.h5")
  write_lna(
    arr,
    archive,
    mask = mask,
    transforms = "embed.shared_basis",
    header = list(
      dims = dim(mask), transform = diag(4), spatial_units = "mm",
      voxel_ordering = "first-axis-fastest"
    ),
    transform_params = list(
      "embed.shared_basis" = list(
        basis_id = "embed-test",
        basis_checksum = written_basis$checksum,
        center = TRUE
      )
    )
  )

  result <- read_lna(archive)
  expected <- matrix(aperm(arr, c(4, 1, 2, 3)), nrow = 3, ncol = 4)
  expect_equal(result$stash$input, expected, tolerance = 1e-8)
})

test_that("embed.shared_basis rejects mask mismatches", {
  dir <- withr::local_tempdir()
  basis_mask <- array(c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE), dim = c(3, 2, 1))
  wrong_mask <- array(c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE), dim = c(3, 2, 1))
  written_basis <- write_shared_basis(
    list(
      loadings = diag(4), mask = basis_mask, kind = "mask-check",
      grid = shared_basis_test_grid(basis_mask), support_type = "masked"
    ),
    file.path(dir, "bases"),
    "mask-check"
  )

  archive_dir <- file.path(dir, "sub-01", "func")
  dir.create(archive_dir, recursive = TRUE)
  archive_path <- file.path(archive_dir, "sub-01_task-rest_bold.lna.h5")
  h5 <- hdf5r::H5File$new(archive_path, mode = "w")
  on.exit(if (h5$is_valid) h5$close_all(), add = TRUE)

  expect_error(
    forward_step.embed.shared_basis(
      "embed.shared_basis",
      list(
        type = "embed.shared_basis",
        params = list(
          basis_id = "mask-check",
          basis_checksum = written_basis$checksum,
          center = TRUE
        ),
        inputs = c("input"),
        outputs = c("coefficients")
      ),
      DataHandle$new(
        initial_stash = list(input = matrix(rnorm(8), nrow = 2, ncol = 4)),
        initial_meta = list(header = list(
          dims = dim(wrong_mask), transform = diag(4), spatial_units = "mm",
          voxel_ordering = "first-axis-fastest"
        )),
        plan = Plan$new(),
        h5 = h5,
        current_run_id = "run-01",
        mask_info = list(
          mask = wrong_mask,
          mask_array = wrong_mask,
          active_voxels = 4
        )
      )
    ),
    "support identity does not match"
  )
})

test_that("masked shared basis fails closed when input support is missing", {
  mask <- array(c(TRUE, FALSE, TRUE, TRUE), dim = c(2, 2, 1))
  basis <- shared_basis_test_identity(mask, support_type = "masked")
  error <- rlang::catch_cnd(
    .validate_shared_basis_mask(
      shared_basis_test_handle(mask = NULL, dims = dim(mask)),
      basis,
      "test:shared-basis"
    )
  )

  expect_s3_class(error, "lna_error_spatial_mismatch")
  expect_identical(error$field, "input_support")
  expect_identical(error$asset, "fixture-basis")
  expect_identical(error$expected_identity, basis$support_id)
})

test_that("shared basis without support identity cannot validate", {
  mask <- array(c(TRUE, FALSE, TRUE, TRUE), dim = c(2, 2, 1))
  basis <- shared_basis_test_identity(mask)
  basis$support_id <- ""
  error <- rlang::catch_cnd(
    .validate_shared_basis_mask(
      shared_basis_test_handle(mask), basis, "test:shared-basis"
    )
  )

  expect_s3_class(error, "lna_error_spatial_mismatch")
  expect_identical(error$field, "basis_support_id")
})

test_that("full-grid operation requires an explicit support descriptor", {
  full <- array(TRUE, dim = c(2, 2, 1))
  explicit <- shared_basis_test_identity(full, support_type = "full-grid")
  implicit <- explicit
  implicit$support_type <- "masked"
  no_mask <- shared_basis_test_handle(mask = NULL, dims = dim(full))

  expect_no_error(
    .validate_shared_basis_mask(no_mask, explicit, "test:shared-basis")
  )
  expect_error(
    .validate_shared_basis_mask(no_mask, implicit, "test:shared-basis"),
    "requires explicit input support",
    class = "lna_error_spatial_mismatch"
  )
  partial <- full
  partial[1] <- FALSE
  expect_error(
    .validate_shared_basis_mask(
      shared_basis_test_handle(partial), explicit, "test:shared-basis"
    ),
    "support identity does not match",
    class = "lna_error_spatial_mismatch"
  )
})

test_that("same dimensions with a different affine reject shared basis", {
  mask <- array(c(TRUE, FALSE, TRUE, TRUE), dim = c(2, 2, 1))
  basis <- shared_basis_test_identity(mask, affine = diag(4))
  shifted <- diag(4)
  shifted[1, 4] <- 10
  error <- rlang::catch_cnd(
    .validate_shared_basis_mask(
      shared_basis_test_handle(mask, affine = shifted),
      basis,
      "test:shared-basis"
    )
  )

  expect_s3_class(error, "lna_error_spatial_mismatch")
  expect_identical(error$field, "grid_id")
  expect_identical(error$expected_identity, basis$grid_id)
  expect_false(identical(error$actual_identity, basis$grid_id))
})

test_that("shared basis guard reports a missing input affine", {
  mask <- array(c(TRUE, FALSE, TRUE, TRUE), dim = c(2, 2, 1))
  basis <- shared_basis_test_identity(mask)
  handle <- DataHandle$new(
    initial_meta = list(header = list(dims = dim(mask))),
    mask_info = list(mask = mask)
  )
  error <- rlang::catch_cnd(
    .validate_shared_basis_mask(handle, basis, "test:shared-basis")
  )

  expect_s3_class(error, "lna_error_spatial_mismatch")
  expect_identical(error$field, "input_transform")
})
