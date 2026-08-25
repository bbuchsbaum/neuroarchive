parcel_mean_fixture <- function(weights = NULL, labels = NULL, header = list()) {
  mask <- array(c(TRUE, FALSE, TRUE, TRUE), dim = c(2, 2, 1))
  if (length(header) == 0L) {
    header <- list(
      dims = dim(mask),
      transform = diag(4),
      spatial_units = "mm",
      voxel_ordering = "first-axis-fastest"
    )
  }
  x <- array(seq_len(12), dim = c(2, 2, 1, 3))
  if (is.null(labels)) {
    labels <- array(c(10, 0, 30, 10), dim = dim(mask))
  }
  if (is.null(weights)) {
    weights <- array(c(1, 0, 2, 3), dim = dim(mask))
  }
  handle <- DataHandle$new(
    initial_stash = list(input = x),
    initial_meta = list(mask = mask, header = header),
    plan = Plan$new(),
    run_ids = "run-01",
    current_run_id = "run-01",
    mask_info = list(mask = mask, mask_array = mask, active_voxels = sum(mask))
  )
  params <- list(
    labels = labels,
    weights = weights,
    fit_scope = "fixed",
    source_id = "fixture-atlas-v1",
    signal_units = "percent-signal-change",
    spatial_units = "mm",
    voxel_ordering = "first-axis-fastest"
  )
  if (!is.null(header$dims) && !is.null(header$transform)) {
    params$map_grid_id <- .lna_grid_id(
      header$dims[1:3],
      header$transform,
      params$spatial_units,
      params$voxel_ordering
    )
  }
  list(handle = handle, params = params, mask = mask, x = x)
}

run_parcel_mean_fixture <- function(fixture, output = "parcel_series") {
  desc <- list(
    type = "parcel.mean",
    params = fixture$params,
    inputs = list("input"),
    outputs = list(output)
  )
  run_transform_step("forward", "parcel.mean", desc, fixture$handle, 0L)
}

test_that("parcel.mean matches an independent weighted matrix oracle", {
  fixture <- parcel_mean_fixture()
  out <- run_parcel_mean_fixture(fixture)
  expected <- cbind(c(13 / 4, 29 / 4, 45 / 4), c(3, 7, 11))
  colnames(expected) <- c("10", "30")
  expect_equal(out$stash$parcel_series, expected)

  desc <- out$plan$descriptors[[1L]]
  expect_identical(desc$params$matrix_order, "time-by-parcel")
  expect_identical(desc$params$label_dtype, "uint16")
  expect_identical(desc$params$weight_policy, "fixed-voxel-weights")
  expect_identical(desc$params$fit_scope, "fixed")
  expect_identical(desc$version, "1.1.0")
  expect_match(desc$params$asset_id, "^sha256:[0-9a-f]{64}$")
  expect_match(desc$params$asset_path, "^/spatial/qpp_assets/[0-9a-f]{64}$")
  expect_match(desc$params$map_checksum, "^sha256:")
  expect_identical(vapply(desc$datasets, `[[`, character(1), "role"),
                   c("parcel_labels", "source_parcel_ids", "parcel_counts",
                     "parcel_weights"))
  expect_false(any(vapply(out$plan$datasets$role, identical, logical(1),
                          "parcel_series")))
  expect_false(any(vapply(out$plan$payloads, function(x) {
    is.matrix(x) && identical(dim(x), dim(expected))
  }, logical(1))))
})

test_that("parcel.mean rejects noncanonical calibration identities before fitting", {
  fixture <- parcel_mean_fixture()
  fixture$params$fit_scope <- "external-calibration"
  fixture$params$calibration_id <- "external-atlas-v1"
  fixture$params$calibration_checksum <- "not-a-sha256"
  error <- rlang::catch_cnd(run_parcel_mean_fixture(fixture))
  expect_s3_class(error, "lna_error_checksum")
  expect_identical(error$code, "noncanonical_sha256")
  expect_identical(error$location, "forward_step.parcel.mean")
  expect_length(fixture$handle$plan$descriptors, 0L)
  expect_length(fixture$handle$plan$payloads, 0L)
})

test_that("parcel.mean unweighted means and gather reconstruction are exact", {
  fixture <- parcel_mean_fixture(weights = NULL)
  fixture$params$weights <- NULL
  out <- run_parcel_mean_fixture(fixture)
  expected <- cbind(c(2.5, 6.5, 10.5), c(3, 7, 11))
  colnames(expected) <- c("10", "30")
  expect_equal(out$stash$parcel_series, expected)
  desc <- out$plan$descriptors[[1L]]
  expect_identical(desc$params$weight_policy, "uniform")

  restored <- run_transform_step("invert", "parcel.mean", desc, out, 0L)
  oracle <- expected[, c(1, 2, 1), drop = FALSE]
  dimnames(oracle) <- NULL
  expect_equal(restored$stash$input, oracle)
  expect_identical(restored$meta$parcel_reconstruction$matrix_order,
                   "time-by-active-voxel")
  expect_identical(restored$mask_info$mask, fixture$mask)
})

test_that("parcel.mean inverse applies time and full-volume ROI selections", {
  fixture <- parcel_mean_fixture()
  out <- run_parcel_mean_fixture(fixture)
  desc <- out$plan$descriptors[[1L]]
  out$subset <- list(
    time_idx = c(1L, 3L),
    roi_mask = array(c(TRUE, FALSE, FALSE, TRUE), dim = dim(fixture$mask))
  )
  restored <- run_transform_step("invert", "parcel.mean", desc, out, 0L)
  expect_equal(unname(restored$stash$input),
               matrix(c(13 / 4, 45 / 4, 13 / 4, 45 / 4), nrow = 2))
  expect_true(restored$subset_applied$time_idx)
  expect_true(restored$subset_applied$roi_mask)
})

test_that("noncontiguous and singleton source labels remap deterministically", {
  fixture <- parcel_mean_fixture(
    labels = array(c(9000, 0, 2, 47), dim = c(2, 2, 1))
  )
  fixture$params$weights <- NULL
  out <- run_parcel_mean_fixture(fixture)
  assets <- out$meta$parcel_mean_assets
  expect_identical(assets$source_ids, c(2, 47, 9000))
  expect_identical(assets$dense_labels, c(3L, 1L, 2L))
  expect_identical(assets$counts, c(1L, 1L, 1L))
  expect_equal(unname(out$stash$parcel_series),
               cbind(c(3, 7, 11), c(4, 8, 12), c(1, 5, 9)))
})

test_that("parcel.mean rejects missing, foreign, empty, and zero-weight maps", {
  missing <- parcel_mean_fixture(
    labels = array(c(10, 0, NA, 10), dim = c(2, 2, 1))
  )
  expect_error(run_parcel_mean_fixture(missing), "missing values")

  foreign <- parcel_mean_fixture(
    labels = array(c(10, 99, 30, 10), dim = c(2, 2, 1))
  )
  expect_error(run_parcel_mean_fixture(foreign), "zero outside")

  empty <- parcel_mean_fixture()
  empty$params$source_labels <- c(10, 20, 30)
  expect_error(run_parcel_mean_fixture(empty), "empty parcels")

  zero <- parcel_mean_fixture(
    weights = array(c(0, 0, 2, 0), dim = c(2, 2, 1))
  )
  expect_error(run_parcel_mean_fixture(zero), "positive total weight")
})

test_that("active-order maps fail closed on reordered or foreign masks", {
  fixture <- parcel_mean_fixture()
  fixture$params$labels <- c(10, 30, 10)
  expect_error(run_parcel_mean_fixture(fixture), "map_mask_checksum")
  fixture$params$map_mask_checksum <- .compute_mask_checksum(fixture$mask)
  expect_no_error(run_parcel_mean_fixture(fixture))
  fixture$params$map_mask_checksum <- paste0("sha256:", strrep("0", 64))
  expect_error(run_parcel_mean_fixture(fixture), "map_mask_checksum")
})

test_that("grid, affine, and fitted-map provenance cannot be omitted", {
  transform <- diag(4)
  header <- list(dims = c(2L, 2L, 1L), transform = transform)
  fixture <- parcel_mean_fixture(header = header)
  fixture$params$map_grid_id <- NULL
  expect_error(run_parcel_mean_fixture(fixture), "grid identity")
  fixture$params$map_grid_id <- .lna_grid_id(
    header$dims, transform, "mm", "first-axis-fastest"
  )
  fixture$params$map_affine_checksum <- .lna_affine_id(transform)
  expect_no_error(run_parcel_mean_fixture(fixture))

  missing_affine <- parcel_mean_fixture(
    header = list(dims = c(2L, 2L, 1L))
  )
  missing_affine$params$map_grid_id <- fixture$params$map_grid_id
  affine_error <- rlang::catch_cnd(run_parcel_mean_fixture(missing_affine))
  expect_s3_class(affine_error, "lna_error_transform_step")
  expect_match(conditionMessage(affine_error), "affine is missing")
  affine_cause <- affine_error$parent_error
  expect_s3_class(affine_cause, "lna_error_spatial_mismatch")
  expect_identical(affine_cause$field, "input_transform")
  expect_identical(affine_cause$asset, "fixture-atlas-v1")
  expect_identical(
    affine_cause$expected_identity,
    fixture$params$map_grid_id
  )
  expect_true(is.na(affine_cause$actual_identity))

  shifted <- diag(4)
  shifted[1, 4] <- 5
  foreign_grid <- parcel_mean_fixture(header = list(
    dims = c(2L, 2L, 1L), transform = shifted
  ))
  foreign_grid$params$map_grid_id <- fixture$params$map_grid_id
  grid_error <- rlang::catch_cnd(run_parcel_mean_fixture(foreign_grid))
  expect_s3_class(grid_error, "lna_error_transform_step")
  grid_cause <- grid_error$parent_error
  expect_s3_class(grid_cause, "lna_error_spatial_mismatch")
  expect_identical(grid_cause$field, "grid_id")
  expect_identical(
    grid_cause$expected_identity,
    fixture$params$map_grid_id
  )
  expect_false(identical(
    grid_cause$actual_identity,
    fixture$params$map_grid_id
  ))

  fitted <- parcel_mean_fixture()
  fitted$params$fit_scope <- "subject"
  expect_error(run_parcel_mean_fixture(fitted), "fit_provenance")
  fitted$params$fit_provenance <- list(subject = "sub-calibration")
  expect_no_error(run_parcel_mean_fixture(fitted))

  external <- parcel_mean_fixture()
  external$params$fit_scope <- "external-calibration"
  expect_error(run_parcel_mean_fixture(external), "calibration_id")
})

test_that("parcel.mean uses uint32 dense labels when K exceeds uint16", {
  k <- 65536L
  mask <- array(rep(TRUE, k), dim = c(k, 1, 1))
  fixture <- parcel_mean_fixture()
  fixture$mask <- mask
  fixture$x <- array(seq_len(k), dim = c(k, 1, 1, 1))
  large_header <- list(dims = dim(mask), transform = diag(4))
  fixture$handle <- DataHandle$new(
    initial_stash = list(input = fixture$x),
    initial_meta = list(mask = mask, header = large_header),
    plan = Plan$new(), run_ids = "run-01", current_run_id = "run-01",
    mask_info = list(mask = mask, mask_array = mask, active_voxels = k)
  )
  fixture$params$labels <- array(seq_len(k), dim = dim(mask))
  fixture$params$weights <- NULL
  fixture$params$map_grid_id <- .lna_grid_id(
    dim(mask), diag(4), "mm", "first-axis-fastest"
  )
  out <- run_parcel_mean_fixture(fixture)
  expect_identical(out$plan$descriptors[[1L]]$params$label_dtype, "uint32")
})

test_that("parcel.mean rejects accidental multi-run list input", {
  fixture <- parcel_mean_fixture()
  fixture$handle$stash$input <- list(run1 = fixture$x, run2 = fixture$x)
  expect_error(run_parcel_mean_fixture(fixture), "one run at a time")
})

test_that("core_write accepts a sparse mask for parcel.mean", {
  fixture <- parcel_mean_fixture()
  result <- core_write(
    fixture$x,
    transforms = "parcel.mean",
    transform_params = list("parcel.mean" = fixture$params),
    mask = fixture$mask,
    header = fixture$handle$meta$header
  )
  expect_true(result$handle$has_key("parcel.mean_s0_out"))
  expect_identical(result$handle$mask_info$active_voxels, 3L)
})

test_that("parcel.mean assets materialize and support HDF5 gather inversion", {
  fixture <- parcel_mean_fixture()
  fixture$params$weights <- NULL
  out <- run_parcel_mean_fixture(fixture)
  expected_series <- out$stash$parcel_series
  tmp <- local_tempfile(fileext = ".h5")
  h5 <- neuroarchive:::open_h5(tmp, mode = "w")
  materialise_plan(
    h5, out$plan, run_ids = "run-01",
    header = list(dims = c(2L, 2L, 1L), spacing = c(2, 2, 3),
                  origin = c(0, 0, 0), transform = diag(c(2, 2, 3, 1)))
  )
  desc <- read_json_descriptor(h5[["transforms"]], "00_parcel.mean.json")
  asset_path <- sub("^/", "", desc$params$asset_path)
  expect_true(h5$exists(paste0(asset_path, "/parcel_labels")))
  expect_true(h5$exists(paste0(asset_path, "/source_parcel_ids")))
  expect_true(h5$exists(paste0(asset_path, "/parcel_counts")))
  expect_false(h5$exists(paste0(asset_path, "/parcel_weights")))

  read_handle <- DataHandle$new(
    initial_stash = list(parcel_series = expected_series),
    initial_meta = list(), h5 = h5, current_run_id = "run-01",
    mask_info = list(mask = fixture$mask, active_voxels = sum(fixture$mask))
  )
  restored <- run_transform_step("invert", "parcel.mean", desc, read_handle, 0L)
  expect_equal(restored$stash$input,
               expected_series[, c(1, 2, 1), drop = FALSE],
               ignore_attr = TRUE)
  vec <- as_neurovec(restored)
  expect_s4_class(vec, "SparseNeuroVec")
  expect_equal(dim(vec), c(2, 2, 1, 3))
  expect_equal(neuroim2::spacing(neuroim2::space(vec)), c(2, 2, 3))
  neuroarchive:::close_h5_safely(h5)
  expect_true(validate_lna(tmp))
})

test_that("parcel.mean schema is valid JSON", {
  schema <- test_inst_file("schemas", "parcel.mean.schema.json")
  expect_true(jsonlite::validate(paste(readLines(schema, warn = FALSE),
                                       collapse = "\n")))
})
