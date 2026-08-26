# Tests for neuroim2 interoperability: as_neurovec, reader$series, as_latent_neurovec

# Helper: create a DenseNeuroVec for testing (only call inside tests that skip_if_not_installed)
make_test_vec <- function(spatial = c(4L, 4L, 2L), n_time = 5L,
                          spacing = c(2, 2, 3), transform = NULL) {
  arr <- array(rnorm(prod(spatial) * n_time), dim = c(spatial, n_time))
  sp <- if (is.null(transform)) {
    neuroim2::NeuroSpace(c(spatial, n_time), spacing = spacing)
  } else {
    neuroim2::NeuroSpace(c(spatial, n_time), trans = transform)
  }
  neuroim2::DenseNeuroVec(arr, sp)
}

test_affine <- function(kind = c("rotation", "shear")) {
  kind <- match.arg(kind)
  affine <- diag(4)
  affine[1:3, 1:3] <- if (kind == "rotation") {
    theta <- pi / 6
    matrix(
      c(2 * cos(theta), -2 * sin(theta), 0,
        2 * sin(theta),  2 * cos(theta), 0,
        0,               0,              3),
      nrow = 3,
      byrow = TRUE
    )
  } else {
    matrix(c(2, 0.35, 0, 0, 2, 0.2, 0, 0, 3), nrow = 3, byrow = TRUE)
  }
  affine[1:3, 4] <- c(-12, 8, 24)
  affine
}

delete_lna_header_attr <- function(file, name) {
  h5 <- hdf5r::H5File$new(file, mode = "r+")
  on.exit(h5$close_all(), add = TRUE)
  paths <- c("header/global", "header/runs/run-01")
  for (path in paths) {
    if (h5$exists(path)) {
      group <- h5[[path]]
      if (group$attr_exists(name)) group$attr_delete(name)
      group$close()
    }
  }
  invisible(file)
}

replace_lna_header_attr <- function(file, name, value) {
  h5 <- hdf5r::H5File$new(file, mode = "r+")
  on.exit(h5$close_all(), add = TRUE)
  paths <- c("header/global", "header/runs/run-01")
  for (path in paths) {
    if (h5$exists(path)) {
      group <- h5[[path]]
      if (group$attr_exists(name)) group$attr_delete(name)
      neuroarchive:::h5_attr_write(group, name, value)
      group$close()
    }
  }
  invisible(file)
}

test_that("as_neurovec round-trips DenseNeuroVec through LNA", {
  vec <- make_test_vec()
  tmp <- tempfile(fileext = ".lna.h5")
  on.exit(unlink(tmp), add = TRUE)
  write_lna(vec, tmp, transforms = "quant",
            transform_params = list(quant = list(bits = 16)))

  result <- as_neurovec(read_lna(tmp))

  expect_s4_class(result, "DenseNeuroVec")
  expect_equal(dim(result), c(4, 4, 2, 5))
})

test_that("as_neurovec preserves spatial metadata", {
  vec <- make_test_vec(spatial = c(3L, 3L, 3L), n_time = 10L, spacing = c(2, 2, 2))
  tmp <- tempfile(fileext = ".lna.h5")
  on.exit(unlink(tmp), add = TRUE)
  write_lna(vec, tmp, transforms = "quant",
            transform_params = list(quant = list(bits = 16)))

  result <- as_neurovec(read_lna(tmp))
  result_sp <- neuroim2::space(result)

  expect_equal(dim(result)[1:3], c(3, 3, 3))
  expect_equal(neuroim2::spacing(result_sp), c(2, 2, 2))
  expect_identical(
    attr(result, "lna_geometry_provenance")$source,
    "lna_header"
  )
})

test_that("rotated affine survives eager handle and file-path reconstruction", {
  affine <- test_affine("rotation")
  vec <- make_test_vec(transform = affine)
  tmp <- tempfile(fileext = ".lna.h5")
  on.exit(unlink(tmp), add = TRUE)
  write_lna(vec, tmp, transforms = "quant")

  from_handle <- as_neurovec(read_lna(tmp), file = tmp)
  from_path <- as_neurovec(tmp)

  expect_equal(neuroim2::trans(neuroim2::space(from_handle)), affine,
               tolerance = 1e-6)
  expect_equal(neuroim2::trans(neuroim2::space(from_path)), affine,
               tolerance = 1e-6)
  expect_equal(
    neuroim2::trans(neuroim2::space(from_handle)),
    neuroim2::trans(neuroim2::space(from_path))
  )
})

test_that("sheared affine survives lazy and reader-method reconstruction", {
  affine <- test_affine("shear")
  vec <- make_test_vec(transform = affine)
  tmp <- tempfile(fileext = ".lna.h5")
  on.exit(unlink(tmp), add = TRUE)
  write_lna(vec, tmp, transforms = "quant")

  reader <- read_lna(tmp, lazy = TRUE)
  on.exit(reader$close(), add = TRUE)
  from_s3 <- as_neurovec(reader)
  from_method <- reader$as_neurovec()

  expect_equal(neuroim2::trans(neuroim2::space(from_s3)), affine,
               tolerance = 1e-6)
  expect_equal(neuroim2::trans(neuroim2::space(from_method)), affine,
               tolerance = 1e-6)
})

test_that("as_neurovec from file path works", {
  vec <- make_test_vec(spatial = c(4L, 4L, 2L), n_time = 3L)
  tmp <- tempfile(fileext = ".lna.h5")
  on.exit(unlink(tmp), add = TRUE)
  write_lna(vec, tmp, transforms = "quant")

  result <- as_neurovec(tmp)
  expect_s4_class(result, "DenseNeuroVec")
  expect_equal(dim(result), c(4, 4, 2, 3))
})

test_that("as_neurovec.lna_reader works", {
  vec <- make_test_vec(spatial = c(4L, 4L, 2L), n_time = 3L)
  tmp <- tempfile(fileext = ".lna.h5")
  on.exit(unlink(tmp), add = TRUE)
  write_lna(vec, tmp, transforms = "quant")

  reader <- read_lna(tmp, lazy = TRUE)
  result <- as_neurovec(reader)
  reader$close()

  expect_s4_class(result, "DenseNeuroVec")
  expect_equal(dim(result), c(4, 4, 2, 3))
})

test_that("lna_reader$as_neurovec works", {
  vec <- make_test_vec(spatial = c(4L, 4L, 2L), n_time = 3L)
  tmp <- tempfile(fileext = ".lna.h5")
  on.exit(unlink(tmp), add = TRUE)
  write_lna(vec, tmp, transforms = "quant")

  reader <- read_lna(tmp, lazy = TRUE)
  result <- reader$as_neurovec()
  reader$close()

  expect_s4_class(result, "DenseNeuroVec")
  expect_equal(dim(result), c(4, 4, 2, 3))
})

test_that("lna_reader$series extracts time series by ROI mask", {
  vec <- make_test_vec(spatial = c(4L, 4L, 2L), n_time = 5L)
  tmp <- tempfile(fileext = ".lna.h5")
  on.exit(unlink(tmp), add = TRUE)
  write_lna(vec, tmp, transforms = "quant",
            transform_params = list(quant = list(bits = 16)))

  reader <- read_lna(tmp, lazy = TRUE)

  # ROI mask selecting a subset of voxels
  roi <- array(FALSE, dim = c(4, 4, 2))
  roi[1, 1, 1] <- TRUE
  roi[2, 2, 1] <- TRUE

  ts <- reader$series(roi)
  reader$close()

  expect_true(is.matrix(ts))
  # stash$input is time x voxels
  expect_equal(nrow(ts), 5)  # 5 time points
  expect_equal(ncol(ts), 2)  # 2 voxels in ROI
})

test_that("lna_reader$series with time subsetting", {
  vec <- make_test_vec(spatial = c(4L, 4L, 2L), n_time = 10L)
  tmp <- tempfile(fileext = ".lna.h5")
  on.exit(unlink(tmp), add = TRUE)
  write_lna(vec, tmp, transforms = "quant",
            transform_params = list(quant = list(bits = 16)))

  reader <- read_lna(tmp, lazy = TRUE)
  roi <- array(TRUE, dim = c(4, 4, 2))  # all voxels
  ts <- reader$series(roi, time_idx = 1:3)
  reader$close()

  expect_true(is.matrix(ts))
  expect_equal(nrow(ts), 3)  # 3 time points
  expect_equal(ncol(ts), 32) # all 32 voxels
})

test_that("neuroim2 series() works on as_neurovec result", {
  vec <- make_test_vec(spatial = c(4L, 4L, 2L), n_time = 5L)
  tmp <- tempfile(fileext = ".lna.h5")
  on.exit(unlink(tmp), add = TRUE)
  write_lna(vec, tmp, transforms = "quant",
            transform_params = list(quant = list(bits = 16)))

  result <- as_neurovec(tmp)

  # series() should work on the DenseNeuroVec
  ts <- neuroim2::series(result, 1, 1, 1)
  expect_length(ts, 5)
})

test_that("ROI subset returns SparseNeuroVec", {
  affine <- test_affine("shear")
  vec <- make_test_vec(spatial = c(4L, 4L, 2L), n_time = 5L,
                       transform = affine)
  tmp <- tempfile(fileext = ".lna.h5")
  on.exit(unlink(tmp), add = TRUE)
  write_lna(vec, tmp, transforms = "quant",
            transform_params = list(quant = list(bits = 16)))

  roi <- array(FALSE, dim = c(4, 4, 2))
  roi[1:2, 1:2, 1] <- TRUE  # 4 voxels

  result <- as_neurovec(read_lna(tmp, roi_mask = roi), file = tmp)

  expect_s4_class(result, "SparseNeuroVec")
  expect_equal(dim(result), c(4, 4, 2, 5))
  expect_equal(neuroim2::trans(neuroim2::space(result)), affine,
               tolerance = 1e-6)
  expect_equal(
    neuroim2::trans(neuroim2::space(methods::slot(result, "mask"))),
    affine,
    tolerance = 1e-6
  )

  # series() works on SparseNeuroVec
  ts <- neuroim2::series(result, 1, 1, 1)
  expect_length(ts, 5)
})

test_that("full read returns DenseNeuroVec", {
  vec <- make_test_vec(spatial = c(4L, 4L, 2L), n_time = 3L)
  tmp <- tempfile(fileext = ".lna.h5")
  on.exit(unlink(tmp), add = TRUE)
  write_lna(vec, tmp, transforms = "quant")

  result <- as_neurovec(read_lna(tmp))
  expect_s4_class(result, "DenseNeuroVec")
})

test_that("missing dimensions fail closed for fMRI reconstruction", {
  vec <- make_test_vec()
  tmp <- tempfile(fileext = ".lna.h5")
  on.exit(unlink(tmp), add = TRUE)
  write_lna(vec, tmp, transforms = "quant")
  delete_lna_header_attr(tmp, "dims")

  expect_error(
    as_neurovec(tmp),
    "spatial dimensions are missing",
    class = "lna_error_geometry"
  )
})

test_that("missing and malformed affines fail closed", {
  vec <- make_test_vec()
  missing <- tempfile(fileext = ".lna.h5")
  malformed <- tempfile(fileext = ".lna.h5")
  on.exit(unlink(c(missing, malformed)), add = TRUE)
  write_lna(vec, missing, transforms = "quant")
  write_lna(vec, malformed, transforms = "quant")
  delete_lna_header_attr(missing, "transform")
  replace_lna_header_attr(malformed, "transform", c(1, 0, 0, 1))

  expect_error(
    as_neurovec(missing),
    "affine is missing",
    class = "lna_error_geometry"
  )
  expect_error(
    as_neurovec(malformed),
    "finite numeric 4 x 4 matrix",
    class = "lna_error_geometry"
  )
})

test_that("identity geometry requires explicit opt-in and records provenance", {
  arr <- array(rnorm(3 * 4 * 2 * 5), dim = c(3, 4, 2, 5))
  tmp <- tempfile(fileext = ".lna.h5")
  on.exit(unlink(tmp), add = TRUE)
  write_lna(arr, tmp, transforms = "quant")

  expect_error(as_neurovec(tmp), class = "lna_error_geometry")

  result <- as_neurovec(tmp, assume_identity = TRUE)
  provenance <- attr(result, "lna_geometry_provenance")
  expect_equal(dim(result), dim(arr))
  expect_equal(neuroim2::trans(neuroim2::space(result)), diag(4))
  expect_identical(provenance$source, "explicit_identity_assumption")
  expect_true(provenance$assumed_identity)
  expect_identical(provenance$schema_version, "1.0.0")
})
