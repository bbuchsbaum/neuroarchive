# Tests for LNANeuroVec lazy subset-dispatch NeuroVec

skip_if_not_installed("neuroim2")

make_lna_file <- function(spatial = c(4L, 4L, 2L), n_time = 5L,
                          spacing = c(2, 2, 3), transforms = "quant",
                          transform_params = list(quant = list(bits = 16))) {
  arr <- array(rnorm(prod(spatial) * n_time), dim = c(spatial, n_time))
  sp <- neuroim2::NeuroSpace(c(spatial, n_time), spacing = spacing)
  vec <- neuroim2::DenseNeuroVec(arr, sp)
  tmp <- tempfile(fileext = ".lna.h5")
  write_lna(vec, tmp, transforms = transforms,
            transform_params = transform_params)
  list(file = tmp, vec = vec, arr = arr)
}

# --- Construction ---

test_that("LNANeuroVec constructs from file path", {
  info <- make_lna_file()
  on.exit(unlink(info$file), add = TRUE)

  lnv <- LNANeuroVec(info$file)
  expect_s4_class(lnv, "LNANeuroVec")
  expect_s4_class(lnv, "SparseNeuroVec")
  expect_equal(dim(lnv), c(4, 4, 2, 5))
  lna_close(lnv)
})

test_that("LNANeuroVec constructs from lna_reader", {
  info <- make_lna_file()
  on.exit(unlink(info$file), add = TRUE)

  reader <- read_lna(info$file, lazy = TRUE)
  lnv <- LNANeuroVec(reader)
  expect_s4_class(lnv, "LNANeuroVec")
  expect_equal(dim(lnv), c(4, 4, 2, 5))
  lna_close(lnv)
})

test_that("as_lna_neurovec works from file path", {
  info <- make_lna_file()
  on.exit(unlink(info$file), add = TRUE)

  lnv <- as_lna_neurovec(info$file)
  expect_s4_class(lnv, "LNANeuroVec")
  lna_close(lnv)
})

# --- Metadata (no data access) ---

test_that("dim and space are correct without materializing", {
  info <- make_lna_file(spatial = c(3L, 3L, 3L), n_time = 10L, spacing = c(1, 1, 1))
  on.exit(unlink(info$file), add = TRUE)

  lnv <- LNANeuroVec(info$file)
  expect_equal(dim(lnv), c(3, 3, 3, 10))
  expect_false(lna_loaded(lnv))  # no full materialization

  sp <- neuroim2::space(lnv)
  expect_s4_class(sp, "NeuroSpace")
  lna_close(lnv)
})

# --- Lazy data access ---

test_that("series() extracts time series lazily", {
  info <- make_lna_file()
  on.exit(unlink(info$file), add = TRUE)

  lnv <- LNANeuroVec(info$file)
  ts <- neuroim2::series(lnv, 1L)
  expect_length(ts, 5)
  # quant-only pipeline supports ROI subsetting — should NOT fully materialize
  expect_false(lna_loaded(lnv))
  lna_close(lnv)
})

test_that("series with 3D coords works", {
  info <- make_lna_file()
  on.exit(unlink(info$file), add = TRUE)

  lnv <- LNANeuroVec(info$file)
  ts <- neuroim2::series(lnv, 1, 1, 1)
  expect_length(ts, 5)
  lna_close(lnv)
})

test_that("[i,j,k,t] indexing returns the correct value from lazy state", {
  info <- make_lna_file()
  on.exit(unlink(info$file), add = TRUE)

  lnv <- LNANeuroVec(info$file)
  val <- lnv[1, 1, 1, 1]
  expect_true(is.numeric(val))
  expect_length(val, 1)
  expect_equal(as.numeric(val), info$arr[1, 1, 1, 1], tolerance = 0.05)
  lna_close(lnv)
})

test_that("[[t]] extracts a single volume lazily", {
  info <- make_lna_file()
  on.exit(unlink(info$file), add = TRUE)

  lnv <- LNANeuroVec(info$file)
  vol <- lnv[[1]]
  expect_s4_class(vol, "SparseNeuroVol")
  expect_equal(dim(vol), c(4, 4, 2))
  expect_false(lna_loaded(lnv))  # single-volume read, not full
  lna_close(lnv)
})

test_that("sub_vector extracts time subset lazily", {
  info <- make_lna_file()
  on.exit(unlink(info$file), add = TRUE)

  lnv <- LNANeuroVec(info$file)
  sv <- neuroim2::sub_vector(lnv, 1:3)
  expect_s4_class(sv, "SparseNeuroVec")
  expect_equal(dim(sv), c(4, 4, 2, 3))
  expect_false(lna_loaded(lnv))
  lna_close(lnv)
})

# --- Caching ---

test_that("repeated access to same voxels hits cache", {
  info <- make_lna_file()
  on.exit(unlink(info$file), add = TRUE)

  lnv <- LNANeuroVec(info$file)
  ts1 <- neuroim2::series(lnv, 1L)
  ts2 <- neuroim2::series(lnv, 1L)
  expect_equal(ts1, ts2)
  # LRU cache should be populated (not full mat for quant-only)
  expect_true(!is.null(lnv@.state$last_roi))
  lna_close(lnv)
})

# --- Lifecycle ---

test_that("lna_close makes further access error", {
  info <- make_lna_file()
  on.exit(unlink(info$file), add = TRUE)

  lnv <- LNANeuroVec(info$file)
  lna_close(lnv)
  expect_error(neuroim2::series(lnv, 1L), "closed")
})

test_that("show method works", {
  info <- make_lna_file()
  on.exit(unlink(info$file), add = TRUE)

  lnv <- LNANeuroVec(info$file)
  output <- capture.output(show(lnv))
  expect_true(any(grepl("LNANeuroVec", output)))
  expect_true(any(grepl("lazy", output)))
  lna_close(lnv)
})

# --- Round-trip accuracy ---

test_that("round-trip preserves values within tolerance", {
  info <- make_lna_file()
  on.exit(unlink(info$file), add = TRUE)

  lnv <- LNANeuroVec(info$file)
  orig_ts <- info$arr[1, 1, 1, ]
  lna_ts <- neuroim2::series(lnv, 1, 1, 1)
  # 16-bit quant has very small error; use generous tolerance
  expect_equal(as.numeric(lna_ts), orig_ts, tolerance = 0.05)
  lna_close(lnv)
})

# --- Multi-transform ---

test_that("works with basis + quant pipeline", {
  info <- make_lna_file(
    spatial = c(4L, 4L, 2L), n_time = 10L,
    transforms = c("basis", "quant"),
    transform_params = list(
      basis = list(k = 3, allow_dense_passthrough = TRUE),
      quant = list(bits = 12)
    )
  )
  on.exit(unlink(info$file), add = TRUE)

  lnv <- LNANeuroVec(info$file)
  expect_equal(dim(lnv)[1:3], c(4, 4, 2))

  # Time-only subsetting works.
  vol <- lnv[[1]]
  expect_s4_class(vol, "SparseNeuroVol")

  sv <- neuroim2::sub_vector(lnv, 1:3)
  expect_equal(dim(sv)[4], 3)

  # Spatial access should stay subsetted and avoid full materialization.
  ts <- expect_no_warning(neuroim2::series(lnv, 1L))
  expect_true(is.numeric(ts))
  expect_equal(length(ts), 10)
  expect_false(lna_loaded(lnv))
  expect_equal(as.numeric(ts), info$arr[1, 1, 1, ], tolerance = 0.05)
  lna_close(lnv)
})

test_that("works with basis + embed + quant pipeline without full materialization", {
  info <- make_lna_file(
    spatial = c(4L, 4L, 2L), n_time = 10L,
    transforms = c("basis", "embed", "quant"),
    transform_params = list(
      basis = list(k = 3),
      embed = list(
        basis_path = "/basis/00_basis/matrix",
        center_data_with = "/basis/00_basis/center"
      ),
      quant = list(bits = 12)
    )
  )
  on.exit(unlink(info$file), add = TRUE)

  ref <- as_dense_mat(read_lna(info$file)$stash$input)

  lnv <- LNANeuroVec(info$file)
  ts <- expect_no_warning(neuroim2::series(lnv, 1L))
  expect_equal(as.numeric(ts), ref[, 1], tolerance = 1e-6)
  expect_false(lna_loaded(lnv))
  lna_close(lnv)
})
