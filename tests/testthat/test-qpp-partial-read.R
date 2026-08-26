partial_qpp_file <- function(y, delta = 1, time_chunk = 4L,
                             parcel_chunk = 2L) {
  y <- as.matrix(y)
  handle <- DataHandle$new(
    initial_stash = list(parcel_series = y),
    initial_meta = list(parcel_mean_assets = list(
      signal_units = "percent-signal-change"
    )),
    plan = Plan$new(), run_ids = "run-01", current_run_id = "run-01"
  )
  desc <- list(
    type = "qpp",
    params = list(
      delta = delta, quantizer_policy = "fixed",
      units = "percent-signal-change",
      time_chunk = time_chunk, parcel_chunk = parcel_chunk
    ),
    inputs = list("parcel_series"), outputs = list("qpp_packed")
  )
  out <- run_transform_step("forward", "qpp", desc, handle, 0L)
  path <- tempfile(fileext = ".h5")
  h5 <- neuroarchive:::open_h5(path, "w")
  materialise_plan(h5, out$plan, run_ids = "run-01")
  neuroarchive:::close_h5_safely(h5)
  list(
    path = path,
    reconstructed = sweep(
      round(sweep(y, 2L, rep(delta, length.out = ncol(y)), "/")),
      2L, rep(delta, length.out = ncol(y)), "*"
    )
  )
}

decode_partial_qpp <- function(file, time_idx = NULL, parcel_idx = NULL,
                               cache = .qpp_new_tile_cache(0)) {
  h5 <- neuroarchive:::open_h5(file, "r")
  on.exit(neuroarchive:::close_h5_safely(h5), add = TRUE)
  group <- h5[["transforms"]]
  on.exit(safe_h5_close(group), add = TRUE)
  transforms <- discover_transforms(group)
  name <- transforms$name[transforms$type == "qpp"][[1L]]
  desc <- read_json_descriptor(group, name)
  handle <- DataHandle$new(
    h5 = h5,
    initial_meta = list(qpp_tile_cache = cache),
    subset = list(time_idx = time_idx, parcel_idx = parcel_idx),
    current_run_id = "run-01"
  )
  run_transform_step("invert", "qpp", desc, handle, 0L)
}

partial_parcel_qpp_file <- function() {
  mask <- array(TRUE, dim = c(8, 1, 1))
  x <- array(seq_len(8 * 12), dim = c(8, 1, 1, 12))
  labels <- array(rep(1:4, each = 2), dim = dim(mask))
  affine <- diag(c(2, 2, 2, 1))
  header <- list(
    dims = dim(mask), spacing = c(2, 2, 2), origin = c(0, 0, 0),
    transform = affine, spatial_units = "mm",
    voxel_ordering = "first-axis-fastest"
  )
  result <- core_write(
    x, transforms = c("parcel.mean", "qpp"), mask = mask, header = header,
    transform_params = list(
      "parcel.mean" = list(
        labels = labels, fit_scope = "fixed", source_id = "roi-atlas-v1",
        signal_units = "percent-signal-change", spatial_units = "mm",
        voxel_ordering = "first-axis-fastest",
        map_grid_id = .lna_grid_id(
          dim(mask), affine, "mm", "first-axis-fastest"
        )
      ),
      qpp = list(
        delta = 0.5, quantizer_policy = "fixed",
        units = "percent-signal-change", time_chunk = 4L,
        parcel_chunk = 1L
      )
    )
  )
  path <- tempfile(fileext = ".h5")
  h5 <- neuroarchive:::open_h5(path, "w")
  materialise_plan(
    h5, result$plan, run_ids = "run-01",
    header = header
  )
  neuroarchive:::close_h5_safely(h5)
  vox_time <- convert_to_masked_vox_time_matrix(x, mask)
  parcel <- vapply(1:4, function(k) {
    colMeans(vox_time[labels[mask] == k, , drop = FALSE])
  }, numeric(12))
  reconstructed <- sweep(round(parcel / 0.5), 2L, rep(0.5, 4), "*")
  list(path = path, mask = mask, x = x, labels = labels,
       parcel = reconstructed)
}

test_that("QPP reads only intersecting tiles, byte ranges, and parcel steps", {
  y <- outer(1:16, 1:8, function(t, k) t * 100 + k)
  fixture <- partial_qpp_file(y, delta = 0.25, time_chunk = 4L,
                              parcel_chunk = 2L)
  on.exit(unlink(fixture$path), add = TRUE)
  time_idx <- c(16L, 1L, 1L, 8L)
  parcel_idx <- c(8L, 1L, 8L)
  out <- decode_partial_qpp(fixture$path, time_idx, parcel_idx)
  expect_equal(
    out$stash$parcel_series,
    fixture$reconstructed[time_idx, parcel_idx, drop = FALSE]
  )
  io <- out$meta$qpp_io_receipt
  expect_identical(io$access_method, "H5D$read(args=...)")
  expect_identical(io$index_rows_read, 6L)
  expect_identical(io$index_rows_total, 16)
  expect_lt(io$payload_bytes_read, io$payload_total_bytes)
  expect_identical(io$payload_read_calls, 6L)
  expect_identical(io$step_values_read, 2L)
  expect_identical(io$step_values_total, 8)
  expect_identical(io$tile_digest_rows_read, 6L)
  expect_identical(io$step_digest_rows_read, 2L)
  expect_true(io$selected_tile_digests_verified)
  expect_true(io$selected_step_digests_verified)
  expect_true(io$physical_time_partial)
  expect_true(io$physical_parcel_partial)
  expect_false(io$full_step_checksum_verified)
  expect_true(all(vapply(
    io$read_calls, function(call) identical(call$access_method, "H5D$read(args=...)"),
    logical(1)
  )))
})

test_that("cross-tile, duplicate, unsorted, empty, and boundary selections work", {
  y <- outer(1:9, 1:5, function(t, k) t * 10 + k)
  fixture <- partial_qpp_file(y, time_chunk = 4L, parcel_chunk = 2L)
  on.exit(unlink(fixture$path), add = TRUE)

  time_idx <- c(9L, 4L, 5L, 1L, 5L)
  parcel_idx <- c(5L, 2L, 3L, 2L)
  out <- decode_partial_qpp(fixture$path, time_idx, parcel_idx)
  expect_equal(
    out$stash$parcel_series,
    fixture$reconstructed[time_idx, parcel_idx, drop = FALSE]
  )

  no_time <- decode_partial_qpp(fixture$path, integer(), c(2L, 1L))
  expect_identical(dim(no_time$stash$parcel_series), c(0L, 2L))
  expect_identical(no_time$meta$qpp_io_receipt$index_rows_read, 0L)
  expect_identical(no_time$meta$qpp_io_receipt$payload_bytes_read, 0)

  no_parcel <- decode_partial_qpp(fixture$path, c(1L, 9L), integer())
  expect_identical(dim(no_parcel$stash$parcel_series), c(2L, 0L))
  expect_identical(no_parcel$meta$qpp_io_receipt$index_rows_read, 0L)
  expect_identical(no_parcel$meta$qpp_io_receipt$step_values_read, 0L)

  expect_error(decode_partial_qpp(fixture$path, 0L, 1L), "time_idx")
  expect_error(decode_partial_qpp(fixture$path, 1L, 6L), "parcel_idx")
})

test_that("reader-owned decoded-tile cache hits and evicts within its byte bound", {
  y <- outer(1:12, 1:4, function(t, k) t + k)
  fixture <- partial_qpp_file(y, time_chunk = 4L, parcel_chunk = 2L)
  on.exit(unlink(fixture$path), add = TRUE)
  cache <- .qpp_new_tile_cache(max_bytes = 400)

  first <- decode_partial_qpp(fixture$path, 1L, 1L, cache)
  expect_identical(first$meta$qpp_io_receipt$tiles_decoded, 1L)
  second <- decode_partial_qpp(fixture$path, 2L, 2L, cache)
  expect_identical(second$meta$qpp_io_receipt$tile_cache_hits, 1L)
  expect_identical(second$meta$qpp_io_receipt$payload_bytes_read, 0)

  third <- decode_partial_qpp(fixture$path, 5L, 1L, cache)
  expect_lte(cache$bytes, cache$max_bytes)
  expect_lte(length(cache$order), 1L)
  expect_gte(cache$evictions, 1L)
  expect_identical(third$meta$qpp_io_receipt$cache_max_bytes, 400)
})

test_that("tile cache identity changes with the stored-byte digest", {
  first_fixture <- partial_qpp_file(
    outer(1:8, 1:4, function(t, k) t + k),
    delta = 0.25, time_chunk = 4L, parcel_chunk = 2L
  )
  second_fixture <- partial_qpp_file(
    outer(1:8, 1:4, function(t, k) 100 + 3 * t + k),
    delta = 0.25, time_chunk = 4L, parcel_chunk = 2L
  )
  on.exit(unlink(c(first_fixture$path, second_fixture$path)), add = TRUE)
  cache <- .qpp_new_tile_cache(max_bytes = 4096)

  first <- decode_partial_qpp(first_fixture$path, 1:2, 1:2, cache)
  expect_identical(first$meta$qpp_io_receipt$tiles_decoded, 1L)
  second <- decode_partial_qpp(second_fixture$path, 1:2, 1:2, cache)
  expect_identical(second$meta$qpp_io_receipt$tile_cache_hits, 0L)
  expect_identical(second$meta$qpp_io_receipt$tiles_decoded, 1L)
  expect_false(identical(
    first$meta$qpp_decoded_q_checksum,
    second$meta$qpp_decoded_q_checksum
  ))
  expect_false(identical(
    first$stash$parcel_series,
    second$stash$parcel_series
  ))
})

test_that("ROI gather reads only requested label positions and parcel tiles", {
  fixture <- partial_parcel_qpp_file()
  on.exit(unlink(fixture$path), add = TRUE)
  roi <- array(FALSE, dim = dim(fixture$mask))
  roi[c(1, 8)] <- TRUE
  time_idx <- c(12L, 1L, 5L, 5L)
  reader <- read_lna(fixture$path, lazy = TRUE, validate = FALSE)
  on.exit(reader$close(), add = TRUE)
  out <- reader$data(roi_mask = roi, time_idx = time_idx)
  expect_equal(
    out$stash$input,
    fixture$parcel[time_idx, c(1L, 4L), drop = FALSE]
  )
  io <- out$meta$qpp_io_receipt
  expect_identical(io$label_values_read, 2L)
  expect_identical(io$label_values_total, 8)
  expect_identical(io$step_values_read, 2L)
  expect_identical(io$step_values_total, 4)
  expect_identical(io$index_rows_read, 6L)
  expect_lt(io$payload_bytes_read, io$payload_total_bytes)
  expect_identical(out$meta$qpp_selected_parcels, c(1L, 4L))
  expect_true(out$subset_applied$roi_mask)
  expect_true(out$subset_applied$time_idx)

  empty_roi <- array(FALSE, dim = dim(fixture$mask))
  empty <- reader$data(roi_mask = empty_roi, time_idx = c(1L, 12L))
  expect_identical(dim(empty$stash$input), c(2L, 0L))
  expect_identical(empty$meta$qpp_io_receipt$index_rows_read, 0L)
  expect_identical(empty$meta$qpp_io_receipt$payload_bytes_read, 0)
  expect_identical(empty$meta$qpp_io_receipt$step_values_read, 0L)
  expect_identical(empty$meta$qpp_io_receipt$label_values_read, 0L)
})

test_that("bytes read and result memory scale with touched QPP tiles", {
  y <- outer(1:32, 1:16, function(t, k) sin(t / 3) + k)
  fixture <- partial_qpp_file(y, delta = 0.01, time_chunk = 4L,
                              parcel_chunk = 2L)
  on.exit(unlink(fixture$path), add = TRUE)
  short <- decode_partial_qpp(fixture$path, 1:2, 1:2)
  full <- decode_partial_qpp(fixture$path)
  short_io <- short$meta$qpp_io_receipt
  full_io <- full$meta$qpp_io_receipt
  expect_lt(short_io$payload_bytes_read, full_io$payload_bytes_read)
  expect_identical(full_io$payload_bytes_read, full_io$payload_total_bytes)
  expect_lt(
    as.numeric(object.size(short$stash$parcel_series)),
    as.numeric(object.size(full$stash$parcel_series))
  )
})

test_that("LNANeuroVec routes joint time and ROI indexing through QPP ranges", {
  skip_if_not_installed("neuroim2")
  fixture <- partial_parcel_qpp_file()
  on.exit(unlink(fixture$path), add = TRUE)
  lnv <- LNANeuroVec(fixture$path)
  on.exit(lna_close(lnv), add = TRUE)
  value <- lnv[1, 1, 1, c(12, 1), drop = FALSE]
  expect_equal(
    as.numeric(value),
    fixture$parcel[c(12L, 1L), 1L],
    tolerance = 1e-12
  )
  io <- lnv@.state$reader$data_cache$meta$qpp_io_receipt
  expect_true(io$physical_time_partial)
  expect_true(io$physical_parcel_partial)
  expect_false(lna_loaded(lnv))
  expect_true(lnv@.state$pipeline_caps$physical_time_idx)
  expect_true(lnv@.state$pipeline_caps$physical_roi_mask)
  expect_identical(lnv@.state$pipeline_caps$roi_strategy, "parcel-gather")
})
