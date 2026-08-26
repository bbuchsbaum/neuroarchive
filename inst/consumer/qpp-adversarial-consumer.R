# Installed-package adversarial court for the experimental analysis-safe QPP
# profile. This script uses neuroarchive's public API plus direct HDF5 mutation;
# it does not load test helpers or reach into the package namespace.

library(neuroarchive)

na_files <- character()
na_track <- function(path) {
  na_files <<- unique(c(na_files, path))
  path
}
on.exit(unlink(na_files), add = TRUE)

na_condition <- function(fun, class, code, location) {
  condition <- tryCatch({
    fun()
    NULL
  }, error = identity)
  if (is.null(condition)) stop("expected an error, but the operation succeeded")
  if (!inherits(condition, class)) {
    stop(
      "expected class ", class, "; got ",
      paste(base::class(condition), collapse = "/")
    )
  }
  if (!identical(condition$code, code)) {
    stop("expected code ", code, "; got ", condition$code %||% "<missing>")
  }
  if (!identical(condition$location, location)) {
    stop("expected location ", location, "; got ", condition$location %||% "<missing>")
  }
  invisible(condition)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

na_read_descriptor <- function(group, name) {
  dataset <- group[[name]]
  on.exit(try(dataset$close(), silent = TRUE), add = TRUE)
  jsonlite::fromJSON(dataset$read(), simplifyVector = FALSE)
}

na_write_descriptor <- function(group, name, descriptor) {
  json <- as.character(jsonlite::toJSON(
    descriptor, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA
  ))
  if (group$exists(name)) group$link_delete(name)
  dtype <- hdf5r::h5types$H5T_C_S1$copy()
  space <- hdf5r::H5S$new("scalar")
  dataset <- NULL
  on.exit({
    if (!is.null(dataset)) try(dataset$close(), silent = TRUE)
    try(space$close(), silent = TRUE)
    try(dtype$close(), silent = TRUE)
  }, add = TRUE)
  dtype$set_size(Inf)
  dtype$set_cset(hdf5r::h5const$H5T_CSET_UTF8)
  dataset <- group$create_dataset(
    name, dtype = dtype, space = space, chunk_dims = NULL
  )
  dataset[] <- json
  invisible(NULL)
}

na_descriptors <- function(group) {
  names <- group$ls()$name
  out <- lapply(names, function(name) na_read_descriptor(group, name))
  names(out) <- names
  out
}

na_role_path <- function(descriptor, role) {
  matches <- Filter(function(x) identical(x$role, role), descriptor$datasets)
  if (length(matches) != 1L) stop("expected one dataset role: ", role)
  matches[[1L]]$path
}

na_copy_mutate <- function(source, mutate) {
  target <- na_track(tempfile(fileext = ".lna.h5"))
  stopifnot(file.copy(source, target, overwrite = TRUE))
  h5 <- hdf5r::H5File$new(target, mode = "r+")
  tryCatch({
    group <- h5[["transforms"]]
    mutate(h5, group, na_descriptors(group))
    h5$flush()
  }, finally = {
    try(h5$close_all(), silent = TRUE)
  })
  target
}

na_mutate_descriptor <- function(source, type = "qpp", index = 1L, mutate) {
  na_copy_mutate(source, function(h5, group, descriptors) {
    names <- names(descriptors)[vapply(
      descriptors, function(x) identical(x$type, type), logical(1)
    )]
    name <- names[[index]]
    na_write_descriptor(group, name, mutate(descriptors[[name]]))
  })
}

# Minimal independent Pure-R QPP1 decoder used only as an installed oracle.
na_ref_unpack <- function(bytes, count, width) {
  count <- as.integer(count)
  width <- as.integer(width)
  if (length(bytes) != ceiling(count * width / 8)) stop("oracle byte length")
  if (width == 0L) return(numeric(count))
  out <- numeric(count)
  bit_cursor <- 0L
  bytes <- as.integer(bytes)
  for (i in seq_len(count)) {
    for (bit in 0:(width - 1L)) {
      byte <- bit_cursor %/% 8L + 1L
      offset <- bit_cursor %% 8L
      out[[i]] <- out[[i]] +
        ((bytes[[byte]] %/% (2^offset)) %% 2) * (2^bit)
      bit_cursor <- bit_cursor + 1L
    }
  }
  out
}

na_ref_inverse <- function(anchors, residuals, predictor) {
  if (identical(predictor, "P0")) return(residuals)
  if (identical(predictor, "P1")) {
    return(c(anchors, anchors + cumsum(residuals)))
  }
  out <- c(anchors, numeric(length(residuals)))
  for (i in seq_along(residuals)) {
    j <- i + 2L
    out[[j]] <- residuals[[i]] + 2 * out[[j - 1L]] - out[[j - 2L]]
  }
  out
}

na_ref_decode_tile <- function(bytes) {
  bytes <- as.raw(bytes)
  pos <- 1L
  take <- function(n) {
    if (pos + n - 1L > length(bytes)) stop("truncated oracle tile")
    value <- if (n == 0L) raw() else bytes[pos:(pos + n - 1L)]
    pos <<- pos + n
    value
  }
  u8 <- function() as.integer(take(1L))
  u16 <- function() {
    value <- as.integer(take(2L))
    value[[1L]] + 256 * value[[2L]]
  }
  u32 <- function() {
    value <- as.integer(take(4L))
    sum(value * c(1, 2^8, 2^16, 2^24))
  }
  i32 <- function() {
    value <- u32()
    if (value >= 2^31) value - 2^32 else value
  }
  if (!identical(take(4L), charToRaw("QPP1"))) stop("oracle magic")
  n_time <- u16()
  n_parcel <- u16()
  q <- matrix(0, nrow = n_time, ncol = n_parcel)
  for (parcel in seq_len(n_parcel)) {
    predictor_id <- u8()
    width <- u8()
    n_anchor <- u8()
    reserved <- u8()
    n_residual <- u16()
    payload_bytes <- u32()
    predictor <- c("P0", "P1", "P2")[[predictor_id + 1L]]
    expected_anchor <- c(P0 = 0L, P1 = 1L, P2 = 2L)[[predictor]]
    if (reserved != 0L || width > 34L || n_anchor != expected_anchor ||
        n_residual != n_time - n_anchor ||
        payload_bytes != ceiling(n_residual * width / 8)) {
      stop("oracle header")
    }
    anchors <- if (n_anchor) {
      vapply(seq_len(n_anchor), function(i) i32(), numeric(1))
    } else {
      numeric()
    }
    unsigned <- na_ref_unpack(take(payload_bytes), n_residual, width)
    residuals <- ifelse(unsigned %% 2 == 0, unsigned / 2, -(unsigned + 1) / 2)
    q[, parcel] <- na_ref_inverse(anchors, residuals, predictor)
  }
  if (pos != length(bytes) + 1L) stop("oracle trailing bytes")
  q
}

na_oracle_q <- function(file, run_id) {
  h5 <- hdf5r::H5File$new(file, mode = "r")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
  group <- h5[["transforms"]]
  descriptors <- na_descriptors(group)
  qpp <- Filter(function(x) {
    identical(x$type, "qpp") && identical(x$run_id, run_id)
  }, descriptors)[[1L]]
  payload_dataset <- h5[[na_role_path(qpp, "qpp_payload")]]
  index_dataset <- h5[[na_role_path(qpp, "qpp_tile_index")]]
  payload <- as.raw(payload_dataset$read())
  index <- index_dataset$read()
  payload_dataset$close()
  index_dataset$close()
  q <- matrix(NA_real_, qpp$params$time_count, qpp$params$parcel_count)
  for (row in seq_len(nrow(index))) {
    stopifnot(as.numeric(index$outer_codec[[row]]) == 0)
    start <- as.numeric(index$byte_offset[[row]]) + 1L
    count <- as.numeric(index$byte_length[[row]])
    tile <- na_ref_decode_tile(payload[seq.int(start, length.out = count)])
    time <- seq.int(
      as.numeric(index$time_start[[row]]) + 1L,
      length.out = as.numeric(index$time_count[[row]])
    )
    parcel <- seq.int(
      as.numeric(index$parcel_start[[row]]) + 1L,
      length.out = as.numeric(index$parcel_count[[row]])
    )
    q[time, parcel] <- tile
  }
  stopifnot(!anyNA(q))
  q
}

mask <- array(TRUE, dim = c(8, 1, 1))
labels <- array(rep(1:4, each = 2), dim = dim(mask))
affine <- diag(c(2, 2, 2, 1))
geometry <- list(
  dims = dim(mask), spacing = c(2, 2, 2), origin = c(0, 0, 0),
  transform = affine, spatial_units = "mm",
  voxel_ordering = "first-axis-fastest"
)
base <- outer(seq_len(8), seq_len(12), function(v, t) v + sin(t / 3))
runs <- list(
  "run-01" = array(base, dim = c(8, 1, 1, 12)),
  "run-02" = array(base + 10 + cos(rep(seq_len(12), each = 8) / 4),
                    dim = c(8, 1, 1, 12))
)
delta <- c(0.01, 0.02, 0.03, 0.04)
normalization <- list(
  method = "none", fit_scope = "fixed", source_id = "court-values-v1",
  calibration_identity = NULL, fit_id = NULL
)
profile <- lna_profile_analysis_safe(
  labels = labels, delta = delta,
  signal_units = "percent-signal-change", source_id = "court-atlas-v1",
  normalization_provenance = normalization,
  map_grid_id = lna_grid_identity(dim(mask), affine),
  time_chunk = 4L, parcel_chunk = 2L
)
archive <- na_track(tempfile(fileext = ".lna.h5"))
write_lna(
  runs, archive, transforms = profile$transforms,
  transform_params = profile$transform_params, mask = mask,
  run_geometry = list("run-01" = geometry, "run-02" = geometry),
  checksum = "none"
)

validation <- validate_lna(
  archive, checksum = FALSE, qpp = "all", receipt = TRUE
)
stopifnot(
  length(validation$qpp) == 2L,
  all(vapply(validation$qpp, function(x) isTRUE(x$q_roundtrip_verified), logical(1))),
  identical(analysis_contract(archive)$allowed_split_unit, "pointwise")
)
folds <- data.frame(observation_id = seq_len(12), fold = rep(c("A", "B"), 6))
stopifnot(isTRUE(assert_split_compatible(archive, folds)$compatible))

eager <- read_lna(
  archive, run_id = "run-*", validate = FALSE, output_dtype = "float64"
)
lazy <- read_lna(
  archive, run_id = "run-*", validate = FALSE, output_dtype = "float64",
  lazy = TRUE
)
stopifnot(identical(names(eager), names(runs)), identical(names(lazy), names(runs)))
roi <- array(FALSE, dim = dim(mask))
roi[c(1L, 8L)] <- TRUE
for (run_id in names(runs)) {
  single <- na_track(tempfile(fileext = ".lna.h5"))
  write_lna(
    runs[[run_id]], single, transforms = profile$transforms,
    transform_params = profile$transform_params, mask = mask,
    header = geometry, run_id = run_id, checksum = "none"
  )
  independent <- read_lna(
    single, run_id = run_id, validate = FALSE, output_dtype = "float64"
  )
  stopifnot(isTRUE(all.equal(
    eager[[run_id]]$stash$input, independent$stash$input, tolerance = 1e-12
  )))
  range <- lazy[[run_id]]$data(roi_mask = roi, time_idx = 2:4)
  independent_range <- read_lna(
    single, run_id = run_id, validate = FALSE, output_dtype = "float64",
    roi_mask = roi, time_idx = 2:4
  )
  io <- range$meta$qpp_io_receipt
  stopifnot(
    isTRUE(all.equal(
      range$stash$input, independent_range$stash$input, tolerance = 1e-12
    )),
    identical(range$meta$qpp_selected_parcels, c(1L, 4L)),
    io$label_values_read == 2L,
    io$tile_digest_rows_read == io$index_rows_read,
    io$step_digest_rows_read == io$step_values_read,
    isTRUE(io$selected_tile_digests_verified),
    isTRUE(io$selected_step_digests_verified),
    io$payload_bytes_read < io$payload_total_bytes,
    length(range$stash$input) == 6L,
    as.numeric(object.size(range$stash$input)) < 4096,
    io$cache_bytes <= io$cache_max_bytes
  )
}
lazy[[1L]]$close()
stopifnot(is.null(lazy[[1L]]$h5), isTRUE(lazy[[2L]]$h5$is_valid))
lazy[[2L]]$close()

vox_time <- matrix(as.numeric(runs[["run-01"]]), nrow = 8L)
parcel_series <- vapply(seq_len(4L), function(parcel) {
  colMeans(vox_time[as.numeric(labels) == parcel, , drop = FALSE])
}, numeric(12L))
expected_q <- round(sweep(parcel_series, 2L, delta, "/"))
stopifnot(identical(na_oracle_q(archive, "run-01"), expected_q))

# Public-constructor checksum rejection.
na_condition(
  function() lna_profile_analysis_safe(
    labels = labels, delta = delta, signal_units = "psc", source_id = "atlas",
    parcellation_fit_scope = "external-calibration",
    parcellation_calibration_id = "bad-calibration",
    parcellation_calibration_checksum = "not-a-sha256",
    normalization_provenance = normalization,
    map_grid_id = lna_grid_identity(dim(mask), affine)
  ),
  "lna_error_checksum", "noncanonical_sha256",
  "lna_profile_analysis_safe:parcellation_calibration_checksum"
)

# Standalone parcel-series QPP cannot be relabelled analysis-safe.
standalone <- na_track(tempfile(fileext = ".lna.h5"))
write_lna(
  array(seq_len(32), dim = c(4, 1, 1, 8)), standalone,
  transforms = "qpp",
  transform_params = list(qpp = list(
    delta = 1, quantizer_policy = "fixed", units = "psc",
    normalization_provenance = normalization,
    time_chunk = 4L, parcel_chunk = 2L
  )),
  checksum = "none"
)
qpp_only_claim <- na_mutate_descriptor(standalone, mutate = function(desc) {
  desc$profile_id <- "analysis-safe"
  desc$params$profile_id <- "analysis-safe"
  desc
})
na_condition(
  function() validate_lna(qpp_only_claim, checksum = FALSE),
  "lna_error_qpp_binding", "qpp_parcel_binding",
  "analysis_contract:qpp_parcel_binding"
)

missing <- na_mutate_descriptor(archive, mutate = function(desc) {
  desc$params$parcel_binding <- NULL
  desc
})
na_condition(
  function() validate_lna(missing, checksum = FALSE),
  "lna_error_qpp_binding", "qpp_parcel_binding",
  "analysis_contract:qpp_parcel_binding"
)

duplicate <- na_copy_mutate(archive, function(h5, group, descriptors) {
  qpp <- Filter(function(x) identical(x$type, "qpp"), descriptors)[[1L]]
  parcel <- descriptors[[qpp$params$parcel_binding$descriptor]]
  na_write_descriptor(group, "04_parcel.mean.json", parcel)
})
na_condition(
  function() validate_lna(duplicate, checksum = FALSE),
  "lna_error_qpp_binding", "qpp_parcel_duplicate",
  "validate_lna:qpp:parcel_binding"
)

wrong_run <- na_mutate_descriptor(archive, mutate = function(desc) {
  desc$params$parcel_binding$run_id <- "run-02"
  desc
})
na_condition(
  function() validate_lna(wrong_run, checksum = FALSE),
  "lna_error_qpp_binding", "qpp_parcel_wrong_run",
  "validate_lna:qpp:parcel_binding"
)

wrong_map <- na_mutate_descriptor(archive, mutate = function(desc) {
  desc$params$parcel_binding$map_checksum <- paste0("sha256:", strrep("0", 64))
  desc
})
na_condition(
  function() validate_lna(wrong_map, checksum = FALSE),
  "lna_error_qpp_binding", "qpp_parcel_missing",
  "validate_lna:qpp:parcel_binding"
)

wrong_support <- na_mutate_descriptor(archive, mutate = function(desc) {
  desc$params$parcel_binding$mask_checksum <- paste0("sha256:", strrep("1", 64))
  desc
})
na_condition(
  function() validate_lna(wrong_support, checksum = FALSE),
  "lna_error_qpp_binding", "qpp_parcel_missing",
  "validate_lna:qpp:parcel_binding"
)

payload_flip <- na_copy_mutate(archive, function(h5, group, descriptors) {
  desc <- Filter(function(x) {
    identical(x$type, "qpp") && identical(x$run_id, "run-01")
  }, descriptors)[[1L]]
  dataset <- h5[[na_role_path(desc, "qpp_payload")]]
  bytes <- as.raw(dataset$read())
  bytes[[1L]] <- as.raw(bitwXor(as.integer(bytes[[1L]]), 1L))
  dataset$write(args = list(seq_along(bytes)), value = as.integer(bytes))
  dataset$close()
})
na_condition(
  function() read_lna(payload_flip, run_id = "run-01", time_idx = 1L,
                      validate = FALSE),
  "lna_error_qpp_integrity", "qpp_tile_digest_mismatch",
  "invert_step.qpp:tile_digest"
)
na_condition(
  function() validate_lna(payload_flip, checksum = FALSE, qpp = "all"),
  "lna_error_qpp_integrity", "qpp_tile_digest_mismatch",
  "validate_lna:qpp:tile_digest"
)

step_change <- na_copy_mutate(archive, function(h5, group, descriptors) {
  desc <- Filter(function(x) {
    identical(x$type, "qpp") && identical(x$run_id, "run-01")
  }, descriptors)[[1L]]
  dataset <- h5[[na_role_path(desc, "qpp_parcel_steps")]]
  steps <- as.numeric(dataset$read())
  steps[[1L]] <- steps[[1L]] * 2
  dataset$write(args = list(seq_along(steps)), value = steps)
  dataset$close()
})
na_condition(
  function() read_lna(step_change, run_id = "run-01", time_idx = 1L,
                      validate = FALSE),
  "lna_error_qpp_integrity", "qpp_step_digest_mismatch",
  "invert_step.qpp:step_digest"
)

index_change <- na_copy_mutate(archive, function(h5, group, descriptors) {
  desc <- Filter(function(x) {
    identical(x$type, "qpp") && identical(x$run_id, "run-01")
  }, descriptors)[[1L]]
  dataset <- h5[[na_role_path(desc, "qpp_tile_index")]]
  index <- dataset$read()
  index$digest_row[[1L]] <- index$digest_row[[1L]] + 1
  dataset$write(args = list(seq_len(nrow(index))), value = index)
  dataset$close()
})
na_condition(
  function() read_lna(index_change, run_id = "run-01", time_idx = 1L,
                      validate = FALSE),
  "lna_error_qpp_integrity", "qpp_digest_index",
  "invert_step.qpp:index_digest"
)

unknown_version <- na_mutate_descriptor(archive, mutate = function(desc) {
  desc$version <- "9.9.9"
  desc
})
na_condition(
  function() validate_lna(unknown_version, checksum = FALSE),
  "lna_error_required_transform", "required_transform",
  "validate_lna:qpp:decoder_version"
)

message("QPP installed adversarial court: PASS")
