# Installed-package smoke test for the experimental analysis-safe QPP profile.
# Run only after installing a built neuroarchive tarball.

library(neuroarchive)

mask <- array(TRUE, dim = c(8, 1, 1))
labels <- array(rep(1:4, each = 2), dim = dim(mask))
affine <- diag(c(2, 2, 2, 1))
x <- array(
  outer(seq_len(8), seq_len(12), function(v, t) v + sin(t / 3)),
  dim = c(8, 1, 1, 12)
)
normalization <- list(
  method = "none", fit_scope = "fixed", source_id = "consumer-values-v1",
  calibration_identity = NULL, fit_id = NULL
)
profile <- lna_profile_analysis_safe(
  labels = labels, delta = 0.01,
  signal_units = "percent-signal-change", source_id = "consumer-atlas-v1",
  normalization_provenance = normalization,
  map_grid_id = lna_grid_identity(dim(mask), affine),
  time_chunk = 4L, parcel_chunk = 2L
)
stopifnot(
  identical(profile$lifecycle, "experimental"),
  isTRUE(profile$available),
  identical(profile$allowed_split_unit, "pointwise"),
  grepl("individual TR", profile$allowed_split_description, fixed = TRUE)
)

file <- tempfile(fileext = ".lna.h5")
on.exit(unlink(file), add = TRUE)
write_lna(
  x, file, transforms = profile$transforms,
  transform_params = profile$transform_params, mask = mask,
  header = list(
    dims = c(8L, 1L, 1L), spacing = c(2, 2, 2), origin = c(0, 0, 0),
    transform = affine
  ),
  run_id = "run-01", checksum = "none"
)

validation <- validate_lna(
  file, checksum = FALSE, qpp = "all", receipt = TRUE
)
stopifnot(
  inherits(validation, "lna_validation_receipt"),
  isTRUE(validation$qpp[[1L]]$q_roundtrip_verified)
)

contract <- analysis_contract(file)
stopifnot(identical(contract$allowed_split_unit, "pointwise"))
folds <- data.frame(
  observation_id = seq_len(12), fold = rep(c("train", "test"), 6)
)
split_receipt <- assert_split_compatible(file, folds)
stopifnot(isTRUE(split_receipt$compatible))

roi <- array(FALSE, dim = dim(mask))
roi[c(1, 8)] <- TRUE
subset <- read_lna(
  file, roi_mask = roi, time_idx = c(11L, 2L, 5L), validate = FALSE
)
stopifnot(
  identical(dim(subset$stash$input), c(3L, 2L)),
  isTRUE(subset$meta$qpp_io_receipt$physical_time_partial),
  isTRUE(subset$meta$qpp_io_receipt$physical_parcel_partial),
  isTRUE(subset$meta$qpp_io_receipt$selected_tile_digests_verified),
  isTRUE(subset$meta$qpp_io_receipt$selected_step_digests_verified)
)

wavelet <- lna_profile_multiresolution_run()
stopifnot(
  inherits(wavelet, "lna_profile_unavailable"),
  !isTRUE(wavelet$available), is.null(wavelet$transforms),
  identical(wavelet$allowed_split_unit, "run")
)

message("QPP installed-consumer check: PASS")
