# Public profile constructors for QPP lifecycle-safe workflows.

.profile_qpp_string <- function(value, field) {
  if (!is.character(value) || length(value) != 1L || !nzchar(value)) {
    abort_lna(paste0(field, " must be one non-empty string"),
              .subclass = "lna_error_validation",
              location = paste0("lna_profile_analysis_safe:", field))
  }
  value
}

.profile_qpp_calibration <- function(scope, id, checksum, prefix) {
  if (!identical(scope, "external-calibration")) {
    return(list(id = NULL, checksum = NULL))
  }
  id <- .profile_qpp_string(id, paste0(prefix, "_calibration_id"))
  checksum <- .profile_qpp_string(
    checksum, paste0(prefix, "_calibration_checksum")
  )
  checksum <- .lna_canonical_sha256(
    checksum, paste0(prefix, "_calibration_checksum"), canonicalize = TRUE,
    location = paste0(
      "lna_profile_analysis_safe:", prefix, "_calibration_checksum"
    )
  )
  list(id = id, checksum = checksum)
}

#' Construct the analysis-safe QPP profile
#'
#' @description
#' Creates write arguments for a pointwise-quantized parcel-series archive.
#' The returned object is directly consumable by [write_lna()] through its
#' `transforms` and `transform_params` fields. It does not fit a parcellation,
#' normalization, or quantizer from the data being encoded.
#'
#' @param labels Full-volume parcel-label array, using zero outside the mask, or
#'   an active-mask-order label vector accompanied by `map_mask_checksum`.
#' @param delta One fixed positive quantizer step or one per represented parcel.
#' @param signal_units Scientific amplitude units, for example
#'   `"percent-signal-change"`.
#' @param source_id Immutable parcellation/atlas identifier.
#' @param spatial_units Spatial coordinate units.
#' @param voxel_ordering Voxel linearization convention used by `map_grid_id`.
#' @param weights Optional fixed voxel weights.
#' @param parcellation_fit_scope Either `"fixed"` or
#'   `"external-calibration"`.
#' @param parcellation_calibration_id,parcellation_calibration_checksum Required
#'   immutable provenance when the parcellation is externally calibrated.
#' @param quantizer_policy Either `"fixed"` or `"external-calibration"`.
#' @param quantizer_calibration_id,quantizer_calibration_checksum Required
#'   immutable provenance for externally calibrated steps.
#' @param normalization_provenance Required fixed or disjoint-external
#'   normalization record. Use `method = "none"` only when the supplied values
#'   genuinely require no fitted normalization.
#' @param map_mask_checksum Required when labels or weights are supplied only in
#'   active-mask order.
#' @param map_grid_id Required canonical identity tying the parcel map to grid
#'   dimensions, affine, spatial units, and voxel ordering.
#' @param map_affine_checksum Optional checksum tying the parcel map to an affine.
#' @param time_chunk,parcel_chunk Independently decodable tile dimensions.
#' @return An `lna_profile` list. Its lifecycle is currently `"experimental"`
#'   because the committed untouched scientific court missed one connectivity
#'   ceiling; this label does not change the pointwise leakage contract.
#'   `allowed_split_unit` uses the machine vocabulary returned by
#'   [analysis_contract()] (`"pointwise"` here), while
#'   `allowed_split_description` gives the corresponding reader-facing rule.
#' @section Estimand and error contract:
#' Input is a four-dimensional fMRI image. The codec reconstructs one fixed
#' parcel mean for every original TR, then gathers those means back to active
#' voxels. Quantization is pointwise with empirical error reported by
#' [validate_lna()]; integer prediction is exactly reversible and does not alter
#' the quantized samples.
#' @section Leakage contract:
#' Arbitrary later time-wise splits are compatible only when parcel labels,
#' weights, steps, and normalization are fixed or fitted on disjoint external
#' data. `fmrilatent` assets can satisfy this contract; package identity is not
#' a safety criterion.
#' @export
lna_profile_analysis_safe <- function(
    labels, delta, signal_units, source_id, spatial_units = "mm",
    voxel_ordering = "first-axis-fastest",
    weights = NULL,
    parcellation_fit_scope = c("fixed", "external-calibration"),
    parcellation_calibration_id = NULL,
    parcellation_calibration_checksum = NULL,
    quantizer_policy = c("fixed", "external-calibration"),
    quantizer_calibration_id = NULL,
    quantizer_calibration_checksum = NULL,
    normalization_provenance = NULL,
    map_mask_checksum = NULL, map_grid_id = NULL, map_affine_checksum = NULL,
    time_chunk = 128L, parcel_chunk = 128L) {
  parcellation_fit_scope <- match.arg(parcellation_fit_scope)
  quantizer_policy <- match.arg(quantizer_policy)
  if (!is.numeric(labels) || !length(labels) || any(!is.finite(labels)) ||
      any(labels != trunc(labels)) || any(labels < 0) || any(labels > 2^53 - 1)) {
    abort_lna("labels must be finite non-negative exact integers",
              .subclass = "lna_error_validation",
              location = "lna_profile_analysis_safe:labels")
  }
  if (!is.numeric(delta) || !length(delta) || any(!is.finite(delta)) ||
      any(delta <= 0)) {
    abort_lna("delta must contain finite positive steps",
              .subclass = "lna_error_validation",
              location = "lna_profile_analysis_safe:delta")
  }
  if (!is.null(weights) &&
      (!is.numeric(weights) || any(!is.finite(weights)) || any(weights < 0))) {
    abort_lna("weights must be finite and non-negative",
              .subclass = "lna_error_validation",
              location = "lna_profile_analysis_safe:weights")
  }
  signal_units <- .profile_qpp_string(signal_units, "signal_units")
  source_id <- .profile_qpp_string(source_id, "source_id")
  spatial_units <- .profile_qpp_string(spatial_units, "spatial_units")
  voxel_ordering <- .profile_qpp_string(voxel_ordering, "voxel_ordering")
  time_chunk <- .qpp_chunk_size(time_chunk, "time_chunk")
  parcel_chunk <- .qpp_chunk_size(parcel_chunk, "parcel_chunk")
  parcel_calibration <- .profile_qpp_calibration(
    parcellation_fit_scope, parcellation_calibration_id,
    parcellation_calibration_checksum, "parcellation"
  )
  quantizer_calibration <- .profile_qpp_calibration(
    quantizer_policy, quantizer_calibration_id,
    quantizer_calibration_checksum, "quantizer"
  )
  if (is.null(normalization_provenance)) {
    abort_lna(
      "normalization_provenance is required; analysis-safe profiles fail closed",
      .subclass = "lna_error_analysis_contract",
      location = "lna_profile_analysis_safe:normalization_provenance"
    )
  }
  normalization_provenance <- .analysis_normalization_record(
    normalization_provenance, canonicalize_checksums = TRUE
  )
  if (!normalization_provenance$fit_scope %in% c(
    "fixed", "external-calibration"
  )) {
    abort_lna(
      paste(
        "analysis-safe normalization fit_scope must be fixed or",
        "external-calibration"
      ),
      .subclass = "lna_error_analysis_contract",
      location = "lna_profile_analysis_safe:normalization_provenance"
    )
  }
  map_grid_id <- .lna_canonical_sha256(
    map_grid_id, "map_grid_id", canonicalize = TRUE,
    location = "lna_profile_analysis_safe:map_grid_id",
    subclass = "lna_error_spatial_identity"
  )
  map_mask_checksum <- .lna_canonical_sha256(
    map_mask_checksum, "map_mask_checksum", allow_null = TRUE,
    canonicalize = TRUE,
    location = "lna_profile_analysis_safe:map_mask_checksum",
    subclass = "lna_error_spatial_identity"
  )
  map_affine_checksum <- .lna_canonical_sha256(
    map_affine_checksum, "map_affine_checksum", allow_null = TRUE,
    canonicalize = TRUE,
    location = "lna_profile_analysis_safe:map_affine_checksum",
    subclass = "lna_error_spatial_identity"
  )
  parcel_params <- list(
    labels = labels, weights = weights, fit_scope = parcellation_fit_scope,
    source_id = source_id, signal_units = signal_units,
    spatial_units = spatial_units, voxel_ordering = voxel_ordering,
    calibration_id = parcel_calibration$id,
    calibration_checksum = parcel_calibration$checksum,
    normalization_provenance = normalization_provenance,
    map_mask_checksum = map_mask_checksum,
    map_grid_id = map_grid_id,
    map_affine_checksum = map_affine_checksum
  )
  qpp_params <- list(
    delta = as.numeric(delta), quantizer_policy = quantizer_policy,
    units = signal_units,
    calibration_id = quantizer_calibration$id,
    calibration_checksum = quantizer_calibration$checksum,
    normalization_provenance = normalization_provenance,
    allow_clip = FALSE, time_chunk = time_chunk, parcel_chunk = parcel_chunk
  )
  structure(list(
    profile_id = "analysis-safe", lifecycle = "experimental",
    available = TRUE,
    lifecycle_reason = paste(
      "implementation gates pass, but the committed untouched scientific court",
      "missed the maximum connectivity-edge ceiling"
    ),
    reconstructed_estimand = "fixed parcel mean at every original TR",
    error_contract = "pointwise quantization; empirical max error and RMSE in validation receipt",
    leakage_contract = "arbitrary TR splits only with fixed or disjoint-external fitted assets",
    allowed_split_unit = "pointwise",
    allowed_split_description = "individual TR observations may be assigned independently",
    transforms = c("parcel.mean", "qpp"),
    transform_params = list("parcel.mean" = parcel_params, qpp = qpp_params),
    gate_receipt = "benchmarks/qpp-gates-v1.json"
  ), class = c("lna_profile", "list"))
}

#' Describe the multiresolution-run profile decision
#'
#' @description
#' Returns the frozen no-go status for the partition-aware integer lifting
#' prototype. It is an executable status constructor, not write arguments: the
#' calibration-selected wavelet missed both the untouched scientific gate and
#' the predeclared rate margin, so no production transform was admitted.
#'
#' @return An unavailable `lna_profile` status object with links to the installed
#'   reference specification and benchmark receipt. Its machine split unit is
#'   `"run"`; the description states that complete runs must remain intact.
#' @section Leakage contract:
#' Any future implementation has complete-run (or larger declared unit) lossy
#' support. Within-run train/test splits are invalid; folds must be unions of
#' complete transformed units.
#' @export
lna_profile_multiresolution_run <- function() {
  structure(list(
    profile_id = "multiresolution-run", lifecycle = "prototype-no-go",
    available = FALSE,
    lifecycle_reason = paste(
      "calibration-selected haar-gentle used 95.317% of QPP payload bytes",
      "(required <=90%) and failed untouched GLM-t agreement"
    ),
    reconstructed_estimand = "parcel series at requested progressive detail bands",
    error_contract = "run-supported wavelet loss; no production bound admitted",
    leakage_contract = "folds must be unions of complete transformed units",
    allowed_split_unit = "run",
    allowed_split_description = "each transformed run must remain in one fold",
    prototype_reference = "prototypes/temporal-lift53-reference.R",
    specification = "spec/temporal-lift53-reference-v1.md",
    gate_receipt = "benchmarks/qpp-wavelet-prototype-v1.json",
    transforms = NULL, transform_params = NULL
  ), class = c("lna_profile_unavailable", "lna_profile", "list"))
}

#' Print an LNA profile decision
#'
#' @param x An object returned by an LNA profile constructor.
#' @param ... Reserved for future display options.
#' @return `x`, invisibly.
#' @export
print.lna_profile <- function(x, ...) {
  cat("LNA profile:", x$profile_id, "\n")
  cat("  lifecycle:", x$lifecycle, "\n")
  cat("  available for writing:", if (isTRUE(x$available)) "yes" else "no", "\n")
  cat("  estimand:", x$reconstructed_estimand, "\n")
  cat("  allowed split unit (machine):", x$allowed_split_unit, "\n")
  cat("  split interpretation:", x$allowed_split_description, "\n")
  cat("  reason:", x$lifecycle_reason, "\n")
  invisible(x)
}
