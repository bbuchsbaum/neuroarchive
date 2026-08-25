# Transform subsetting capability registry
.transform_subset_capabilities <- list(
  quant = list(roi_mask = TRUE, time_idx = TRUE),
  basis = list(roi_mask = TRUE, time_idx = TRUE),
  embed = list(roi_mask = TRUE, time_idx = TRUE),
  "embed.shared_basis" = list(roi_mask = TRUE, time_idx = TRUE),
  "parcel.mean" = list(
    roi_mask = TRUE, time_idx = TRUE, parcel_idx = TRUE,
    semantic = list(roi_mask = TRUE, time_idx = TRUE, parcel_idx = TRUE),
    physical = list(roi_mask = TRUE, time_idx = TRUE, parcel_idx = TRUE),
    temporal_granularity = "pass-through",
    parcel_granularity = "selected-labels",
    roi_strategy = "parcel-gather"
  ),
  qpp = list(
    roi_mask = FALSE, time_idx = TRUE, parcel_idx = TRUE,
    semantic = list(roi_mask = FALSE, time_idx = TRUE, parcel_idx = TRUE),
    physical = list(roi_mask = FALSE, time_idx = TRUE, parcel_idx = TRUE),
    temporal_granularity = "qpp-time-tile",
    parcel_granularity = "qpp-parcel-tile",
    roi_strategy = "parcel-index"
  ),
  delta = list(roi_mask = TRUE, time_idx = TRUE),
  temporal = list(roi_mask = TRUE, time_idx = TRUE),
  "spat.hrbf" = list(roi_mask = TRUE, time_idx = TRUE),
  "spat.hrbf_project" = list(roi_mask = TRUE, time_idx = TRUE),
  "spat.haar_octwave" = list(roi_mask = TRUE, time_idx = TRUE),
  "embed.transfer_hrbf_basis" = list(roi_mask = TRUE, time_idx = TRUE),
  sparsepca = list(roi_mask = TRUE, time_idx = TRUE),
  aggregate_runs = list(roi_mask = FALSE, time_idx = FALSE),
  "basis.empirical_hrbf_compressed" = list(roi_mask = FALSE, time_idx = FALSE)
)

.transform_lifecycle <- list(
  stable = c("quant", "delta", "basis", "embed", "embed.shared_basis", "temporal"),
  experimental = c(
    "parcel.mean",
    "qpp",
    "spat.hrbf",
    "spat.hrbf_project",
    "spat.haar_octwave",
    "embed.transfer_hrbf_basis",
    "sparsepca",
    "basis.empirical_hrbf_compressed"
  ),
  developer = c("aggregate_runs")
)

.transform_supports_subset <- function(type) {
  caps <- .transform_subset_capabilities[[type]] %||%
    list(roi_mask = FALSE, time_idx = FALSE)
  semantic <- caps$semantic %||% list(
    roi_mask = isTRUE(caps$roi_mask),
    time_idx = isTRUE(caps$time_idx),
    parcel_idx = isTRUE(caps$parcel_idx)
  )
  physical <- caps$physical %||% list(
    roi_mask = FALSE, time_idx = FALSE, parcel_idx = FALSE
  )
  caps$semantic <- semantic
  caps$physical <- physical
  caps$parcel_idx <- isTRUE(semantic$parcel_idx)
  caps$temporal_granularity <- caps$temporal_granularity %||% "full-dataset"
  caps$parcel_granularity <- caps$parcel_granularity %||% "full-dataset"
  caps$roi_strategy <- caps$roi_strategy %||% "direct"
  caps
}

.transform_lifecycle_status <- function(type, metadata = list()) {
  if (is.character(metadata$lifecycle) && length(metadata$lifecycle) == 1L) {
    return(metadata$lifecycle)
  }
  if (type %in% .transform_lifecycle$stable) {
    return("stable")
  }
  if (type %in% .transform_lifecycle$experimental) {
    return("experimental")
  }
  if (type %in% .transform_lifecycle$developer) {
    return("developer")
  }
  "custom"
}

.pipeline_roi_strategy <- function(types) {
  if (length(types) == 0) {
    return("direct")
  }

  if (all(c("parcel.mean", "qpp") %in% types)) {
    return("parcel-gather")
  }

  basis_idx <- which(types == "basis")
  if (length(basis_idx) == 0) {
    if (any(types == "embed.shared_basis")) {
      return("shared-coefficient")
    }
    if (any(startsWith(types, "embed"))) {
      return("coefficient")
    }
    return("direct")
  }

  segment_modes <- vapply(basis_idx, function(idx) {
    if (.basis_step_has_following_embed(types, idx)) "coefficient" else "dense-pass-through"
  }, character(1))

  if (length(unique(segment_modes)) == 1L) {
    segment_modes[[1]]
  } else {
    "mixed"
  }
}

.pipeline_supports_subset <- function(types) {
  caps <- lapply(types, .transform_supports_subset)
  roi_strategy <- .pipeline_roi_strategy(types)
  parcel_gather <- all(c("parcel.mean", "qpp") %in% types)
  temporal_granularities <- unique(vapply(caps, `[[`, character(1),
                                           "temporal_granularity"))
  temporal_granularities <- setdiff(temporal_granularities, "pass-through")
  parcel_granularities <- unique(vapply(caps, `[[`, character(1),
                                         "parcel_granularity"))
  parcel_granularities <- setdiff(parcel_granularities, "pass-through")
  list(
    roi_mask = parcel_gather ||
      all(vapply(caps, function(c) isTRUE(c$semantic$roi_mask), logical(1))),
    time_idx = all(vapply(caps, function(c) isTRUE(c$semantic$time_idx), logical(1))),
    parcel_idx = all(vapply(caps, function(c) isTRUE(c$semantic$parcel_idx), logical(1))),
    physical_roi_mask = parcel_gather ||
      all(vapply(caps, function(c) isTRUE(c$physical$roi_mask), logical(1))),
    physical_time_idx = all(vapply(caps, function(c) isTRUE(c$physical$time_idx), logical(1))),
    physical_parcel_idx = all(vapply(caps, function(c) isTRUE(c$physical$parcel_idx), logical(1))),
    temporal_granularity = if (length(temporal_granularities)) {
      paste(temporal_granularities, collapse = "+")
    } else {
      "pass-through"
    },
    parcel_granularity = if (length(parcel_granularities)) {
      paste(parcel_granularities, collapse = "+")
    } else {
      "pass-through"
    },
    roi_strategy = roi_strategy,
    coefficient_pipeline = roi_strategy %in% c("coefficient", "shared-coefficient"),
    dense_passthrough_pipeline = roi_strategy %in% c("dense-pass-through", "mixed")
  )
}
.basis_step_has_following_embed <- function(transforms, idx) {
  next_basis <- which(transforms == "basis" & seq_along(transforms) > idx)
  stop_idx <- if (length(next_basis) > 0L) {
    next_basis[[1L]] - 1L
  } else {
    length(transforms)
  }
  if (stop_idx <= idx) return(FALSE)
  any(startsWith(transforms[seq.int(idx + 1L, stop_idx)], "embed"))
}
