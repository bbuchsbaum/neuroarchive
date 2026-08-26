# Pure, bidser-independent planning for BIDS-shaped LNA exchange datasets.

.lna_dataset_plan_schema_version <- "1.0.0"
.lna_bids_entity_order <- c(
  "sub", "ses", "task", "acq", "ce", "rec", "dir", "run", "echo",
  "part", "chunk", "space", "res", "den", "desc", "hemi", "label",
  "from", "to", "mode", "recording", "proc", "split"
)

.dataset_plan_abort <- function(message, field, ...) {
  abort_lna(
    message,
    ...,
    field = field,
    .subclass = "lna_error_dataset_plan",
    location = paste0("plan_lna_dataset:", field)
  )
}

.dataset_plan_scalar_string <- function(x, field) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    .dataset_plan_abort(paste0(field, " must be one non-empty string"), field)
  }
  enc2utf8(x)
}

.dataset_plan_path_component <- function(x, field, bids_label = FALSE) {
  x <- .dataset_plan_scalar_string(x, field)
  pattern <- if (bids_label) {
    "^[A-Za-z0-9][A-Za-z0-9+.-]*$"
  } else {
    "^[A-Za-z0-9][A-Za-z0-9+.-]*$"
  }
  if (!grepl(pattern, x)) {
    .dataset_plan_abort(
      paste0(field, " is not a safe portable path component: ", x),
      field
    )
  }
  x
}

.dataset_plan_join <- function(...) {
  parts <- unlist(list(...), use.names = FALSE)
  parts <- parts[!is.na(parts) & nzchar(parts)]
  paste(parts, collapse = "/")
}

.dataset_plan_entity_names <- function(entities) {
  known <- intersect(.lna_bids_entity_order, names(entities))
  extra <- sort(setdiff(names(entities), .lna_bids_entity_order))
  c(known, extra)
}

.dataset_plan_bids_filename <- function(record) {
  entities <- record$entities
  if (is.null(entities$sub)) {
    .dataset_plan_abort(
      paste0("record ", record$record_id, " is missing the sub entity"),
      "entities.sub",
      record_id = record$record_id
    )
  }
  reserved <- intersect(
    names(entities), c("suffix", "datatype", "extension", "scope", "pipeline")
  )
  if (length(reserved)) {
    .dataset_plan_abort(
      paste0("reserved layout field used as entity: ", paste(reserved, collapse = ", ")),
      "entities",
      record_id = record$record_id
    )
  }
  order <- .dataset_plan_entity_names(entities)
  segments <- vapply(order, function(entity) {
    value <- .dataset_plan_path_component(
      entities[[entity]], paste0("entities.", entity), bids_label = TRUE
    )
    paste0(entity, "-", value)
  }, character(1))
  suffix <- .dataset_plan_path_component(record$bids$suffix, "bids.suffix")
  paste0(paste(c(segments, suffix), collapse = "_"), ".lna.h5")
}

.dataset_plan_relative_path <- function(record) {
  entities <- record$entities
  subject <- paste0(
    "sub-",
    .dataset_plan_path_component(entities$sub, "entities.sub", bids_label = TRUE)
  )
  session <- if (is.null(entities$ses)) {
    character()
  } else {
    paste0(
      "ses-",
      .dataset_plan_path_component(entities$ses, "entities.ses", bids_label = TRUE)
    )
  }
  datatype <- .dataset_plan_path_component(record$bids$datatype, "bids.datatype")
  filename <- .dataset_plan_bids_filename(record)
  scope_prefix <- if (identical(record$bids$scope, "derivative")) {
    pipeline <- record$bids$pipeline %||% "unspecified-pipeline"
    c(
      "derivatives",
      .dataset_plan_path_component(pipeline, "bids.pipeline")
    )
  } else {
    character()
  }
  .dataset_plan_join(scope_prefix, subject, session, datatype, filename)
}

.dataset_plan_codec <- function(codec_profile) {
  if (is.character(codec_profile) && length(codec_profile) == 1L) {
    codec_profile <- list(
      profile_id = codec_profile,
      loss_mode = "unspecified",
      transforms = character(),
      transform_params = list(),
      error_contract = lna_error_contract("unspecified")
    )
  }
  if (!is.list(codec_profile) || is.data.frame(codec_profile)) {
    .dataset_plan_abort(
      "codec_profile must be a profile ID or named list",
      "codec_profile"
    )
  }
  profile_id <- .dataset_plan_scalar_string(
    codec_profile$profile_id, "codec_profile.profile_id"
  )
  loss_mode <- .dataset_plan_scalar_string(
    codec_profile$loss_mode %||% "unspecified",
    "codec_profile.loss_mode"
  )
  if (!loss_mode %in% c("lossless", "lossy", "unspecified")) {
    .dataset_plan_abort(
      "codec_profile.loss_mode must be lossless, lossy, or unspecified",
      "codec_profile.loss_mode"
    )
  }
  transforms <- codec_profile$transforms %||% character()
  if (is.list(transforms) && is.null(names(transforms))) {
    transforms <- as.character(unlist(
      transforms, recursive = TRUE, use.names = FALSE
    ))
  }
  if (!is.character(transforms) || anyNA(transforms) || any(!nzchar(transforms))) {
    .dataset_plan_abort(
      "codec_profile.transforms must be a character vector",
      "codec_profile.transforms"
    )
  }
  transform_params <- codec_profile$transform_params %||% list()
  if (!is.list(transform_params)) {
    .dataset_plan_abort(
      "codec_profile.transform_params must be a list",
      "codec_profile.transform_params"
    )
  }
  tryCatch(
    .source_record_json_value(
      transform_params, "codec_profile.transform_params"
    ),
    lna_error_source_record = function(e) {
      .dataset_plan_abort(conditionMessage(e), "codec_profile.transform_params")
    }
  )
  error_contract <- codec_profile$error_contract
  if (is.null(error_contract)) {
    if (identical(loss_mode, "lossy")) {
      .dataset_plan_abort(
        "lossy codec profiles require an error_contract before writing",
        "codec_profile.error_contract"
      )
    }
    error_contract <- lna_error_contract(loss_mode)
  } else {
    error_contract <- tryCatch(
      validate_lna_error_contract(error_contract),
      lna_error_numerical_restoration = function(e) {
        .dataset_plan_abort(conditionMessage(e), "codec_profile.error_contract")
      }
    )
  }
  if (!identical(error_contract$loss_mode, loss_mode)) {
    .dataset_plan_abort(
      "codec loss_mode disagrees with its error_contract",
      "codec_profile.error_contract"
    )
  }
  list(
    profile_id = profile_id,
    loss_mode = loss_mode,
    transforms = transforms,
    transform_params = transform_params,
    error_contract = unclass(error_contract)
  )
}

.dataset_plan_companions <- function(record) {
  companions <- lapply(record$associations, function(association) {
    list(
      association_id = association$association_id,
      role = association$role,
      target_record_id = association$target_record_id,
      source_path = association$path,
      source_uri = association$uri,
      source_checksum = association$checksum,
      destination_relative = association$path,
      metadata = association$metadata
    )
  })
  if (!length(companions)) return(list())
  ids <- vapply(companions, `[[`, character(1), "association_id")
  companions[order(ids)]
}

.dataset_plan_entry <- function(record, output_root, codec_profile) {
  destination_relative <- .dataset_plan_relative_path(record)
  geometry <- record$geometry
  list(
    artifact_id = record$record_id,
    record_id = record$record_id,
    source = record$source,
    destination_relative = destination_relative,
    destination = file.path(
      output_root,
      do.call(file.path, as.list(strsplit(destination_relative, "/", fixed = TRUE)[[1L]]))
    ),
    bids = record$bids,
    entities = record$entities,
    resolved_metadata = record$resolved_metadata,
    expected_companions = .dataset_plan_companions(record),
    geometry = list(
      grid_id = geometry$grid_id,
      support_id = geometry$support_id,
      dims = geometry$dims,
      affine = geometry$affine,
      n_volumes = geometry$n_volumes,
      spatial_units = geometry$units,
      voxel_ordering = geometry$voxel_ordering
    ),
    reference_assertion = record$reference_assertion,
    codec_profile_id = codec_profile$profile_id,
    codec_transform_params = codec_profile$transform_params,
    shared_basis = list(status = "none")
  )
}

.dataset_plan_hash <- function(payload) {
  json <- jsonlite::toJSON(
    payload,
    auto_unbox = TRUE,
    null = "null",
    digits = NA,
    matrix = "rowmajor",
    pretty = FALSE
  )
  paste0(
    "sha256:",
    digest::digest(charToRaw(enc2utf8(json)), algo = "sha256", serialize = FALSE)
  )
}

.dataset_plan_check_collisions <- function(entries) {
  destinations <- vapply(entries, `[[`, character(1), "destination_relative")
  if (anyDuplicated(destinations)) {
    collision_paths <- unique(destinations[duplicated(destinations)])
    record_ids <- vapply(entries, `[[`, character(1), "record_id")
    involved <- record_ids[destinations %in% collision_paths]
    .dataset_plan_abort(
      paste0(
        "destination collision before writing: ",
        paste(collision_paths, collapse = ", "),
        " (records: ", paste(involved, collapse = ", "), ")"
      ),
      "destination_relative",
      collision_paths = collision_paths,
      record_ids = involved
    )
  }
  invisible(TRUE)
}

#' Plan a BIDS-shaped LNA exchange dataset
#'
#' This pure planning function accepts only neuroarchive-owned source records;
#' it does not parse filenames or call bidser. Every source image maps to one
#' `.lna.h5` artifact. All entity comparisons are exact and all destination
#' collisions are rejected before any filesystem mutation.
#'
#' @param records One `lna_source_record` or a list of source records.
#' @param output_root Destination dataset root. The planner does not create it.
#' @param codec_profile Profile ID or named list with `profile_id`, optional
#'   `loss_mode`, `transforms`, and `transform_params`.
#' @param fit_resources Optional shared-resource fit specifications. The fit
#'   phase is explicit even when this list is empty.
#' @param metadata_policy Portable metadata, association, and privacy policy.
#'   Pure planning permits unrecorded provenance by default; the transactional
#'   writer defaults to the strict portable profile.
#' @return An immutable `lna_dataset_plan` with deterministic artifacts and
#'   `plan_hash`.
#' @section Low-level multi-run mode:
#' Dataset planning deliberately uses one source image per LNA artifact. The
#' existing [write_lna()] list-input API remains the explicit, non-default
#' low-level route for a multi-run-in-one-file archive.
#' @export
plan_lna_dataset <- function(
    records,
    output_root,
    codec_profile = list(
      profile_id = "lna-independent-lossless-v1",
      loss_mode = "lossless",
      transforms = character(),
      transform_params = list()
    ),
    fit_resources = list(),
    metadata_policy = lna_metadata_policy()) {
  output_root <- .dataset_plan_scalar_string(output_root, "output_root")
  normalized <- validate_lna_source_records(records)
  if (inherits(normalized, "lna_source_record")) normalized <- list(normalized)
  codec_profile <- .dataset_plan_codec(codec_profile)

  entries <- lapply(normalized, function(record) {
    .dataset_plan_entry(record, output_root, codec_profile)
  })
  order <- order(
    vapply(entries, `[[`, character(1), "destination_relative"),
    vapply(entries, `[[`, character(1), "record_id")
  )
  entries <- entries[order]
  .dataset_plan_check_collisions(entries)
  metadata <- .metadata_apply_plan(entries, metadata_policy)
  entries <- metadata$entries
  entries <- .dataset_apply_shared_basis_resources(
    entries, fit_resources, codec_profile
  )

  payload <- list(
    schema_version = .lna_dataset_plan_schema_version,
    archive_mode = "one-source-one-lna",
    output_root = output_root,
    codec_profile = codec_profile,
    metadata_policy = metadata$policy,
    artifacts = entries,
    execution = .build_lna_execution(entries, fit_resources),
    summary = list(
      n_records = length(normalized),
      n_artifacts = length(entries),
      n_companions = sum(vapply(
        entries, function(entry) length(entry$expected_companions), integer(1)
      ))
    )
  )
  payload$plan_hash <- .dataset_plan_hash(payload)
  validate_lna_dataset_plan(
    structure(payload, class = c("lna_dataset_plan", "list"))
  )
}

#' @export
as.data.frame.lna_dataset_plan <- function(x, ...) {
  do.call(rbind, lapply(x$artifacts, function(entry) {
    data.frame(
      record_id = entry$record_id,
      source_path = entry$source$path %||% NA_character_,
      source_uri = entry$source$uri %||% NA_character_,
      source_checksum = entry$source$checksum,
      destination_relative = entry$destination_relative,
      destination = entry$destination,
      datatype = entry$bids$datatype,
      suffix = entry$bids$suffix,
      scope = entry$bids$scope,
      pipeline = entry$bids$pipeline %||% NA_character_,
      grid_id = entry$geometry$grid_id,
      support_id = entry$geometry$support_id %||% NA_character_,
      codec_profile_id = entry$codec_profile_id,
      stringsAsFactors = FALSE
    )
  }))
}

#' @export
print.lna_dataset_plan <- function(x, ...) {
  cat("<lna_dataset_plan>", x$plan_hash, "\n")
  cat("  root:", x$output_root, "\n")
  cat("  records:", x$summary$n_records, "artifacts:", x$summary$n_artifacts, "\n")
  cat("  fit resources:", length(x$execution$fit_resources), "\n")
  print(as.data.frame(x)[, c(
    "record_id", "destination_relative", "grid_id", "codec_profile_id"
  ), drop = FALSE], row.names = FALSE)
  invisible(x)
}
