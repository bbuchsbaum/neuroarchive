# Transactional generation updates for BIDS-shaped LNA datasets.

.lna_dataset_update_schema_version <- "1.0.0"

.update_abort <- function(message, field, ...) {
  abort_lna(
    message,
    ...,
    field = field,
    .subclass = "lna_error_dataset_update",
    location = paste0("update_lna_dataset:", field)
  )
}

.update_normalize_root <- function(path, field, must_work = TRUE) {
  path <- .execution_scalar_string(path, field)
  normalizePath(path, winslash = "/", mustWork = must_work)
}

.update_entry_record <- function(entry) {
  dims <- as.integer(unlist(entry$geometry$dims, recursive = TRUE, use.names = FALSE))
  affine <- entry$geometry$affine
  if (is.list(affine)) {
    rows <- length(affine)
    affine <- matrix(
      as.numeric(unlist(affine, recursive = TRUE, use.names = FALSE)),
      nrow = rows, byrow = TRUE
    )
  }
  n_volumes <- entry$geometry$n_volumes
  if (is.list(n_volumes)) {
    n_volumes <- unlist(n_volumes, recursive = TRUE, use.names = FALSE)
  }
  associations <- lapply(entry$expected_companions, function(item) {
    list(
      association_id = item$association_id,
      role = item$role,
      target_record_id = item$target_record_id,
      path = item$source_path,
      uri = item$source_uri,
      checksum = item$source_checksum,
      metadata = item$metadata
    )
  })
  lna_source_record(
    record_id = entry$record_id,
    source_path = entry$source$path,
    source_uri = entry$source$uri,
    source_checksum = entry$source$checksum,
    entities = entry$entities,
    datatype = entry$bids$datatype,
    suffix = entry$bids$suffix,
    extension = entry$bids$extension,
    scope = entry$bids$scope,
    pipeline = entry$bids$pipeline,
    resolved_metadata = entry$resolved_metadata,
    associations = associations,
    geometry = list(
      dims = dims,
      affine = affine,
      units = entry$geometry$spatial_units,
      voxel_ordering = entry$geometry$voxel_ordering,
      n_volumes = n_volumes,
      grid_id = entry$geometry$grid_id,
      support_id = entry$geometry$support_id
    ),
    reference_assertion = entry$reference_assertion
  )
}

.update_records_by_id <- function(records) {
  stats::setNames(records, vapply(records, `[[`, character(1), "record_id"))
}

.update_manifest_rows <- function(manifest) {
  rows <- split(manifest$data, seq_len(nrow(manifest$data)))
  stats::setNames(rows, vapply(rows, function(row) row$record_id[[1L]], character(1)))
}

.update_scalar <- function(row, field) {
  value <- row[[field]][[1L]]
  if (!length(value) || is.na(value)) NULL else value
}

.update_artifact_integrity <- function(root, row) {
  relative <- .update_scalar(row, "lna_path")
  expected_checksum <- .update_scalar(row, "lna_checksum")
  expected_bytes <- .update_scalar(row, "lna_bytes")
  path <- if (is.null(relative)) NULL else .writer_full_path(root, relative)
  if (is.null(path) || !file.exists(path)) {
    return(list(status = "corrupt", reason = "artifact-missing"))
  }
  actual_bytes <- unname(file.info(path)$size)
  if (!is.null(expected_bytes) && !identical(as.numeric(actual_bytes), as.numeric(expected_bytes))) {
    return(list(
      status = "corrupt", reason = "artifact-size-mismatch",
      expected_bytes = as.numeric(expected_bytes), actual_bytes = as.numeric(actual_bytes)
    ))
  }
  actual_checksum <- .journal_file_checksum(path)
  if (is.null(expected_checksum) || !identical(actual_checksum, expected_checksum)) {
    return(list(
      status = "corrupt", reason = "artifact-checksum-mismatch",
      expected_checksum = expected_checksum, actual_checksum = actual_checksum
    ))
  }
  list(
    status = "valid", reason = NULL, path = path,
    checksum = actual_checksum, bytes = as.numeric(actual_bytes)
  )
}

.update_semantic_payload <- function(entry, codec_profile, metadata_policy) {
  list(
    source_path = entry$source$path,
    source_uri = entry$source$uri,
    bids = entry$bids,
    entities = entry$entities,
    resolved_metadata = entry$resolved_metadata,
    expected_companions = entry$expected_companions,
    geometry = entry$geometry,
    reference_assertion = entry$reference_assertion,
    codec_profile = codec_profile,
    metadata_policy = metadata_policy
  )
}

.update_same_semantics <- function(current, incoming, current_plan, incoming_plan) {
  identical(
    .dataset_plan_hash(.update_semantic_payload(
      current, current_plan$codec_profile, current_plan$metadata_policy
    )),
    .dataset_plan_hash(.update_semantic_payload(
      incoming, incoming_plan$codec_profile, incoming_plan$metadata_policy
    ))
  )
}

.update_revision <- function(root) {
  path <- file.path(root, "lna_update.json")
  if (!file.exists(path)) return(1L)
  receipt <- read_lna_dataset_update(path)
  as.integer(receipt$revision)
}

.update_plan_payload <- function(plan) {
  payload <- unclass(plan)
  payload$update_plan_hash <- NULL
  payload$target_plan <- unclass(payload$target_plan)
  payload
}

.update_target_records <- function(current, incoming, classifications, mode, prune) {
  current_map <- .update_records_by_id(current)
  incoming_map <- .update_records_by_id(incoming)
  ids <- sort(unique(c(names(current_map), names(incoming_map))))
  out <- list()
  for (id in ids) {
    status <- classifications[[id]]$status
    if (identical(status, "added")) {
      out[[id]] <- incoming_map[[id]]
    } else if (identical(status, "source-missing")) {
      if (!identical(mode, "mirror") || !isTRUE(prune)) out[[id]] <- current_map[[id]]
    } else if (identical(mode, "append")) {
      out[[id]] <- current_map[[id]]
    } else {
      out[[id]] <- incoming_map[[id]]
    }
  }
  unname(out)
}

.update_action <- function(status, mode, prune) {
  if (identical(status, "added")) return("add")
  if (identical(status, "source-missing")) {
    if (identical(mode, "mirror") && isTRUE(prune)) return("remove")
    return("retain")
  }
  if (identical(mode, "append")) {
    if (identical(status, "corrupt")) return("blocked-corrupt")
    return("retain")
  }
  if (status %in% c("changed", "stale", "corrupt")) return("rebuild")
  "retain"
}

#' Plan an incremental LNA dataset update
#'
#' The returned preview is immutable and classifies every current or incoming
#' record as `unchanged`, `added`, `changed`, `source-missing`, `stale`, or
#' `corrupt`. `append` retains changed records and only adds new ones;
#' `refresh` rebuilds changed, stale, or corrupt records; `mirror` additionally
#' plans removals only when `prune = TRUE`.
#'
#' @param dataset_root Existing archive root.
#' @param records Complete current source inventory as source records.
#' @param source_root Root used to resolve the incoming records.
#' @param mode Update policy: `append`, `refresh`, or `mirror`.
#' @param prune Whether a mirror preview should include removals.
#' @param codec_profile Optional replacement codec profile. The current profile
#'   is retained when omitted.
#' @param metadata_policy Optional replacement metadata policy.
#' @param fit_resources Optional target fit resources. When omitted, current
#'   receipts are reused exactly; changed shared-basis membership or source
#'   checksums therefore require an explicit refit before refresh or mirror.
#' @return An `lna_dataset_update_plan` suitable for explicit review and
#'   [update_lna_dataset()].
#' @export
plan_lna_dataset_update <- function(
    dataset_root,
    records,
    source_root = ".",
    mode = c("append", "refresh", "mirror"),
    prune = FALSE,
    codec_profile = NULL,
    metadata_policy = NULL,
    fit_resources = NULL) {
  mode <- match.arg(mode)
  if (!is.logical(prune) || length(prune) != 1L || is.na(prune)) {
    .update_abort("prune must be TRUE or FALSE", "prune")
  }
  if (isTRUE(prune) && !identical(mode, "mirror")) {
    .update_abort("prune is only valid for mirror previews", "prune")
  }
  dataset_root <- .update_root_path(dataset_root, "dataset_root")
  .update_assert_plannable(dataset_root)
  source_root <- .update_normalize_root(source_root, "source_root")
  current_plan <- read_lna_dataset_plan(file.path(dataset_root, ".lna_plan.json"))
  manifest <- read_lna_manifest(dataset_root)
  declaration <- read_lna_dataset_declaration(dataset_root)
  validate_lna_dataset_declaration(declaration, current_plan$plan_hash)
  current_ids <- vapply(current_plan$artifacts, `[[`, character(1), "record_id")
  manifest_rows <- .update_manifest_rows(manifest)
  if (!setequal(current_ids, names(manifest_rows))) {
    .update_abort("manifest and serialized plan disagree on record identities", "manifest")
  }
  normalized <- validate_lna_source_records(records)
  if (inherits(normalized, "lna_source_record")) normalized <- list(normalized)
  codec_profile <- codec_profile %||% current_plan$codec_profile
  metadata_policy <- metadata_policy %||% current_plan$metadata_policy
  fit_resources <- fit_resources %||% current_plan$execution$fit_resources
  incoming_plan <- plan_lna_dataset(
    normalized, dataset_root, codec_profile = codec_profile,
    fit_resources = fit_resources, metadata_policy = metadata_policy
  )
  current_records <- lapply(current_plan$artifacts, .update_entry_record)
  current_map <- stats::setNames(
    current_plan$artifacts,
    vapply(current_plan$artifacts, `[[`, character(1), "record_id")
  )
  incoming_map <- stats::setNames(
    incoming_plan$artifacts,
    vapply(incoming_plan$artifacts, `[[`, character(1), "record_id")
  )
  ids <- sort(unique(c(names(current_map), names(incoming_map))))
  classifications <- stats::setNames(lapply(ids, function(id) {
    current <- current_map[[id]]
    incoming <- incoming_map[[id]]
    integrity <- if (is.null(current)) {
      list(status = "not-applicable", reason = NULL)
    } else {
      .update_artifact_integrity(dataset_root, manifest_rows[[id]])
    }
    status <- if (is.null(current)) {
      "added"
    } else if (is.null(incoming)) {
      "source-missing"
    } else if (identical(integrity$status, "corrupt")) {
      "corrupt"
    } else if (!identical(current$source$checksum, incoming$source$checksum)) {
      "changed"
    } else if (!.update_same_semantics(current, incoming, current_plan, incoming_plan)) {
      "stale"
    } else {
      "unchanged"
    }
    list(
      record_id = id,
      status = status,
      action = .update_action(status, mode, prune),
      current_source_checksum = if (is.null(current)) NULL else current$source$checksum,
      incoming_source_checksum = if (is.null(incoming)) NULL else incoming$source$checksum,
      integrity_status = integrity$status,
      integrity_reason = integrity$reason
    )
  }), ids)
  if (identical(mode, "append") &&
      !identical(current_plan$codec_profile, incoming_plan$codec_profile)) {
    .update_abort("append cannot change the archive codec profile; use refresh", "codec_profile")
  }
  source_missing <- names(classifications)[vapply(classifications, function(x) {
    identical(x$status, "source-missing")
  }, logical(1))]
  profile_changed <- !identical(current_plan$codec_profile, incoming_plan$codec_profile) ||
    !identical(current_plan$metadata_policy, incoming_plan$metadata_policy)
  if (identical(mode, "refresh") && profile_changed && length(source_missing)) {
    .update_abort(
      "refreshing a codec or metadata profile requires a complete source inventory",
      "records", source_missing_record_ids = source_missing
    )
  }
  target_records <- .update_target_records(
    current_records, normalized, classifications, mode, prune
  )
  target_plan <- plan_lna_dataset(
    target_records, dataset_root, codec_profile = codec_profile,
    fit_resources = fit_resources, metadata_policy = metadata_policy
  )
  revision <- .update_revision(dataset_root) + 1L
  payload <- list(
    schema_version = .lna_dataset_update_schema_version,
    dataset_root = dataset_root,
    source_root = source_root,
    revision = revision,
    mode = mode,
    prune = isTRUE(prune),
    base_plan_hash = current_plan$plan_hash,
    target_plan_hash = target_plan$plan_hash,
    classifications = unname(classifications),
    target_plan = unclass(target_plan)
  )
  payload$update_plan_hash <- .dataset_plan_hash(payload)
  structure(payload, class = c("lna_dataset_update_plan", "list"))
}

#' @export
as.data.frame.lna_dataset_update_plan <- function(x, ...) {
  do.call(rbind, lapply(x$classifications, function(item) {
    data.frame(
      record_id = item$record_id,
      status = item$status,
      action = item$action,
      integrity_status = item$integrity_status,
      stringsAsFactors = FALSE
    )
  }))
}

#' @export
print.lna_dataset_update_plan <- function(x, ...) {
  cat("<lna_dataset_update_plan>", x$update_plan_hash, "\n")
  cat("  revision:", x$revision, "mode:", x$mode, "prune:", x$prune, "\n")
  print(as.data.frame(x), row.names = FALSE)
  invisible(x)
}

.update_validate_plan <- function(plan) {
  if (!inherits(plan, "lna_dataset_update_plan") || !is.list(plan)) {
    .update_abort("plan must be an lna_dataset_update_plan preview", "plan")
  }
  required <- c(
    "schema_version", "dataset_root", "source_root", "revision", "mode", "prune",
    "base_plan_hash", "target_plan_hash", "classifications", "target_plan",
    "update_plan_hash"
  )
  if (!setequal(names(plan), required) ||
      !identical(plan$schema_version, .lna_dataset_update_schema_version)) {
    .update_abort("update plan schema is incompatible", "plan")
  }
  target <- validate_lna_dataset_plan(plan$target_plan)
  if (!identical(target$plan_hash, plan$target_plan_hash) ||
      !identical(.dataset_plan_hash(.update_plan_payload(plan)), plan$update_plan_hash)) {
    .update_abort("update plan hash verification failed", "plan")
  }
  plan$target_plan <- target
  structure(plan, class = c("lna_dataset_update_plan", "list"))
}

.update_journal_append <- function(root, plan, state, details = list()) {
  path <- file.path(root, ".lna_update_journal.jsonl")
  events <- if (file.exists(path)) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    lapply(lines[nzchar(lines)], jsonlite::fromJSON, simplifyVector = FALSE)
  } else list()
  sequence <- length(events) + 1L
  event <- list(
    schema_version = "1.0.0",
    sequence = sequence,
    update_plan_hash = plan$update_plan_hash,
    archive_plan_hash = plan$target_plan_hash,
    revision = plan$revision,
    mode = plan$mode,
    state = state,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
    details = details
  )
  event$event_id <- .dataset_plan_hash(event)
  line <- jsonlite::toJSON(
    event, auto_unbox = TRUE, null = "null", digits = NA, pretty = FALSE
  )
  cat(line, "\n", file = path, append = TRUE, sep = "")
  invisible(event)
}

.update_write_receipt <- function(root, plan, previous_generation) {
  counts <- table(vapply(plan$classifications, `[[`, character(1), "status"))
  receipt <- list(
    schema_version = .lna_dataset_update_schema_version,
    revision = as.numeric(plan$revision),
    update_plan_hash = plan$update_plan_hash,
    base_plan_hash = plan$base_plan_hash,
    archive_plan_hash = plan$target_plan_hash,
    mode = plan$mode,
    prune = plan$prune,
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
    classifications = as.list(as.integer(counts)),
    previous_generation = previous_generation
  )
  names(receipt$classifications) <- names(counts)
  .update_atomic_json(receipt, file.path(root, "lna_update.json"))
  receipt
}

#' Read the latest LNA dataset update receipt
#'
#' @param path Archive root or `lna_update.json` path.
#' @return The validated update receipt.
#' @export
read_lna_dataset_update <- function(path) {
  path <- .execution_scalar_string(path, "path")
  if (!identical(basename(path), "lna_update.json")) path <- file.path(path, "lna_update.json")
  if (!file.exists(path)) .update_abort("update receipt does not exist", "path")
  receipt <- tryCatch(
    jsonlite::read_json(path, simplifyVector = FALSE),
    error = function(e) .update_abort("update receipt is invalid JSON", "path")
  )
  required <- c(
    "schema_version", "revision", "update_plan_hash", "base_plan_hash",
    "archive_plan_hash", "mode", "prune", "completed_at", "classifications",
    "previous_generation"
  )
  if (!is.list(receipt) || !setequal(names(receipt), required) ||
      !identical(receipt$schema_version, .lna_dataset_update_schema_version)) {
    .update_abort("update receipt schema is incompatible", "path")
  }
  receipt$update_plan_hash <- .source_record_checksum(
    receipt$update_plan_hash, "receipt.update_plan_hash"
  )
  receipt$base_plan_hash <- .source_record_checksum(
    receipt$base_plan_hash, "receipt.base_plan_hash"
  )
  receipt$archive_plan_hash <- .source_record_checksum(
    receipt$archive_plan_hash, "receipt.archive_plan_hash"
  )
  receipt
}

.update_hooks <- function(plan, current_manifest, user_hooks) {
  hooks <- .writer_hooks(user_hooks)
  dataset_root <- plan$dataset_root
  source_root <- plan$source_root
  actions <- stats::setNames(
    vapply(plan$classifications, `[[`, character(1), "action"),
    vapply(plan$classifications, `[[`, character(1), "record_id")
  )
  rows <- .update_manifest_rows(current_manifest)
  reuse <- function(entry) identical(actions[[entry$record_id]], "retain")
  user_verify <- hooks$verify_source
  user_fit <- hooks$fit_resource
  user_load <- hooks$load_source
  user_write <- hooks$write_artifact
  user_participants <- hooks$write_participants
  hooks$verify_source <- function(entry, ignored_source_root) {
    if (!reuse(entry)) return(user_verify(entry, source_root))
    row <- rows[[entry$record_id]]
    integrity <- .update_artifact_integrity(dataset_root, row)
    if (!identical(integrity$status, "valid")) {
      .update_abort(
        paste0("retained artifact is not valid: ", integrity$reason),
        "artifact", record_id = entry$record_id
      )
    }
    integrity$path
  }
  hooks$fit_resource <- function(resource, target_plan, archive_root) {
    if (!identical(resource$resource_type, "shared-basis")) {
      return(user_fit(resource, target_plan, archive_root))
    }
    p <- .dataset_shared_basis_parameters(resource)
    destination_dir <- file.path(archive_root, "bases")
    destination <- file.path(destination_dir, p$basis_filename)
    if (!file.exists(destination)) {
      candidates <- c(
        file.path(dataset_root, "bases", p$basis_filename),
        file.path(source_root, "bases", p$basis_filename)
      )
      candidate <- candidates[file.exists(candidates)][1L]
      if (is.na(candidate)) {
        .update_abort("shared-basis update resource is unavailable", "fit_resources")
      }
      dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)
      if (!isTRUE(file.copy(candidate, destination, overwrite = FALSE, copy.mode = TRUE))) {
        .update_abort("could not stage shared-basis update resource", "fit_resources")
      }
      basis <- read_shared_basis(destination, verify_checksum = TRUE)
      .update_basis_registry(destination_dir, p$basis_id, list(
        checksum = p$basis_checksum,
        filename = p$basis_filename,
        kind = p$basis_kind,
        mask_checksum = basis$mask_checksum,
        grid_id = p$grid_id,
        support_id = p$support_id,
        support_type = p$support_type,
        n_atoms = as.integer(ncol(basis$loadings)),
        n_voxels = as.integer(nrow(basis$loadings)),
        created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
      ))
    }
    user_fit(resource, target_plan, archive_root)
  }
  hooks$load_source <- function(entry, source_path) {
    if (reuse(entry)) {
      return(structure(list(path = source_path), class = "lna_reused_artifact"))
    }
    user_load(entry, source_path)
  }
  hooks$write_artifact <- function(data, temp_path, entry, target_plan) {
    if (inherits(data, "lna_reused_artifact")) {
      if (!isTRUE(file.copy(data$path, temp_path, overwrite = FALSE, copy.mode = TRUE))) {
        .update_abort("could not stage retained artifact", "artifact", record_id = entry$record_id)
      }
      return(invisible(temp_path))
    }
    user_write(data, temp_path, entry, target_plan)
  }
  hooks$resolve_companion <- function(companion, ignored_source_root, archive_root) {
    archived <- .writer_full_path(dataset_root, companion$destination_relative)
    if (file.exists(archived) &&
        identical(.journal_file_checksum(archived), companion$source_checksum)) {
      return(archived)
    }
    .writer_full_path(source_root, companion$source_path)
  }
  hooks$write_participants <- function(target_plan, ignored_source_root, output_root, overwrite) {
    policy <- target_plan$metadata_policy$participants
    if (identical(policy$mode, "exclude")) return(NULL)
    source <- .writer_full_path(source_root, policy$source_path)
    if (file.exists(source)) {
      return(user_participants(target_plan, source_root, output_root, overwrite))
    }
    declaration <- read_lna_dataset_declaration(dataset_root)
    receipt <- declaration$participants_receipt
    archived <- file.path(dataset_root, "participants.tsv")
    if (is.null(receipt) || !file.exists(archived) ||
        !identical(receipt$source_checksum, policy$source_checksum) ||
        !identical(.journal_file_checksum(archived), receipt$output_checksum)) {
      .update_abort("participants source is unavailable and no valid archived subset exists", "participants")
    }
    destination <- file.path(output_root, "participants.tsv")
    if (!isTRUE(file.copy(archived, destination, overwrite = FALSE, copy.mode = TRUE))) {
      .update_abort("could not stage archived participants subset", "participants")
    }
    receipt
  }
  hooks
}

.update_timestamp <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC")
}

.update_root_path <- function(path, field = "dataset_root") {
  path <- .execution_scalar_string(path, field)
  if (dir.exists(path)) {
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }
  parent <- normalizePath(dirname(path), winslash = "/", mustWork = TRUE)
  file.path(parent, basename(path))
}

.update_parent_paths <- function(root) {
  parent <- dirname(root)
  stem <- paste0(".", basename(root), ".lna-update")
  list(
    lock = file.path(parent, paste0(stem, ".lock")),
    state = file.path(parent, paste0(stem, "-state")),
    audit = file.path(parent, paste0(stem, "-recovery.jsonl"))
  )
}

.update_assert_plannable <- function(root) {
  paths <- .update_parent_paths(root)
  if (dir.exists(paths$state) || dir.exists(paths$lock)) {
    .update_abort(
      paste0(
        "an unfinished update transaction must be recovered first; retry the ",
        "reviewed update_lna_dataset() call and set takeover_stale = TRUE only ",
        "after auditing a stale owner"
      ),
      "recovery"
    )
  }
  if (!dir.exists(root)) {
    .update_abort("dataset root is missing and no recoverable transaction was found", "dataset_root")
  }
  invisible(root)
}

.update_atomic_json <- function(value, path) {
  temporary <- tempfile(
    pattern = paste0(".", basename(path), ".tmp-"), tmpdir = dirname(path)
  )
  on.exit(if (file.exists(temporary)) unlink(temporary, force = TRUE), add = TRUE)
  jsonlite::write_json(
    value, temporary, auto_unbox = TRUE, null = "null", digits = NA,
    pretty = TRUE
  )
  if (!file.rename(temporary, path)) {
    .update_abort("could not atomically publish update state", "transaction_state")
  }
  invisible(path)
}

.update_audit_append <- function(root, action, details = list(), plan = NULL) {
  path <- .update_parent_paths(root)$audit
  event <- list(
    schema_version = "1.0.0",
    action = action,
    timestamp = .update_timestamp(),
    update_plan_hash = if (is.null(plan)) NULL else plan$update_plan_hash,
    details = details
  )
  event$event_id <- .dataset_plan_hash(event)
  cat(
    jsonlite::toJSON(event, auto_unbox = TRUE, null = "null", digits = NA),
    "\n", file = path, append = TRUE, sep = ""
  )
  invisible(event)
}

.update_lock_owner <- function(owner, field = "lock.owner") {
  required <- c(
    "schema_version", "update_plan_hash", "revision", "pid", "host",
    "token", "acquired_at"
  )
  if (!is.list(owner) || is.data.frame(owner) ||
      !setequal(names(owner), required) ||
      !identical(owner$schema_version, "1.0.0")) {
    .update_abort("update lock owner record is malformed", field)
  }
  owner$update_plan_hash <- tryCatch(
    .source_record_checksum(owner$update_plan_hash, paste0(field, ".update_plan_hash")),
    error = function(e) .update_abort("update lock plan hash is malformed", field)
  )
  for (name in c("revision", "pid")) {
    value <- suppressWarnings(as.numeric(unlist(owner[[name]], use.names = FALSE)))
    if (length(value) != 1L || !is.finite(value) || value < 1 || value != floor(value)) {
      .update_abort(paste0("update lock ", name, " is malformed"), field)
    }
    owner[[name]] <- value
  }
  owner$host <- .execution_scalar_string(owner$host, paste0(field, ".host"))
  owner$token <- tryCatch(
    .source_record_checksum(owner$token, paste0(field, ".token")),
    error = function(e) .update_abort("update lock token is malformed", field)
  )
  for (name in "acquired_at") {
    value <- .execution_scalar_string(owner[[name]], paste0(field, ".", name))
    parsed <- suppressWarnings(as.POSIXct(
      sub("Z$", "", value), format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC"
    ))
    if (is.na(parsed)) .update_abort("update lock timestamp is malformed", field)
    owner[[name]] <- value
  }
  owner
}

.update_read_lock_owner <- function(lock) {
  path <- file.path(lock, "owner.json")
  if (!file.exists(path)) .update_abort("update lock has no owner record", "lock")
  decoded <- tryCatch(
    jsonlite::read_json(path, simplifyVector = FALSE),
    error = function(e) .update_abort("update lock owner is invalid JSON", "lock")
  )
  .update_lock_owner(decoded)
}

.update_pid_alive <- function(pid, host) {
  current_host <- unname(Sys.info()[["nodename"]])
  if (!identical(host, current_host)) return(NA)
  if (.Platform$OS.type == "windows") {
    return(tryCatch({
      tools::pskill(as.integer(pid), signal = 0L)
      TRUE
    }, error = function(e) FALSE))
  }
  status <- suppressWarnings(system2(
    "kill", c("-0", as.character(as.integer(pid))),
    stdout = FALSE, stderr = FALSE
  ))
  identical(as.integer(status), 0L)
}

.update_lock_age <- function(owner) {
  acquired <- as.POSIXct(
    sub("Z$", "", owner$acquired_at),
    format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC"
  )
  as.numeric(difftime(Sys.time(), acquired, units = "secs"))
}

.update_validate_stale_after <- function(value) {
  value <- suppressWarnings(as.numeric(value))
  if (length(value) != 1L || !is.finite(value) || value < 0) {
    .update_abort("stale_after_seconds must be one non-negative number", "stale_after_seconds")
  }
  value
}

.update_new_lock_owner <- function(plan) {
  acquired <- .update_timestamp()
  host <- unname(Sys.info()[["nodename"]])
  token <- .dataset_plan_hash(list(
    pid = Sys.getpid(), host = host, acquired_at = acquired,
    nonce = paste(sample(c(letters, LETTERS, 0:9), 32L, replace = TRUE), collapse = "")
  ))
  list(
    schema_version = "1.0.0",
    update_plan_hash = plan$update_plan_hash,
    revision = as.integer(plan$revision),
    pid = as.integer(Sys.getpid()),
    host = host,
    token = token,
    acquired_at = acquired
  )
}

.update_lock <- function(
    root,
    plan,
    takeover_stale = FALSE,
    stale_after_seconds = 3600) {
  if (!is.logical(takeover_stale) || length(takeover_stale) != 1L ||
      is.na(takeover_stale)) {
    .update_abort("takeover_stale must be TRUE or FALSE", "takeover_stale")
  }
  stale_after_seconds <- .update_validate_stale_after(stale_after_seconds)
  lock <- .update_parent_paths(root)$lock
  if (dir.exists(lock)) {
    prior <- .update_read_lock_owner(lock)
    alive <- .update_pid_alive(prior$pid, prior$host)
    age <- .update_lock_age(prior)
    if (isTRUE(alive)) {
      .update_abort(
        "the update lock owner is demonstrably live; takeover is forbidden",
        "lock", owner = prior
      )
    }
    if (!isTRUE(takeover_stale) || age < stale_after_seconds) {
      .update_abort(
        paste0(
          "an update lock remains from a non-live or unverified owner; inspect ",
          "the owner and retry with takeover_stale = TRUE after the stale threshold"
        ),
        "lock", owner = prior, age_seconds = age
      )
    }
    abandoned <- paste0(
      lock, ".abandoned-", format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"),
      "-", substr(sub("^sha256:", "", prior$token), 1L, 8L)
    )
    if (!file.rename(lock, abandoned)) {
      .update_abort("could not quarantine the stale update lock", "lock")
    }
    .update_audit_append(root, "stale-lock-takeover", list(
      prior_owner = prior,
      age_seconds = age,
      liveness = if (is.na(alive)) "unverified-remote-host" else "not-live",
      quarantined_lock = basename(abandoned)
    ), plan)
  }
  if (!dir.create(lock, recursive = FALSE, showWarnings = FALSE)) {
    .update_abort("another dataset update holds the advisory lock", "lock")
  }
  owner <- .update_new_lock_owner(plan)
  .update_atomic_json(owner, file.path(lock, "owner.json"))
  structure(list(path = lock, owner = owner), class = "lna_update_lock")
}

.update_release_lock <- function(lock) {
  if (!dir.exists(lock$path)) return(invisible(FALSE))
  owner <- .update_read_lock_owner(lock$path)
  if (!identical(owner$token, lock$owner$token)) {
    .update_abort("refusing to release an update lock owned by another process", "lock")
  }
  unlink(lock$path, recursive = TRUE, force = TRUE)
  invisible(!dir.exists(lock$path))
}

.lna_update_transaction_phases <- c(
  "locked", "staging", "staged", "preserving-old", "old-preserved",
  "installing-new", "new-installed", "before-completion-receipt",
  "receipt-written", "completed", "recovery-resume", "recovery-forward",
  "recovery-rollback"
)

.update_state_payload <- function(event) {
  event$event_id <- NULL
  event
}

.update_state_event <- function(event, field = "transaction_state") {
  required <- c(
    "schema_version", "sequence", "dataset_basename", "stage_basename",
    "backup_basename", "update_plan_hash", "base_plan_hash",
    "target_plan_hash", "revision", "phase", "timestamp", "owner_token",
    "previous_event_id", "event_id"
  )
  if (!is.list(event) || is.data.frame(event) ||
      !setequal(names(event), required) ||
      !identical(event$schema_version, "1.0.0")) {
    .update_abort("update transaction state is malformed", field)
  }
  for (name in c(
    "dataset_basename", "stage_basename", "backup_basename", "phase", "timestamp"
  )) {
    event[[name]] <- .execution_scalar_string(
      event[[name]], paste0(field, ".", name)
    )
  }
  for (name in c("dataset_basename", "stage_basename", "backup_basename")) {
    value <- event[[name]]
    if (!identical(basename(value), value) || value %in% c(".", "..") ||
        grepl("[\\\\/]", value)) {
      .update_abort("update transaction contains an unsafe sibling name", field)
    }
  }
  if (!event$phase %in% .lna_update_transaction_phases) {
    .update_abort("update transaction phase is unsupported", field)
  }
  parsed_timestamp <- suppressWarnings(as.POSIXct(
    sub("Z$", "", event$timestamp),
    format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC"
  ))
  if (!grepl("Z$", event$timestamp) || is.na(parsed_timestamp)) {
    .update_abort("update transaction timestamp is malformed", field)
  }
  for (name in c(
    "update_plan_hash", "base_plan_hash", "target_plan_hash", "owner_token",
    "event_id"
  )) {
    event[[name]] <- tryCatch(
      .source_record_checksum(event[[name]], paste0(field, ".", name)),
      error = function(e) .update_abort("update transaction checksum is malformed", field)
    )
  }
  if (!is.null(event$previous_event_id)) {
    event$previous_event_id <- tryCatch(
      .source_record_checksum(
        event$previous_event_id, paste0(field, ".previous_event_id")
      ),
      error = function(e) .update_abort("update transaction chain is malformed", field)
    )
  }
  for (name in c("sequence", "revision")) {
    value <- suppressWarnings(as.numeric(unlist(event[[name]], use.names = FALSE)))
    if (length(value) != 1L || !is.finite(value) || value < 1 || value != floor(value)) {
      .update_abort("update transaction sequence or revision is malformed", field)
    }
    event[[name]] <- value
  }
  expected <- .dataset_plan_hash(.update_state_payload(event))
  if (!identical(expected, event$event_id)) {
    .update_abort("update transaction event hash verification failed", field)
  }
  event
}

.update_read_state <- function(root, required = FALSE) {
  directory <- .update_parent_paths(root)$state
  if (!dir.exists(directory)) {
    if (isTRUE(required)) .update_abort("update transaction state is missing", "transaction_state")
    return(NULL)
  }
  files <- sort(list.files(
    directory, pattern = "^[0-9]{6}-[a-z-]+[.]json$", full.names = TRUE
  ))
  entries <- list.files(directory, all.files = TRUE, no.. = TRUE)
  unexpected <- setdiff(entries[!startsWith(entries, ".")], basename(files))
  if (length(unexpected)) {
    .update_abort("update transaction state directory is malformed", "transaction_state")
  }
  if (!length(files)) {
    return(list(directory = directory, events = list(), latest = NULL))
  }
  events <- lapply(seq_along(files), function(i) {
    decoded <- tryCatch(
      jsonlite::read_json(files[[i]], simplifyVector = FALSE),
      error = function(e) .update_abort("update transaction event is invalid JSON", "transaction_state")
    )
    .update_state_event(decoded)
  })
  sequences <- vapply(events, `[[`, numeric(1), "sequence")
  if (!identical(sequences, as.numeric(seq_along(events)))) {
    .update_abort("update transaction sequence is not contiguous", "transaction_state")
  }
  expected_files <- vapply(events, function(event) {
    sprintf("%06d-%s.json", event$sequence, event$phase)
  }, character(1))
  if (!identical(basename(files), expected_files)) {
    .update_abort("update transaction event filenames disagree with their content", "transaction_state")
  }
  identity_fields <- c(
    "dataset_basename", "stage_basename", "backup_basename", "update_plan_hash",
    "base_plan_hash", "target_plan_hash", "revision"
  )
  first <- events[[1L]]
  for (i in seq_along(events)) {
    event <- events[[i]]
    if (!all(vapply(identity_fields, function(name) {
      identical(event[[name]], first[[name]])
    }, logical(1)))) {
      .update_abort("update transaction identity changes between events", "transaction_state")
    }
    previous <- if (i == 1L) NULL else events[[i - 1L]]$event_id
    if (!identical(event$previous_event_id, previous)) {
      .update_abort("update transaction event chain is broken", "transaction_state")
    }
  }
  list(directory = directory, events = events, latest = events[[length(events)]])
}

.update_checkpoint <- function(root, plan, stage, backup, phase, lock, hook = NULL) {
  if (!phase %in% .lna_update_transaction_phases) {
    .update_abort("unknown update transaction checkpoint", "transaction_state")
  }
  state <- .update_read_state(root, required = FALSE)
  if (is.null(state)) {
    directory <- .update_parent_paths(root)$state
    if (!dir.create(directory, recursive = FALSE, showWarnings = FALSE)) {
      .update_abort("could not create parent update transaction state", "transaction_state")
    }
    state <- list(directory = directory, events = list(), latest = NULL)
  }
  previous <- state$latest
  identity <- list(
    dataset_basename = basename(root),
    stage_basename = basename(stage),
    backup_basename = basename(backup),
    update_plan_hash = plan$update_plan_hash,
    base_plan_hash = plan$base_plan_hash,
    target_plan_hash = plan$target_plan_hash,
    revision = as.numeric(plan$revision),
    owner_token = lock$owner$token
  )
  stable_identity <- setdiff(names(identity), "owner_token")
  if (!is.null(previous) && !all(vapply(stable_identity, function(name) {
    identical(previous[[name]], identity[[name]])
  }, logical(1)))) {
    .update_abort("existing update transaction belongs to another owner or plan", "transaction_state")
  }
  event <- c(list(
    schema_version = "1.0.0",
    sequence = length(state$events) + 1L
  ), identity[c("dataset_basename", "stage_basename", "backup_basename",
    "update_plan_hash", "base_plan_hash", "target_plan_hash", "revision")], list(
    phase = phase,
    timestamp = .update_timestamp(),
    owner_token = identity$owner_token,
    previous_event_id = if (is.null(previous)) NULL else previous$event_id
  ))
  event$event_id <- .dataset_plan_hash(.update_state_payload(event))
  filename <- sprintf("%06d-%s.json", event$sequence, phase)
  .update_atomic_json(event, file.path(state$directory, filename))
  if (!is.null(hook)) hook(phase = phase, state = event)
  invisible(event)
}

.update_generation_status <- function(path, expected_plan_hash) {
  if (!dir.exists(path)) return(list(exists = FALSE, valid = FALSE, reason = "missing"))
  tryCatch({
    plan <- read_lna_dataset_plan(file.path(path, ".lna_plan.json"))
    declaration <- read_lna_dataset_declaration(path)
    manifest <- read_lna_manifest(path)
    if (!identical(plan$plan_hash, expected_plan_hash)) {
      return(list(
        exists = TRUE, valid = FALSE, reason = "plan-hash-mismatch",
        plan_hash = plan$plan_hash
      ))
    }
    validate_lna_dataset_declaration(declaration, expected_plan_hash)
    if (!identical(
        declaration$manifest$schema_version,
        manifest$schema$manifest_schema_version)) {
      return(list(exists = TRUE, valid = FALSE, reason = "manifest-version-mismatch"))
    }
    if (!all(manifest$data$status == "complete")) {
      return(list(exists = TRUE, valid = FALSE, reason = "manifest-incomplete"))
    }
    rows <- .update_manifest_rows(manifest)
    integrity <- lapply(rows, function(row) .update_artifact_integrity(path, row))
    invalid <- names(integrity)[vapply(integrity, function(item) {
      !identical(item$status, "valid")
    }, logical(1))]
    if (length(invalid)) {
      return(list(
        exists = TRUE, valid = FALSE, reason = "artifact-integrity-failure",
        record_ids = invalid
      ))
    }
    list(exists = TRUE, valid = TRUE, reason = NULL, plan_hash = plan$plan_hash)
  }, error = function(e) list(
    exists = TRUE, valid = FALSE, reason = conditionMessage(e), error = e
  ))
}

.update_journal_completed <- function(root, plan_hash) {
  path <- file.path(root, ".lna_update_journal.jsonl")
  if (!file.exists(path)) return(FALSE)
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  events <- lapply(lines[nzchar(lines)], function(line) tryCatch(
    jsonlite::fromJSON(line, simplifyVector = FALSE), error = function(e) NULL
  ))
  any(vapply(events, function(event) {
    !is.null(event) && identical(event$update_plan_hash, plan_hash) &&
      identical(event$state, "completed")
  }, logical(1)))
}

.update_existing_receipt <- function(root, plan) {
  path <- file.path(root, "lna_update.json")
  if (!file.exists(path)) return(NULL)
  receipt <- read_lna_dataset_update(path)
  if (!identical(receipt$update_plan_hash, plan$update_plan_hash) ||
      !identical(receipt$base_plan_hash, plan$base_plan_hash) ||
      !identical(receipt$archive_plan_hash, plan$target_plan_hash)) {
    .update_abort("installed update receipt belongs to another transaction", "recovery")
  }
  receipt
}

.update_finish_installed <- function(root, plan, backup, state_lock = NULL) {
  status <- .update_generation_status(root, plan$target_plan_hash)
  if (!isTRUE(status$valid)) {
    .update_abort("installed target generation is not valid", "recovery", status = status)
  }
  receipt <- tryCatch(
    .update_existing_receipt(root, plan),
    error = function(e) {
      path <- file.path(root, "lna_update.json")
      if (file.exists(path)) {
        quarantine <- paste0(
          path, ".invalid-", format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC")
        )
        if (!file.rename(path, quarantine)) stop(e)
      }
      NULL
    }
  )
  if (is.null(receipt)) {
    receipt <- .update_write_receipt(root, plan, basename(backup))
  }
  if (!.update_journal_completed(root, plan$update_plan_hash)) {
    .update_journal_append(
      root, plan, "completed", list(previous_generation = basename(backup))
    )
  }
  receipt
}

.update_quarantine_path <- function(path, label) {
  if (!dir.exists(path)) return(NULL)
  destination <- paste0(
    path, ".", label, "-", format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC")
  )
  if (!file.rename(path, destination)) {
    .update_abort("could not quarantine an invalid update generation", "recovery")
  }
  destination
}

.update_recover_transaction <- function(
    root,
    plan,
    lock,
    prefer_rollback = FALSE) {
  state <- .update_read_state(root, required = FALSE)
  if (is.null(state)) {
    suffix <- substr(sub("^sha256:", "", plan$update_plan_hash), 1L, 16L)
    backup <- file.path(
      dirname(root),
      paste0(".", basename(root), ".lna-previous-r", plan$revision, "-", suffix)
    )
    root_target <- .update_generation_status(root, plan$target_plan_hash)
    preserved <- .update_generation_status(backup, plan$base_plan_hash)
    receipt <- tryCatch(
      .update_existing_receipt(root, plan),
      error = function(e) NULL
    )
    completed_receipt <- !is.null(receipt) &&
      .update_journal_completed(root, plan$update_plan_hash)
    if (isTRUE(root_target$valid) && isTRUE(completed_receipt)) {
      if (!isTRUE(preserved$valid)) {
        .update_abort(
          "installed target has no valid deterministic previous generation",
          "recovery", root = root_target, backup = preserved
        )
      }
      .update_audit_append(root, "confirmed-completed-transaction", list(
        backup = basename(backup)
      ), plan)
      return(list(
        completed = TRUE, recovered = TRUE, receipt = receipt,
        previous_generation = backup
      ))
    }
    if (isTRUE(root_target$valid) &&
        !identical(plan$base_plan_hash, plan$target_plan_hash)) {
      .update_abort(
        "target generation lacks completed transaction state or receipt",
        "recovery"
      )
    }
    if (!dir.exists(root)) {
      .update_abort("live root is missing and no transaction state exists", "recovery")
    }
    return(list(completed = FALSE, recovered = FALSE))
  }
  latest <- state$latest
  if (is.null(latest)) {
    root_base <- .update_generation_status(root, plan$base_plan_hash)
    if (!isTRUE(root_base$valid)) {
      .update_abort(
        "empty transaction state is only recoverable while the base root is valid",
        "recovery", root = root_base
      )
    }
    .update_audit_append(root, "discarded-empty-transaction-state", list(), plan)
    unlink(state$directory, recursive = TRUE, force = TRUE)
    return(list(completed = FALSE, recovered = TRUE))
  }
  expected <- c(
    update_plan_hash = plan$update_plan_hash,
    base_plan_hash = plan$base_plan_hash,
    target_plan_hash = plan$target_plan_hash
  )
  if (!all(vapply(names(expected), function(name) {
    identical(latest[[name]], expected[[name]])
  }, logical(1))) || !identical(as.integer(latest$revision), as.integer(plan$revision))) {
    .update_abort(
      "unfinished update transaction belongs to a different reviewed plan",
      "recovery", transaction = latest
    )
  }
  suffix <- substr(sub("^sha256:", "", plan$update_plan_hash), 1L, 16L)
  expected_stage <- paste0(".", basename(root), ".lna-stage-", suffix)
  expected_backup <- paste0(
    ".", basename(root), ".lna-previous-r", plan$revision, "-", suffix
  )
  if (!identical(latest$dataset_basename, basename(root)) ||
      !identical(latest$stage_basename, expected_stage) ||
      !identical(latest$backup_basename, expected_backup)) {
    .update_abort(
      "update transaction sibling paths are not deterministic for the reviewed plan",
      "recovery", transaction = latest
    )
  }
  parent <- dirname(root)
  stage <- file.path(parent, latest$stage_basename)
  backup <- file.path(parent, latest$backup_basename)
  root_base <- .update_generation_status(root, plan$base_plan_hash)
  root_target <- .update_generation_status(root, plan$target_plan_hash)
  staged <- .update_generation_status(stage, plan$target_plan_hash)
  preserved <- .update_generation_status(backup, plan$base_plan_hash)

  if (isTRUE(root_target$valid)) {
    if (isTRUE(staged$exists) || !isTRUE(preserved$valid)) {
      .update_abort(
        "ambiguous installed transaction: target root requires one valid backup and no stage",
        "recovery", root = root_target, stage = staged, backup = preserved
      )
    }
    if (isTRUE(prefer_rollback)) {
      if (!file.rename(root, stage) || !file.rename(backup, root)) {
        .update_abort("could not roll back the installed generation after an R error", "recovery")
      }
      .update_audit_append(root, "rolled-back-r-error", list(
        prior_phase = latest$phase, retained_stage = basename(stage)
      ), plan)
      unlink(state$directory, recursive = TRUE, force = TRUE)
      return(list(completed = FALSE, recovered = TRUE))
    }
    receipt <- .update_finish_installed(root, plan, backup)
    .update_audit_append(root, "completed-installed-transaction", list(
      prior_phase = latest$phase, backup = basename(backup)
    ), plan)
    unlink(state$directory, recursive = TRUE, force = TRUE)
    return(list(
      completed = TRUE, recovered = TRUE, receipt = receipt,
      previous_generation = backup
    ))
  }

  if (!isTRUE(root_base$exists)) {
    if (isTRUE(preserved$valid) && isTRUE(staged$valid) &&
        !isTRUE(prefer_rollback)) {
      .update_checkpoint(
        root, plan, stage, backup, "recovery-forward", lock
      )
      if (!file.rename(stage, root)) {
        .update_abort("could not resume installation of the reviewed generation", "recovery")
      }
      receipt <- .update_finish_installed(root, plan, backup)
      .update_audit_append(root, "resumed-reviewed-generation", list(
        prior_phase = latest$phase, backup = basename(backup)
      ), plan)
      unlink(state$directory, recursive = TRUE, force = TRUE)
      return(list(
        completed = TRUE, recovered = TRUE, receipt = receipt,
        previous_generation = backup
      ))
    }
    if (isTRUE(preserved$valid)) {
      if (isTRUE(root_base$exists) || isTRUE(root_target$exists)) {
        .update_quarantine_path(root, "invalid-root")
      }
      .update_checkpoint(
        root, plan, stage, backup, "recovery-rollback", lock
      )
      if (!file.rename(backup, root)) {
        .update_abort("could not restore the previous valid generation", "recovery")
      }
      .update_audit_append(root, "restored-previous-generation", list(
        prior_phase = latest$phase,
        invalid_stage = if (isTRUE(staged$exists)) basename(stage) else NULL
      ), plan)
      unlink(state$directory, recursive = TRUE, force = TRUE)
      return(list(completed = FALSE, recovered = TRUE))
    }
    .update_abort(
      "ambiguous transaction: no valid live root or previous generation is available",
      "recovery", root = root_base, stage = staged, backup = preserved
    )
  }

  if (!isTRUE(root_base$valid)) {
    if (isTRUE(preserved$valid)) {
      .update_quarantine_path(root, "invalid-root")
      .update_checkpoint(root, plan, stage, backup, "recovery-rollback", lock)
      if (!file.rename(backup, root)) {
        .update_abort("could not restore the previous valid generation", "recovery")
      }
      unlink(state$directory, recursive = TRUE, force = TRUE)
      return(list(completed = FALSE, recovered = TRUE))
    }
    .update_abort("live generation is invalid and no valid backup exists", "recovery")
  }
  if (isTRUE(preserved$exists)) {
    .update_abort(
      "ambiguous transaction: both live base and transaction backup exist",
      "recovery", root = root_base, backup = preserved
    )
  }
  if (isTRUE(prefer_rollback)) {
    .update_audit_append(root, "abandoned-staging-after-r-error", list(
      prior_phase = latest$phase,
      retained_stage = if (isTRUE(staged$exists)) basename(stage) else NULL
    ), plan)
    unlink(state$directory, recursive = TRUE, force = TRUE)
    return(list(completed = FALSE, recovered = TRUE))
  }
  .update_checkpoint(root, plan, stage, backup, "recovery-resume", lock)
  list(completed = FALSE, recovered = TRUE)
}

#' Execute a reviewed LNA dataset update preview
#'
#' A complete replacement generation is built beside the live archive and then
#' swapped into place. The former generation is retained as a recoverable
#' sibling directory. Mirror removals require `confirm_prune` to equal the
#' preview's `update_plan_hash`.
#'
#' @param plan A preview from [plan_lna_dataset_update()].
#' @param confirm_prune Exact preview hash authorizing planned mirror removals.
#' @param hooks Optional writer hooks, primarily for integration and fault tests.
#'   The update-only `after_update_checkpoint` hook receives each durable
#'   transaction phase and is not forwarded to the artifact writer.
#' @param owner_id Optional writer identity.
#' @param takeover_stale Whether an audited takeover of a non-live stale lock is
#'   permitted. A lock owned by a demonstrably live local process is never
#'   overridden.
#' @param stale_after_seconds Minimum lock age before an explicit stale takeover.
#' @return An `lna_dataset_update_result` with receipt, classifications, writer
#'   result, and recoverable previous-generation path.
#' @export
update_lna_dataset <- function(
    plan,
    confirm_prune = NULL,
    hooks = NULL,
    owner_id = NULL,
    takeover_stale = FALSE,
    stale_after_seconds = 3600) {
  plan <- .update_validate_plan(plan)
  root <- .update_root_path(plan$dataset_root, "dataset_root")
  if (!is.list(hooks) && !is.null(hooks)) {
    .update_abort("hooks must be a list", "hooks")
  }
  update_checkpoint_hook <- hooks$after_update_checkpoint %||%
    function(...) invisible(TRUE)
  if (!is.function(update_checkpoint_hook)) {
    .update_abort("after_update_checkpoint must be a function", "hooks")
  }
  writer_user_hooks <- hooks
  writer_user_hooks$after_update_checkpoint <- NULL
  removals <- vapply(plan$classifications, function(x) {
    identical(x$action, "remove")
  }, logical(1))
  if (any(removals) && !identical(confirm_prune, plan$update_plan_hash)) {
    .update_abort(
      "mirror removals require confirm_prune equal to the reviewed update_plan_hash",
      "confirm_prune"
    )
  }
  blocked <- vapply(plan$classifications, function(x) {
    identical(x$action, "blocked-corrupt") ||
      (identical(x$action, "retain") && identical(x$integrity_status, "corrupt"))
  }, logical(1))
  if (any(blocked)) {
    .update_abort(
      "the selected policy would retain corrupt artifacts; use refresh with available sources",
      "classifications",
      record_ids = vapply(plan$classifications[blocked], `[[`, character(1), "record_id")
    )
  }
  suffix <- substr(sub("^sha256:", "", plan$update_plan_hash), 1L, 16L)
  stage <- file.path(dirname(root), paste0(".", basename(root), ".lna-stage-", suffix))
  backup <- file.path(
    dirname(root),
    paste0(".", basename(root), ".lna-previous-r", plan$revision, "-", suffix)
  )
  lock <- .update_lock(
    root, plan,
    takeover_stale = takeover_stale,
    stale_after_seconds = stale_after_seconds
  )
  finished <- FALSE
  on.exit({
    if (!finished) {
      try(
        .update_recover_transaction(
          root, plan, lock, prefer_rollback = TRUE
        ),
        silent = TRUE
      )
    }
    try(.update_release_lock(lock), silent = TRUE)
  }, add = TRUE)
  recovery <- .update_recover_transaction(root, plan, lock)
  if (isTRUE(recovery$completed)) {
    finished <- TRUE
    .update_release_lock(lock)
    return(structure(
      list(
        plan = plan,
        receipt = recovery$receipt,
        classifications = as.data.frame(plan),
        writer_result = NULL,
        previous_generation = recovery$previous_generation,
        recovered = TRUE
      ),
      class = c("lna_dataset_update_result", "list")
    ))
  }
  .update_checkpoint(
    root, plan, stage, backup, "locked", lock, update_checkpoint_hook
  )
  current_plan <- read_lna_dataset_plan(file.path(root, ".lna_plan.json"))
  if (!identical(current_plan$plan_hash, plan$base_plan_hash)) {
    .update_abort("archive changed after the update preview was created", "base_plan_hash")
  }
  if (dir.exists(backup)) {
    .update_abort("recoverable previous-generation path already exists", "backup")
  }
  dir.create(stage, recursive = TRUE, showWarnings = FALSE)
  .update_checkpoint(
    root, plan, stage, backup, "staging", lock, update_checkpoint_hook
  )
  .update_journal_append(stage, plan, "started", list(base_plan_hash = plan$base_plan_hash))
  current_manifest <- read_lna_manifest(root)
  writer_hooks <- .update_hooks(plan, current_manifest, writer_user_hooks)
  writer_result <- tryCatch(
    write_lna_dataset(
      plan$target_plan,
      source_root = plan$source_root,
      overwrite = "error",
      hooks = writer_hooks,
      owner_id = owner_id,
      recover_truncated = "truncate",
      plan_mismatch = "error",
      physical_root = stage
    ),
    error = function(e) {
      .update_journal_append(stage, plan, "failed", list(message = conditionMessage(e)))
      stop(e)
    }
  )
  validate_lna_dataset_declaration(
    read_lna_dataset_declaration(stage), plan_hash = plan$target_plan_hash
  )
  staged_manifest <- read_lna_manifest(stage)
  if (!setequal(
      staged_manifest$data$record_id,
      vapply(plan$target_plan$artifacts, `[[`, character(1), "record_id"))) {
    .update_abort("staged generation manifest is incomplete", "manifest")
  }
  .update_journal_append(stage, plan, "ready-to-swap")
  .update_checkpoint(
    root, plan, stage, backup, "staged", lock, update_checkpoint_hook
  )
  .update_checkpoint(
    root, plan, stage, backup, "preserving-old", lock,
    update_checkpoint_hook
  )
  if (!file.rename(root, backup)) {
    .update_abort("could not preserve the current generation", "backup")
  }
  .update_checkpoint(
    root, plan, stage, backup, "old-preserved", lock,
    update_checkpoint_hook
  )
  .update_checkpoint(
    root, plan, stage, backup, "installing-new", lock,
    update_checkpoint_hook
  )
  if (!file.rename(stage, root)) {
    .update_abort("could not install the staged generation", "stage")
  }
  .update_checkpoint(
    root, plan, stage, backup, "new-installed", lock,
    update_checkpoint_hook
  )
  installed_status <- .update_generation_status(root, plan$target_plan_hash)
  if (!isTRUE(installed_status$valid)) {
    .update_abort("installed target generation failed validation", "stage")
  }
  .update_checkpoint(
    root, plan, stage, backup, "before-completion-receipt", lock,
    update_checkpoint_hook
  )
  receipt <- .update_existing_receipt(root, plan)
  if (is.null(receipt)) {
    receipt <- .update_write_receipt(root, plan, basename(backup))
  }
  .update_checkpoint(
    root, plan, stage, backup, "receipt-written", lock,
    update_checkpoint_hook
  )
  .update_journal_append(
    root, plan, "completed", list(previous_generation = basename(backup))
  )
  .update_checkpoint(
    root, plan, stage, backup, "completed", lock, update_checkpoint_hook
  )
  state <- .update_read_state(root, required = TRUE)
  unlink(state$directory, recursive = TRUE, force = TRUE)
  finished <- TRUE
  .update_release_lock(lock)
  structure(
    list(
      plan = plan,
      receipt = receipt,
      classifications = as.data.frame(plan),
      writer_result = writer_result,
      previous_generation = backup,
      recovered = isTRUE(recovery$recovered)
    ),
    class = c("lna_dataset_update_result", "list")
  )
}
