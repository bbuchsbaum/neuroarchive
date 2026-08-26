# Portable metadata, association, and privacy policy for LNA datasets.

.lna_metadata_policy_schema_version <- "1.0.0"
.lna_dataset_declaration_schema_version <- "1.0.0"
.lna_dataset_profile <- "NeuroarchiveBIDSExchange"
.lna_dataset_profile_version <- "1.0.0"
.lna_companion_roles <- c("events", "sessions", "scans", "confounds", "physio")

.metadata_role_family <- function(role) {
  aliases <- c(
    "events-json" = "events",
    "sessions-json" = "sessions",
    "scans-json" = "scans",
    "confounds-json" = "confounds",
    "physio-json" = "physio"
  )
  if (role %in% names(aliases)) unname(aliases[[role]]) else role
}

.metadata_abort <- function(message, field, record_id = NULL, ...) {
  abort_lna(
    message,
    ...,
    field = field,
    record_id = record_id,
    .subclass = "lna_error_metadata_policy",
    location = paste0("lna_metadata_policy:", field)
  )
}

.metadata_scalar_string <- function(x, field, allow_null = FALSE) {
  if (allow_null && is.null(x)) return(NULL)
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    .metadata_abort(paste0(field, " must be one non-empty string"), field)
  }
  enc2utf8(x)
}

.metadata_choice <- function(x, choices, field) {
  x <- .metadata_scalar_string(x, field)
  if (!x %in% choices) {
    .metadata_abort(
      paste0(field, " must be one of: ", paste(choices, collapse = ", ")),
      field
    )
  }
  x
}

.metadata_string_vector <- function(x, field, allow_empty = TRUE) {
  if (is.null(x) && allow_empty) return(character())
  if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
  if (!is.character(x) || anyNA(x) || any(!nzchar(x)) || anyDuplicated(x) ||
      (!allow_empty && !length(x))) {
    .metadata_abort(
      paste0(field, " must contain unique non-empty strings"), field
    )
  }
  enc2utf8(x)
}

.metadata_participants_policy <- function(x) {
  if (is.character(x) && length(x) == 1L) x <- list(mode = x)
  if (!is.list(x) || is.data.frame(x) || is.null(names(x))) {
    .metadata_abort("participants must be a mode or named list", "participants")
  }
  allowed <- c(
    "mode", "source_path", "source_checksum", "columns",
    "selected_participant_ids"
  )
  unknown <- setdiff(names(x), allowed)
  if (length(unknown)) {
    .metadata_abort(
      paste0("participants has unknown field(s): ", paste(unknown, collapse = ", ")),
      "participants"
    )
  }
  mode <- .metadata_choice(
    x$mode %||% "exclude", c("exclude", "include-subset"),
    "participants.mode"
  )
  if (identical(mode, "exclude")) {
    return(list(
      mode = mode, source_path = NULL, source_checksum = NULL,
      columns = character(), selected_participant_ids = character()
    ))
  }
  source_path <- tryCatch(
    .manifest_validate_path(x$source_path, "participants.source_path"),
    lna_error_manifest = function(e) {
      .metadata_abort(conditionMessage(e), "participants.source_path")
    }
  )
  checksum <- tryCatch(
    .source_record_checksum(x$source_checksum, "participants.source_checksum"),
    lna_error_source_record = function(e) {
      .metadata_abort(conditionMessage(e), "participants.source_checksum")
    }
  )
  columns <- .metadata_string_vector(
    x$columns %||% "participant_id", "participants.columns", allow_empty = FALSE
  )
  if (!"participant_id" %in% columns) {
    .metadata_abort(
      "participants.columns must include participant_id", "participants.columns"
    )
  }
  list(
    mode = mode,
    source_path = source_path,
    source_checksum = checksum,
    columns = c("participant_id", setdiff(columns, "participant_id")),
    selected_participant_ids = .metadata_string_vector(
      x$selected_participant_ids %||% character(),
      "participants.selected_participant_ids"
    )
  )
}

#' Define portable BIDS metadata and privacy policy
#'
#' Neuroarchive receives already-resolved metadata from an ingestion adapter.
#' It records the sidecars used to resolve each field, copies declared
#' companions, and applies an explicit privacy policy. It does not rediscover
#' BIDS inheritance or interpret event conditions.
#'
#' @param portable If `TRUE`, require field-level metadata provenance. The
#'   default is permissive so pure planning can represent incomplete adapter
#'   output; [write_lna_dataset()] uses the portable policy by default.
#' @param metadata_provenance Whether every resolved metadata field must name a
#'   `metadata-source` association.
#' @param events,sessions,scans,confounds,physio Inclusion modes: `"exclude"`,
#'   `"copy-if-present"`, or `"require"`.
#' @param participants `"exclude"` (the privacy-safe default), or a named list
#'   with `mode = "include-subset"`, checksummed `source_path`, and selected
#'   `columns`. Rows are restricted to subjects represented by the plan.
#' @return A validated `lna_metadata_policy`.
#' @export
lna_metadata_policy <- function(
    portable = FALSE,
    metadata_provenance = NULL,
    events = "copy-if-present",
    sessions = "exclude",
    scans = "exclude",
    confounds = "copy-if-present",
    physio = "copy-if-present",
    participants = "exclude") {
  if (is.null(metadata_provenance)) {
    metadata_provenance <- if (isTRUE(portable)) "required" else "allow-unrecorded"
  }
  validate_lna_metadata_policy(list(
    schema_version = .lna_metadata_policy_schema_version,
    inheritance = "flatten-effective",
    metadata_provenance = metadata_provenance,
    companions = list(
      events = events,
      sessions = sessions,
      scans = scans,
      confounds = confounds,
      physio = physio
    ),
    participants = participants,
    task_semantics = list(
      task_entity = "acquisition-label",
      within_run_conditions = "opaque-events-content",
      events_normalization = "none-byte-preserved"
    )
  ))
}

#' @rdname lna_metadata_policy
#' @param policy A metadata policy or decoded JSON object.
#' @export
validate_lna_metadata_policy <- function(policy) {
  if (inherits(policy, "lna_metadata_policy")) policy <- unclass(policy)
  if (!is.list(policy) || is.data.frame(policy) || is.null(names(policy))) {
    .metadata_abort("policy must be a named object", "policy")
  }
  required <- c(
    "schema_version", "inheritance", "metadata_provenance", "companions",
    "participants", "task_semantics"
  )
  if (!setequal(names(policy), required)) {
    .metadata_abort(
      paste0("policy must contain exactly: ", paste(required, collapse = ", ")),
      "policy"
    )
  }
  version <- .metadata_scalar_string(policy$schema_version, "schema_version")
  if (!grepl("^1[.][0-9]+[.][0-9]+$", version)) {
    .metadata_abort("unsupported metadata policy schema version", "schema_version")
  }
  inheritance <- .metadata_choice(
    policy$inheritance, "flatten-effective", "inheritance"
  )
  provenance <- .metadata_choice(
    policy$metadata_provenance, c("required", "allow-unrecorded"),
    "metadata_provenance"
  )
  companions <- policy$companions
  if (!is.list(companions) || !setequal(names(companions), .lna_companion_roles)) {
    .metadata_abort(
      paste0("companions must contain exactly: ",
        paste(.lna_companion_roles, collapse = ", ")),
      "companions"
    )
  }
  companions <- stats::setNames(lapply(.lna_companion_roles, function(role) {
    .metadata_choice(
      companions[[role]], c("exclude", "copy-if-present", "require"),
      paste0("companions.", role)
    )
  }), .lna_companion_roles)
  expected_task <- list(
    task_entity = "acquisition-label",
    within_run_conditions = "opaque-events-content",
    events_normalization = "none-byte-preserved"
  )
  if (!identical(policy$task_semantics, expected_task)) {
    .metadata_abort("task_semantics is incompatible", "task_semantics")
  }
  structure(list(
    schema_version = version,
    inheritance = inheritance,
    metadata_provenance = provenance,
    companions = companions,
    participants = .metadata_participants_policy(policy$participants),
    task_semantics = expected_task
  ), class = c("lna_metadata_policy", "list"))
}

.metadata_provenance_source <- function(companion, record_id) {
  field <- paste0("artifacts.", record_id, ".metadata_provenance")
  metadata <- companion$metadata
  allowed <- c("inheritance_level", "precedence", "fields", "media_type")
  unknown <- setdiff(names(metadata), allowed)
  required <- c("inheritance_level", "precedence", "fields")
  if (length(unknown) || length(setdiff(required, names(metadata)))) {
    .metadata_abort(
      paste0("metadata-source association must contain ",
        paste(required, collapse = ", "), " and optional media_type"),
      field, record_id
    )
  }
  level <- .metadata_choice(
    metadata$inheritance_level,
    c("dataset", "subject", "session", "datatype", "task", "acquisition", "run"),
    paste0(field, ".inheritance_level")
  )
  precedence <- if (is.list(metadata$precedence)) {
    unlist(metadata$precedence, recursive = TRUE, use.names = FALSE)
  } else metadata$precedence
  integer_precedence <- suppressWarnings(as.integer(precedence))
  if (!is.numeric(precedence) || length(precedence) != 1L ||
      !is.finite(precedence) || precedence != integer_precedence ||
      integer_precedence < 0L) {
    .metadata_abort(
      "metadata-source precedence must be one non-negative integer",
      paste0(field, ".precedence"), record_id
    )
  }
  fields <- .metadata_string_vector(
    metadata$fields, paste0(field, ".fields"), allow_empty = FALSE
  )
  list(
    association_id = companion$association_id,
    source_path = companion$source_path,
    source_uri = companion$source_uri,
    source_checksum = companion$source_checksum,
    inheritance_level = level,
    precedence = integer_precedence,
    fields = sort(fields),
    media_type = .metadata_scalar_string(
      metadata$media_type %||% "application/json",
      paste0(field, ".media_type")
    )
  )
}

.metadata_entry_closure <- function(entry, policy) {
  record_id <- entry$record_id
  companions <- entry$expected_companions
  roles <- if (length(companions)) {
    vapply(companions, `[[`, character(1), "role")
  } else character()
  role_families <- vapply(roles, .metadata_role_family, character(1))
  keep <- rep(TRUE, length(companions))
  for (role in .lna_companion_roles) {
    mode <- policy$companions[[role]]
    present <- which(role_families == role)
    if (identical(mode, "require") && !length(present)) {
      .metadata_abort(
        paste0("required ", role, " companion is missing"),
        paste0("companions.", role), record_id
      )
    }
    if (identical(mode, "exclude") && length(present)) keep[present] <- FALSE
  }
  participant_items <- which(roles == "participants")
  if (length(participant_items)) keep[participant_items] <- FALSE
  companions <- companions[keep]

  provenance_items <- companions[vapply(companions, function(item) {
    identical(item$role, "metadata-source")
  }, logical(1))]
  sources <- lapply(provenance_items, .metadata_provenance_source, record_id = record_id)
  if (length(sources)) {
    source_order <- order(
      vapply(sources, `[[`, integer(1), "precedence"),
      vapply(sources, `[[`, character(1), "association_id")
    )
    sources <- sources[source_order]
  }
  metadata_fields <- names(entry$resolved_metadata)
  field_sources <- stats::setNames(vector("list", length(metadata_fields)), metadata_fields)
  for (field in metadata_fields) {
    candidates <- sources[vapply(sources, function(source) {
      field %in% source$fields
    }, logical(1))]
    if (!length(candidates)) {
      if (identical(policy$metadata_provenance, "required")) {
        .metadata_abort(
          paste0("resolved metadata field lacks provenance: ", field),
          "metadata_provenance", record_id
        )
      }
      field_sources[[field]] <- list(status = "unrecorded", association_id = NULL)
      next
    }
    precedence <- vapply(candidates, `[[`, integer(1), "precedence")
    winners <- candidates[precedence == max(precedence)]
    if (length(winners) != 1L) {
      .metadata_abort(
        paste0("resolved metadata field has ambiguous highest-precedence sources: ", field),
        "metadata_provenance", record_id
      )
    }
    winner <- winners[[1L]]
    field_sources[[field]] <- list(
      status = "recorded",
      association_id = winner$association_id,
      source_checksum = winner$source_checksum,
      precedence = winner$precedence
    )
  }
  covered <- unique(unlist(lapply(sources, `[[`, "fields"), use.names = FALSE))
  unknown_fields <- setdiff(covered, metadata_fields)
  if (length(unknown_fields)) {
    .metadata_abort(
      paste0("metadata provenance names unresolved field(s): ",
        paste(unknown_fields, collapse = ", ")),
      "metadata_provenance", record_id
    )
  }
  status <- if (!length(metadata_fields) || all(vapply(field_sources, function(x) {
    identical(x$status, "recorded")
  }, logical(1)))) "complete" else "unrecorded-fields"
  list(
    companions = companions,
    closure = list(
      inheritance = policy$inheritance,
      provenance_requirement = policy$metadata_provenance,
      status = status,
      sources = sources,
      field_sources = field_sources
    )
  )
}

.metadata_apply_plan <- function(entries, policy) {
  policy <- validate_lna_metadata_policy(policy)
  entries <- lapply(entries, function(entry) {
    closure <- .metadata_entry_closure(entry, policy)
    entry$expected_companions <- closure$companions
    entry$metadata_provenance <- closure$closure
    entry
  })
  participant_ids <- sort(unique(vapply(entries, function(entry) {
    paste0("sub-", entry$entities$sub)
  }, character(1))))
  policy$participants$selected_participant_ids <- if (
    identical(policy$participants$mode, "include-subset")
  ) participant_ids else character()
  list(entries = entries, policy = unclass(policy))
}

.metadata_json_identical <- function(x, y) {
  identical(
    jsonlite::toJSON(
      x, auto_unbox = TRUE, null = "null", digits = NA,
      matrix = "rowmajor", pretty = FALSE
    ),
    jsonlite::toJSON(
      y, auto_unbox = TRUE, null = "null", digits = NA,
      matrix = "rowmajor", pretty = FALSE
    )
  )
}

#' Validate metadata closure in an LNA dataset plan
#'
#' @param plan A dataset plan.
#' @return The validated plan, invisibly.
#' @export
validate_lna_metadata_closure <- function(plan) {
  policy <- validate_lna_metadata_policy(plan$metadata_policy)
  applied <- .metadata_apply_plan(plan$artifacts, policy)
  comparable <- lapply(plan$artifacts, function(entry) {
    entry[c("expected_companions", "metadata_provenance")]
  })
  expected <- lapply(applied$entries, function(entry) {
    entry[c("expected_companions", "metadata_provenance")]
  })
  if (!.metadata_json_identical(comparable, expected) ||
      !.metadata_json_identical(plan$metadata_policy, applied$policy)) {
    .metadata_abort(
      "plan metadata closure does not match its declared policy", "plan"
    )
  }
  invisible(plan)
}

#' Create a machine-readable Neuroarchive dataset declaration
#'
#' @param plan A validated dataset plan.
#' @param participants_receipt Optional receipt from privacy-safe participant
#'   subsetting.
#' @return A JSON-safe declaration object.
#' @export
lna_dataset_declaration <- function(plan, participants_receipt = NULL) {
  plan <- validate_lna_dataset_plan(plan)
  declaration <- list(
    schema_version = .lna_dataset_declaration_schema_version,
    profile = .lna_dataset_profile,
    profile_version = .lna_dataset_profile_version,
    plan_hash = plan$plan_hash,
    manifest = list(
      table = "lna_manifest.tsv",
      sidecar = "lna_manifest.json",
      schema_version = .lna_manifest_schema_version,
      authority = "manifest-plus-declaration"
    ),
    extension_registry = list(list(
      extension = ".lna.h5",
      media_type = "application/vnd.neuroarchive.lna+hdf5",
      row_semantics = "one-source-image-one-artifact",
      bids_recognized_extension = FALSE
    )),
    bids_envelope = list(
      mode = "BIDS-shaped-envelope",
      divergence = paste0(
        ".lna.h5 is a Neuroarchive extension and is not a recognized BIDS ",
        "imaging extension"
      ),
      entity_comparison = "exact",
      inheritance_owner = "source-record-adapter"
    ),
    metadata_policy = plan$metadata_policy,
    fit_resources = plan$execution$fit_resources,
    participants_receipt = participants_receipt
  )
  validate_lna_dataset_declaration(declaration, plan_hash = plan$plan_hash)
}

#' @rdname lna_dataset_declaration
#' @param declaration A declaration object.
#' @param plan_hash Optional expected plan hash.
#' @export
validate_lna_dataset_declaration <- function(declaration, plan_hash = NULL) {
  if (!is.list(declaration) || is.data.frame(declaration)) {
    .metadata_abort("declaration must be an object", "declaration")
  }
  required <- c(
    "schema_version", "profile", "profile_version", "plan_hash", "manifest",
    "extension_registry", "bids_envelope", "metadata_policy", "fit_resources",
    "participants_receipt"
  )
  if (!setequal(names(declaration), required)) {
    .metadata_abort(
      paste0("declaration must contain exactly: ", paste(required, collapse = ", ")),
      "declaration"
    )
  }
  if (!identical(declaration$schema_version, .lna_dataset_declaration_schema_version) ||
      !identical(declaration$profile, .lna_dataset_profile) ||
      !identical(declaration$profile_version, .lna_dataset_profile_version)) {
    .metadata_abort("unsupported dataset declaration profile", "declaration")
  }
  declared_hash <- tryCatch(
    .source_record_checksum(declaration$plan_hash, "declaration.plan_hash"),
    lna_error_source_record = function(e) {
      .metadata_abort(conditionMessage(e), "declaration.plan_hash")
    }
  )
  if (!is.null(plan_hash) && !identical(declared_hash, plan_hash)) {
    .metadata_abort("declaration belongs to a different plan", "declaration.plan_hash")
  }
  manifest <- declaration$manifest
  manifest_schema_supported <- tryCatch(
    .manifest_schema_version_supported(manifest$schema_version),
    error = function(e) FALSE
  )
  if (!is.list(manifest) ||
      !identical(manifest$table, "lna_manifest.tsv") ||
      !identical(manifest$sidecar, "lna_manifest.json") ||
      !isTRUE(manifest_schema_supported) ||
      !identical(manifest$authority, "manifest-plus-declaration")) {
    .metadata_abort("declaration manifest contract is incompatible", "declaration.manifest")
  }
  registry <- declaration$extension_registry
  if (!is.list(registry) || length(registry) != 1L ||
      !identical(registry[[1L]]$extension, ".lna.h5") ||
      !identical(registry[[1L]]$media_type,
        "application/vnd.neuroarchive.lna+hdf5") ||
      !identical(registry[[1L]]$bids_recognized_extension, FALSE)) {
    .metadata_abort("declaration extension registry is incompatible", "extension_registry")
  }
  envelope <- declaration$bids_envelope
  if (!is.list(envelope) ||
      !identical(envelope$mode, "BIDS-shaped-envelope") ||
      !identical(envelope$entity_comparison, "exact") ||
      !identical(envelope$inheritance_owner, "source-record-adapter")) {
    .metadata_abort("declaration BIDS envelope is incompatible", "bids_envelope")
  }
  declaration$metadata_policy <- unclass(validate_lna_metadata_policy(
    declaration$metadata_policy
  ))
  resources <- declaration$fit_resources
  if (!is.list(resources) || is.data.frame(resources)) {
    .metadata_abort("declaration fit_resources must be a list", "fit_resources")
  }
  if (length(resources)) {
    resource_record_ids <- unique(unlist(lapply(resources, function(resource) {
      c(
        unlist(resource$input_record_ids, recursive = TRUE, use.names = FALSE),
        unlist(resource$required_by_record_ids, recursive = TRUE, use.names = FALSE)
      )
    }), recursive = TRUE, use.names = FALSE))
    resources <- .execution_fit_resources(resources, resource_record_ids)
    resources <- lapply(resources, function(resource) {
      if (identical(resource$resource_type, "shared-basis")) {
        resource$parameters <- .dataset_shared_basis_parameters(resource)
      }
      resource
    })
  }
  declaration$fit_resources <- resources
  receipt <- declaration$participants_receipt
  if (!is.null(receipt)) {
    receipt_fields <- c(
      "mode", "source_path", "source_checksum", "selected_participant_ids",
      "selected_columns", "destination", "output_checksum", "output_bytes"
    )
    if (!is.list(receipt) || !setequal(names(receipt), receipt_fields) ||
        !identical(receipt$mode, "include-subset")) {
      .metadata_abort(
        "participants_receipt is incompatible", "participants_receipt"
      )
    }
    bytes <- if (is.list(receipt$output_bytes)) {
      unlist(receipt$output_bytes, recursive = TRUE, use.names = FALSE)
    } else receipt$output_bytes
    if (!is.numeric(bytes) || length(bytes) != 1L || !is.finite(bytes) ||
        bytes < 0 || bytes != floor(bytes)) {
      .metadata_abort(
        "participants_receipt.output_bytes must be a non-negative integer",
        "participants_receipt.output_bytes"
      )
    }
    declaration$participants_receipt <- list(
      mode = "include-subset",
      source_path = .manifest_validate_path(
        receipt$source_path, "participants_receipt.source_path"
      ),
      source_checksum = .source_record_checksum(
        receipt$source_checksum, "participants_receipt.source_checksum"
      ),
      selected_participant_ids = .metadata_string_vector(
        receipt$selected_participant_ids,
        "participants_receipt.selected_participant_ids", allow_empty = FALSE
      ),
      selected_columns = .metadata_string_vector(
        receipt$selected_columns,
        "participants_receipt.selected_columns", allow_empty = FALSE
      ),
      destination = .manifest_validate_path(
        receipt$destination, "participants_receipt.destination"
      ),
      output_checksum = .source_record_checksum(
        receipt$output_checksum, "participants_receipt.output_checksum"
      ),
      output_bytes = as.numeric(bytes)
    )
  }
  .source_record_json_value(declaration, "declaration")
  declaration
}

#' @rdname lna_dataset_declaration
#' @param path Archive root or `lna_dataset.json` path.
#' @param overwrite Whether to replace an existing declaration.
#' @export
write_lna_dataset_declaration <- function(
    declaration,
    path,
    overwrite = FALSE) {
  declaration <- validate_lna_dataset_declaration(declaration)
  path <- .metadata_scalar_string(path, "path")
  if (!identical(basename(path), "lna_dataset.json")) {
    path <- file.path(path, "lna_dataset.json")
  }
  if (file.exists(path) && !isTRUE(overwrite)) {
    .metadata_abort("dataset declaration already exists", "path")
  }
  if (!dir.exists(dirname(path)) &&
      !dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)) {
    .metadata_abort("could not create declaration directory", "path")
  }
  jsonlite::write_json(
    declaration, path, auto_unbox = TRUE, null = "null", digits = NA,
    pretty = TRUE
  )
  invisible(path)
}

#' @rdname lna_dataset_declaration
#' @export
read_lna_dataset_declaration <- function(path) {
  path <- .metadata_scalar_string(path, "path")
  if (!identical(basename(path), "lna_dataset.json")) {
    path <- file.path(path, "lna_dataset.json")
  }
  if (!file.exists(path)) {
    .metadata_abort("dataset declaration does not exist", "path")
  }
  decoded <- tryCatch(
    jsonlite::read_json(path, simplifyVector = FALSE),
    error = function(e) .metadata_abort("dataset declaration is invalid JSON", "path")
  )
  validate_lna_dataset_declaration(decoded)
}

.metadata_write_participants <- function(plan, source_root, output_root, overwrite) {
  policy <- plan$metadata_policy$participants
  if (identical(policy$mode, "exclude")) return(NULL)
  source <- .writer_full_path(source_root, policy$source_path)
  if (!file.exists(source) ||
      !identical(.journal_file_checksum(source), policy$source_checksum)) {
    .metadata_abort(
      "participants source is missing or has changed", "participants.source_path"
    )
  }
  table <- tryCatch(
    utils::read.delim(
      source, colClasses = "character", check.names = FALSE,
      stringsAsFactors = FALSE, na.strings = character()
    ),
    error = function(e) .metadata_abort(
      paste0("could not read participants table: ", conditionMessage(e)),
      "participants.source_path"
    )
  )
  missing_columns <- setdiff(policy$columns, names(table))
  if (length(missing_columns)) {
    .metadata_abort(
      paste0("participants table is missing selected column(s): ",
        paste(missing_columns, collapse = ", ")),
      "participants.columns"
    )
  }
  selected <- policy$selected_participant_ids
  if (anyDuplicated(table$participant_id)) {
    .metadata_abort("participants participant_id values must be unique", "participants")
  }
  missing_rows <- setdiff(selected, table$participant_id)
  if (length(missing_rows)) {
    .metadata_abort(
      paste0("participants table lacks selected subject(s): ",
        paste(missing_rows, collapse = ", ")),
      "participants"
    )
  }
  subset <- table[match(selected, table$participant_id), policy$columns, drop = FALSE]
  destination <- file.path(output_root, "participants.tsv")
  temp <- tempfile(".participants-", tmpdir = output_root)
  on.exit(if (file.exists(temp)) unlink(temp, force = TRUE), add = TRUE)
  utils::write.table(
    subset, temp, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = "n/a", fileEncoding = "UTF-8"
  )
  receipt <- function(path) list(
    mode = "include-subset",
    source_path = policy$source_path,
    source_checksum = policy$source_checksum,
    selected_participant_ids = selected,
    selected_columns = policy$columns,
    destination = "participants.tsv",
    output_checksum = .journal_file_checksum(path),
    output_bytes = unname(file.info(path)$size)
  )
  if (file.exists(destination) &&
      identical(.journal_file_checksum(destination), .journal_file_checksum(temp))) {
    unlink(temp, force = TRUE)
    return(receipt(destination))
  }
  if (file.exists(destination) && !identical(overwrite, "rebuild")) {
    .metadata_abort("participants.tsv already exists", "participants")
  }
  backup <- NULL
  if (file.exists(destination)) {
    backup <- tempfile(".participants-backup-", tmpdir = output_root)
    if (!file.rename(destination, backup)) {
      .metadata_abort("could not preserve existing participants.tsv", "participants")
    }
  }
  installed <- FALSE
  on.exit({
    if (!installed && !is.null(backup) && file.exists(backup)) {
      if (file.exists(destination)) unlink(destination, force = TRUE)
      file.rename(backup, destination)
    }
  }, add = TRUE)
  if (!file.rename(temp, destination)) {
    .metadata_abort("could not install participants.tsv", "participants")
  }
  installed <- TRUE
  if (!is.null(backup) && file.exists(backup)) unlink(backup, force = TRUE)
  receipt(destination)
}
