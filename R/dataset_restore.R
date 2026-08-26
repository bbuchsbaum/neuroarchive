# Closed and explicitly partial restoration of LNA exchange datasets to BIDS.

.lna_restoration_schema_version <- "1.0.0"
.lna_restoration_bids_version <- "1.11.1"

.restore_abort <- function(message, field, ...) {
  abort_lna(
    message,
    ...,
    field = field,
    .subclass = "lna_error_restoration",
    location = paste0("restore_bids:", field)
  )
}

.restore_hooks <- function(hooks) {
  defaults <- list(
    reconstruct = function(lna_path, nifti_path, row, dataset) {
      contract <- validate_lna_error_contract(
        row$codec_profile_json[[1L]]$error_contract
      )
      vec <- as_neurovec(
        lna_path,
        output_dtype = if (identical(contract$output_dtype, "float64")) {
          "float64"
        } else "float32"
      )
      tr <- row$repetition_time_seconds[[1L]]
      if (!is.na(tr)) {
        header <- methods::slot(vec, "header")
        header$xyzt_units <- 10L
        header$pixdim <- c(1, neuroim2::spacing(vec)[1:3], tr, 0, 0, 0)
        methods::slot(vec, "header") <- header
      }
      nifti_dtype <- switch(
        contract$output_dtype,
        float64 = "DOUBLE",
        float32 = "FLOAT",
        integer32 = "INT"
      )
      neuroim2::write_vec(vec, nifti_path, data_type = nifti_dtype)
      invisible(nifti_path)
    },
    validate_nifti = function(path, row) {
      neuroim2::read_vec(path)
      invisible(TRUE)
    },
    copy_companion = function(source, destination) {
      file.copy(source, destination, overwrite = FALSE, copy.mode = TRUE)
    }
  )
  if (is.null(hooks)) return(defaults)
  if (!is.list(hooks) || is.data.frame(hooks) || is.null(names(hooks)) ||
      any(!nzchar(names(hooks))) || anyDuplicated(names(hooks))) {
    .restore_abort("hooks must be a uniquely named list", "hooks")
  }
  unknown <- setdiff(names(hooks), names(defaults))
  if (length(unknown)) {
    .restore_abort(
      paste0("unknown restoration hook(s): ", paste(unknown, collapse = ", ")),
      "hooks"
    )
  }
  for (name in names(hooks)) {
    if (!is.function(hooks[[name]])) {
      .restore_abort(paste0("hooks.", name, " must be a function"), "hooks")
    }
    defaults[[name]] <- hooks[[name]]
  }
  defaults
}

.restore_selection <- function(dataset, selection, match) {
  if (!is.list(selection) || is.data.frame(selection) ||
      (length(selection) && (is.null(names(selection)) ||
        any(!nzchar(names(selection))) || anyDuplicated(names(selection))))) {
    .restore_abort("selection must be a uniquely named list", "selection")
  }
  args <- c(list(dataset = dataset), selection, list(match = match, view = "all"))
  selected <- do.call(query_lna, args)
  if (!nrow(selected)) {
    .restore_abort("selection matched no records", "selection")
  }
  selected
}

.restore_associations <- function(row) {
  associations <- row$associations_json[[1L]]
  if (is.null(associations) || !length(associations)) return(list())
  associations
}

.restore_transitive_records <- function(data, requested_ids, closure) {
  available <- data$record_id
  selected <- requested_ids
  missing <- list()
  repeat {
    added <- character()
    rows <- data[data$record_id %in% selected, , drop = FALSE]
    for (i in seq_len(nrow(rows))) {
      associations <- .restore_associations(rows[i, , drop = FALSE])
      for (association in associations) {
        target <- association$target_record_id
        if (is.null(target)) next
        if (!target %in% available) {
          missing[[length(missing) + 1L]] <- list(
            kind = "record-association",
            record_id = rows$record_id[[i]],
            association_id = association$association_id,
            requirement = target,
            reason = "target record is absent from the archive"
          )
        } else if (!target %in% selected) {
          if (identical(closure, "closed")) {
            added <- c(added, target)
          } else {
            missing[[length(missing) + 1L]] <- list(
              kind = "record-association",
              record_id = rows$record_id[[i]],
              association_id = association$association_id,
              requirement = target,
              reason = "target record was outside the partial selection"
            )
          }
        }
      }
    }
    added <- setdiff(unique(added), selected)
    if (!length(added)) break
    selected <- c(selected, added)
  }
  list(record_ids = unique(selected), missing = missing)
}

.restore_nifti_path <- function(lna_path) {
  if (!endsWith(lna_path, ".lna.h5")) {
    .restore_abort("manifest LNA path lacks .lna.h5 extension", "lna_path")
  }
  sub("[.]lna[.]h5$", ".nii.gz", lna_path)
}

.restore_sidecar_path <- function(nifti_path) {
  sub("[.]nii[.]gz$", ".json", nifti_path)
}

.restore_archive_path <- function(dataset, relative) {
  .writer_full_path(dataset$root, relative)
}

.restore_missing_dependency <- function(
    kind,
    record_id,
    requirement,
    reason,
    association_id = NULL) {
  list(
    kind = kind,
    record_id = record_id,
    association_id = association_id,
    requirement = requirement,
    reason = reason
  )
}

.restore_preflight <- function(dataset, selected, closure) {
  missing <- list()
  restorable <- rep(TRUE, nrow(selected))
  for (i in seq_len(nrow(selected))) {
    row <- selected[i, , drop = FALSE]
    record_id <- row$record_id[[1L]]
    lna_path <- .restore_archive_path(dataset, row$lna_path[[1L]])
    complete <- identical(row$status[[1L]], "complete") &&
      !is.na(row$lna_checksum[[1L]]) && !is.na(row$lna_bytes[[1L]])
    valid <- complete && file.exists(lna_path) && identical(
      .journal_file_checksum(lna_path), row$lna_checksum[[1L]]
    )
    if (!valid) {
      restorable[[i]] <- FALSE
      missing[[length(missing) + 1L]] <- .restore_missing_dependency(
        "lna-artifact", record_id, row$lna_path[[1L]],
        "artifact is absent, incomplete, or checksum-invalid"
      )
    }
    provenance <- row$metadata_provenance_json[[1L]]
    if (identical(closure, "closed") &&
        !identical(provenance$status, "complete")) {
      missing[[length(missing) + 1L]] <- .restore_missing_dependency(
        "metadata-provenance", record_id, "complete field provenance",
        "effective metadata provenance is incomplete"
      )
    }
    if (identical(row$datatype[[1L]], "func") &&
        identical(row$suffix[[1L]], "bold")) {
      tr <- row$repetition_time_seconds[[1L]]
      if (is.na(tr) || !is.finite(tr) || tr <= 0) {
        missing[[length(missing) + 1L]] <- .restore_missing_dependency(
          "metadata", record_id, "positive RepetitionTime",
          "BOLD restoration requires RepetitionTime in seconds"
        )
      }
    }
    for (association in .restore_associations(row)) {
      relative <- association$destination_relative
      if (is.null(relative)) next
      source <- .restore_archive_path(dataset, relative)
      valid_companion <- file.exists(source) && identical(
        .journal_file_checksum(source), association$source_checksum
      )
      if (!valid_companion) {
        missing[[length(missing) + 1L]] <- .restore_missing_dependency(
          "companion", record_id, relative,
          "declared companion is absent or checksum-invalid",
          association$association_id
        )
      }
    }
  }
  if (identical(closure, "closed") && length(missing)) {
    .restore_abort(
      paste0("closed restoration has ", length(missing), " missing requirement(s)"),
      "closure", missing_dependencies = missing
    )
  }
  list(
    selected = selected[restorable, , drop = FALSE],
    missing = missing
  )
}

.restore_companion_plan <- function(dataset, selected, closure, missing) {
  skip_roles <- c("metadata-source", "participants")
  items <- list()
  for (i in seq_len(nrow(selected))) {
    row <- selected[i, , drop = FALSE]
    for (association in .restore_associations(row)) {
      relative <- association$destination_relative
      if (is.null(relative) || association$role %in% skip_roles) next
      source <- .restore_archive_path(dataset, relative)
      if (!file.exists(source) || !identical(
        .journal_file_checksum(source), association$source_checksum
      )) next
      items[[length(items) + 1L]] <- list(
        record_id = row$record_id[[1L]],
        association_id = association$association_id,
        role = association$role,
        source = source,
        destination_relative = relative,
        checksum = association$source_checksum
      )
    }
  }
  if (!length(items)) return(list())
  destinations <- vapply(items, `[[`, character(1), "destination_relative")
  groups <- split(seq_along(items), destinations)
  out <- lapply(names(groups), function(destination) {
    group <- items[groups[[destination]]]
    checksums <- unique(vapply(group, `[[`, character(1), "checksum"))
    if (length(checksums) != 1L) {
      .restore_abort(
        paste0("selected companions conflict at ", destination), "closure"
      )
    }
    group[[1L]]
  })
  out[order(vapply(out, `[[`, character(1), "destination_relative"))]
}

.restore_output_inventory <- function(selected, companions, participants) {
  nifti <- vapply(selected$lna_path, .restore_nifti_path, character(1))
  sidecars <- vapply(nifti, .restore_sidecar_path, character(1))
  companion_paths <- if (length(companions)) {
    vapply(companions, `[[`, character(1), "destination_relative")
  } else character()
  participant_path <- if (is.null(participants)) character() else "participants.tsv"
  generated <- c(
    nifti, sidecars, companion_paths, participant_path,
    "dataset_description.json", "README", ".bidsignore", "lna_restoration.json"
  )
  pipelines <- sort(unique(selected$pipeline[
    selected$scope == "derivative" & !is.na(selected$pipeline)
  ]))
  generated <- c(
    generated,
    file.path("derivatives", pipelines, "dataset_description.json")
  )
  duplicates <- unique(generated[duplicated(generated)])
  if (length(duplicates)) {
    .restore_abort(
      paste0("restoration output collision(s): ", paste(duplicates, collapse = ", ")),
      "output"
    )
  }
  list(
    nifti = nifti,
    sidecars = sidecars,
    pipelines = pipelines
  )
}

.restore_package_version <- function() {
  tryCatch(as.character(utils::packageVersion("neuroarchive")),
    error = function(e) "development")
}

.restore_dataset_description <- function(type = "raw", pipeline = NULL, dataset = NULL) {
  if (identical(type, "raw")) {
    return(list(
      Name = "Neuroarchive restored dataset",
      BIDSVersion = .lna_restoration_bids_version,
      DatasetType = "raw",
      Authors = list("Neuroarchive restoration"),
      GeneratedBy = list(list(
        Name = "neuroarchive", Version = .restore_package_version(),
        Description = "Regenerated NIfTI images from an LNA exchange archive"
      ))
    ))
  }
  list(
    Name = paste0(pipeline, " Neuroarchive restoration"),
    BIDSVersion = .lna_restoration_bids_version,
    DatasetType = "derivative",
    GeneratedBy = list(
      list(
        Name = pipeline,
        Description = "Original derivative pipeline identity preserved by Neuroarchive"
      ),
      list(
        Name = "neuroarchive", Version = .restore_package_version(),
        Description = "Regenerated NIfTI images from LNA artifacts"
      )
    ),
    SourceDatasets = list(list(
      URL = paste0("file://", dataset$root),
      Version = dataset$declaration$profile_version
    )),
    DatasetLinks = list(neuroarchive = paste0("file://", dataset$root))
  )
}

.restore_write_json <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(
    x, path, auto_unbox = TRUE, null = "null", digits = NA, pretty = TRUE
  )
  invisible(path)
}

.restore_participants <- function(dataset, selected, stage) {
  receipt <- dataset$declaration$participants_receipt
  source <- file.path(dataset$root, "participants.tsv")
  if (is.null(receipt) || !file.exists(source) || !identical(
    .journal_file_checksum(source), receipt$output_checksum
  )) return(NULL)
  table <- utils::read.delim(
    source, colClasses = "character", check.names = FALSE,
    stringsAsFactors = FALSE, na.strings = character()
  )
  ids <- sort(unique(selected$entity_sub))
  ids <- paste0("sub-", ids)
  keep <- match(ids, table$participant_id)
  if (anyNA(keep)) {
    .restore_abort("archive participants table lacks selected subjects", "participants")
  }
  subset <- table[keep, , drop = FALSE]
  destination <- file.path(stage, "participants.tsv")
  utils::write.table(
    subset, destination, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = "n/a", fileEncoding = "UTF-8"
  )
  list(
    destination = "participants.tsv",
    selected_participant_ids = ids,
    columns = names(subset),
    checksum = .journal_file_checksum(destination)
  )
}

.restore_receipt <- function(
    dataset,
    closure,
    match,
    selection,
    requested_ids,
    selected_ids,
    inventory,
    companions,
    participants,
    missing) {
  list(
    schema_version = .lna_restoration_schema_version,
    status = closure,
    bids_validity = if (identical(closure, "closed")) {
      "closed-structure-intended-for-validation"
    } else "not-claimed-partial-dataset",
    source_archive = list(
      plan_hash = dataset$declaration$plan_hash,
      profile = dataset$declaration$profile,
      profile_version = dataset$declaration$profile_version
    ),
    selection = list(
      match = match,
      filters = selection,
      requested_record_ids = requested_ids,
      restored_record_ids = selected_ids
    ),
    outputs = list(
      nifti_paths = unname(inventory$nifti),
      metadata_sidecars = unname(inventory$sidecars),
      copied_associations = if (length(companions)) {
        vapply(companions, `[[`, character(1), "destination_relative")
      } else character(),
      participants = participants
    ),
    transformation_rules = list(
      imaging = "regenerated-from-lna-not-byte-identical-to-source",
      metadata = "effective-metadata-flattened-to-image-sidecar",
      events = "copied-byte-for-byte-when-present",
      entities = "preserved-exactly",
      original_file_byte_identity = "not-claimed"
    ),
    missing_dependencies = missing
  )
}

.restore_validate_receipt <- function(receipt) {
  required <- c(
    "schema_version", "status", "bids_validity", "source_archive", "selection",
    "outputs", "transformation_rules", "missing_dependencies"
  )
  if (!is.list(receipt) || !setequal(names(receipt), required) ||
      !identical(receipt$schema_version, .lna_restoration_schema_version) ||
      !receipt$status %in% c("closed", "partial")) {
    .restore_abort("restoration receipt is incompatible", "receipt")
  }
  if (identical(receipt$status, "closed") && length(receipt$missing_dependencies)) {
    .restore_abort("closed receipt may not contain missing dependencies", "receipt")
  }
  rules <- receipt$transformation_rules
  if (!identical(rules$original_file_byte_identity, "not-claimed")) {
    .restore_abort("restoration receipt makes an invalid byte-identity claim", "receipt")
  }
  .source_record_json_value(receipt, "receipt")
  receipt
}

#' Restore an LNA exchange archive to a BIDS-shaped NIfTI dataset
#'
#' @param dataset Dataset root or open handle.
#' @param output_root New restoration destination.
#' @param selection Named exact/regex/glob filters accepted by [query_lna()].
#' @param match Explicit filter matching mode.
#' @param closure `"closed"` computes transitive record and companion closure
#'   and refuses missing requirements before writing. `"partial"` restores the
#'   available explicit selection and marks that no BIDS-validity claim is made.
#' @param overwrite Whether an existing destination may be transactionally
#'   replaced.
#' @param hooks Optional reconstruction, validation, and companion-copy hooks.
#' @return A validated restoration receipt, invisibly.
#' @export
restore_bids <- function(
    dataset,
    output_root,
    selection = list(),
    match = c("exact", "regex", "glob"),
    closure = c("closed", "partial"),
    overwrite = FALSE,
    hooks = NULL) {
  match <- match.arg(match)
  closure <- match.arg(closure)
  dataset <- .query_dataset(dataset)
  hooks <- .restore_hooks(hooks)
  output_root <- .execution_scalar_string(output_root, "output_root")
  if (file.exists(output_root) || dir.exists(output_root)) {
    if (!isTRUE(overwrite)) {
      .restore_abort("output_root already exists", "output_root")
    }
  }
  requested <- .restore_selection(dataset, selection, match)
  transitive <- .restore_transitive_records(
    dataset$manifest$data, requested$record_id, closure
  )
  if (identical(closure, "closed") && length(transitive$missing)) {
    .restore_abort(
      "closed restoration has missing record dependencies", "closure",
      missing_dependencies = transitive$missing
    )
  }
  selected <- dataset$manifest$data[
    dataset$manifest$data$record_id %in% transitive$record_ids, , drop = FALSE
  ]
  preflight <- .restore_preflight(dataset, selected, closure)
  selected <- preflight$selected
  missing <- c(transitive$missing, preflight$missing)
  if (!nrow(selected)) {
    .restore_abort("no selected imaging artifacts are restorable", "closure")
  }
  companions <- .restore_companion_plan(dataset, selected, closure, missing)
  participants_expected <- if (
    !is.null(dataset$declaration$participants_receipt)
  ) TRUE else NULL
  inventory <- .restore_output_inventory(
    selected, companions, participants_expected
  )

  parent <- dirname(output_root)
  if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
    .restore_abort("could not create output parent directory", "output_root")
  }
  stage <- tempfile(".lna-restore-", tmpdir = parent)
  if (!dir.create(stage, recursive = FALSE, showWarnings = FALSE)) {
    .restore_abort("could not create restoration staging directory", "output_root")
  }
  installed <- FALSE
  on.exit(if (!installed && dir.exists(stage)) {
    unlink(stage, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  for (i in seq_len(nrow(selected))) {
    row <- selected[i, , drop = FALSE]
    source <- .restore_archive_path(dataset, row$lna_path[[1L]])
    nifti_relative <- .restore_nifti_path(row$lna_path[[1L]])
    destination <- .writer_full_path(stage, nifti_relative)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    temporary <- tempfile(
      pattern = paste0(".", basename(destination), "."),
      tmpdir = dirname(destination), fileext = ".nii.gz"
    )
    on.exit(if (file.exists(temporary)) unlink(temporary, force = TRUE), add = TRUE)
    hooks$reconstruct(source, temporary, row, dataset)
    if (!file.exists(temporary)) {
      .restore_abort("reconstruction did not produce a NIfTI file", "nifti")
    }
    hooks$validate_nifti(temporary, row)
    if (!file.rename(temporary, destination)) {
      .restore_abort("could not install regenerated NIfTI", "nifti")
    }
    metadata <- row$metadata_json[[1L]]
    if (identical(row$scope[[1L]], "derivative")) {
      metadata$Sources <- paste0("bids:neuroarchive:", row$lna_path[[1L]])
    }
    .restore_write_json(metadata, .writer_full_path(
      stage, .restore_sidecar_path(nifti_relative)
    ))
  }

  for (companion in companions) {
    destination <- .writer_full_path(stage, companion$destination_relative)
    if (file.exists(destination)) next
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    if (!isTRUE(hooks$copy_companion(companion$source, destination)) ||
        !file.exists(destination) || !identical(
          .journal_file_checksum(destination), companion$checksum
        )) {
      .restore_abort("companion copy failed checksum validation", "companions")
    }
  }

  participants <- .restore_participants(dataset, selected, stage)
  .restore_write_json(
    .restore_dataset_description("raw", dataset = dataset),
    file.path(stage, "dataset_description.json")
  )
  for (pipeline in inventory$pipelines) {
    .restore_write_json(
      .restore_dataset_description("derivative", pipeline, dataset),
      file.path(stage, "derivatives", pipeline, "dataset_description.json")
    )
  }
  readme <- c(
    "Neuroarchive restoration",
    "",
    paste0("Status: ", toupper(closure)),
    if (identical(closure, "partial")) {
      "This is an explicitly partial extraction and carries no BIDS-validity claim."
    } else {
      "This restoration passed Neuroarchive dependency-closure checks."
    },
    paste0(
      "NIfTI files were regenerated from LNA artifacts; original compressed-file ",
      "byte identity is not claimed."
    )
  )
  writeLines(readme, file.path(stage, "README"), useBytes = TRUE)
  writeLines("lna_restoration.json", file.path(stage, ".bidsignore"),
    useBytes = TRUE)
  receipt <- .restore_validate_receipt(.restore_receipt(
    dataset, closure, match, selection, requested$record_id, selected$record_id,
    inventory, companions, participants, missing
  ))
  .restore_write_json(receipt, file.path(stage, "lna_restoration.json"))

  backup <- NULL
  if (dir.exists(output_root) || file.exists(output_root)) {
    backup <- tempfile(".lna-restore-backup-", tmpdir = parent)
    if (!file.rename(output_root, backup)) {
      .restore_abort("could not preserve existing output_root", "output_root")
    }
  }
  on.exit({
    if (!installed && !is.null(backup) && file.exists(backup)) {
      if (file.exists(output_root)) unlink(output_root, recursive = TRUE, force = TRUE)
      file.rename(backup, output_root)
    }
  }, add = TRUE)
  if (!file.rename(stage, output_root)) {
    .restore_abort("could not atomically install restoration", "output_root")
  }
  installed <- TRUE
  if (!is.null(backup) && file.exists(backup)) {
    unlink(backup, recursive = TRUE, force = TRUE)
  }
  invisible(receipt)
}

#' Read a Neuroarchive restoration receipt
#'
#' @param path Restoration root or `lna_restoration.json` path.
#' @return A validated restoration receipt.
#' @export
read_lna_restoration <- function(path) {
  path <- .execution_scalar_string(path, "path")
  if (!identical(basename(path), "lna_restoration.json")) {
    path <- file.path(path, "lna_restoration.json")
  }
  if (!file.exists(path)) .restore_abort("restoration receipt is absent", "path")
  decoded <- tryCatch(
    jsonlite::read_json(path, simplifyVector = FALSE),
    error = function(e) .restore_abort("restoration receipt is invalid JSON", "path")
  )
  .restore_validate_receipt(decoded)
}
