# Transactional per-record writer built on plans, journals, and manifests.

.writer_abort <- function(message, field, stage = NULL, record_id = NULL, ...) {
  abort_lna(
    message,
    ...,
    field = field,
    stage = stage,
    record_id = record_id,
    .subclass = "lna_error_dataset_writer",
    location = paste0("write_lna_dataset:", field)
  )
}

.writer_default_hooks <- function() {
  list(
    check_space = function(plan, root) invisible(TRUE),
    verify_source = function(entry, source_root) {
      .writer_verify_source(entry, source_root)
    },
    fit_resource = function(resource, plan, archive_root) {
      .dataset_verify_shared_basis_resource(resource, plan, archive_root)
    },
    load_source = function(entry, source_path) {
      if (is.null(source_path)) {
        .writer_abort(
          "default source loader requires a local source path", "source",
          stage = "load", record_id = entry$record_id
        )
      }
      neuroim2::read_vec(source_path)
    },
    write_artifact = function(data, temp_path, entry, plan) {
      geometry <- entry$geometry
      transform_params <- entry$codec_transform_params %||%
        plan$codec_profile$transform_params
      compression_mask <- NULL
      if (identical(entry$shared_basis$status %||% "none", "shared")) {
        basis_file <- resolve_basis(
          entry$shared_basis$basis_id,
          entry$shared_basis$basis_checksum,
          locator = entry$shared_basis$basis_locator,
          archive_path = temp_path
        )
        basis <- read_shared_basis(basis_file, verify_checksum = TRUE)
        compression_mask <- if (identical(basis$support_type, "masked")) {
          basis$mask
        } else NULL
      }
      write_lna(
        data,
        temp_path,
        transforms = unlist(plan$codec_profile$transforms, use.names = FALSE),
        transform_params = transform_params,
        mask = compression_mask,
        run_id = paste0("run-", entry$entities$run %||% "01"),
        checksum = "sha256",
        run_geometry = list(
          dims = geometry$dims,
          transform = geometry$affine,
          spatial_units = geometry$spatial_units,
          voxel_ordering = geometry$voxel_ordering,
          n_volumes = geometry$n_volumes,
          grid_id = geometry$grid_id
        ),
        geometry_policy = "fmri"
      )
    },
    validate_artifact = function(path, entry, stage) {
      validate_lna(path, strict = TRUE, checksum = TRUE)
    },
    install_artifact = function(temp_path, destination, entry) {
      file.rename(temp_path, destination)
    },
    write_manifest = function(manifest, staging_root) {
      write_lna_manifest(manifest, staging_root)
    },
    resolve_companion = function(companion, source_root, archive_root) {
      .writer_full_path(source_root, companion$source_path)
    },
    copy_companion = function(source, temp_path, companion) {
      file.copy(source, temp_path, overwrite = FALSE, copy.mode = TRUE)
    },
    write_participants = function(plan, source_root, output_root, overwrite) {
      .metadata_write_participants(plan, source_root, output_root, overwrite)
    }
  )
}

.writer_hooks <- function(hooks) {
  defaults <- .writer_default_hooks()
  if (is.null(hooks)) return(defaults)
  if (!is.list(hooks) || is.data.frame(hooks) || is.null(names(hooks)) ||
      any(!nzchar(names(hooks))) || anyDuplicated(names(hooks))) {
    .writer_abort("hooks must be a uniquely named list", "hooks")
  }
  unknown <- setdiff(names(hooks), names(defaults))
  if (length(unknown)) {
    .writer_abort(
      paste0("unknown writer hook(s): ", paste(unknown, collapse = ", ")),
      "hooks"
    )
  }
  for (name in names(hooks)) {
    if (!is.function(hooks[[name]])) {
      .writer_abort(paste0("hooks.", name, " must be a function"), "hooks")
    }
    defaults[[name]] <- hooks[[name]]
  }
  defaults
}

.writer_resolve_source <- function(entry, source_root) {
  path <- entry$source$path
  if (is.null(path)) return(NULL)
  if (grepl("^/|^[A-Za-z]:", path)) return(path)
  file.path(
    source_root,
    do.call(file.path, as.list(strsplit(path, "/", fixed = TRUE)[[1L]]))
  )
}

.writer_verify_source <- function(entry, source_root) {
  path <- .writer_resolve_source(entry, source_root)
  if (is.null(path)) {
    if (is.null(entry$source$uri)) {
      .writer_abort(
        "source has neither path nor URI", "source", stage = "preflight",
        record_id = entry$record_id
      )
    }
    return(NULL)
  }
  if (!file.exists(path)) {
    .writer_abort(
      paste0("source file does not exist: ", path), "source",
      stage = "preflight", record_id = entry$record_id
    )
  }
  checksum <- .journal_file_checksum(path)
  if (!identical(checksum, entry$source$checksum)) {
    .writer_abort(
      "source checksum does not match the dataset plan", "source_checksum",
      stage = "preflight", record_id = entry$record_id,
      expected_checksum = entry$source$checksum, actual_checksum = checksum
    )
  }
  path
}

.writer_full_path <- function(root, relative) {
  file.path(
    root,
    do.call(file.path, as.list(strsplit(relative, "/", fixed = TRUE)[[1L]]))
  )
}

.writer_relative_path <- function(root, path) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  prefix <- paste0(root, "/")
  if (!startsWith(path, prefix)) {
    .writer_abort("temporary path escaped archive root", "temporary_path")
  }
  substring(path, nchar(prefix) + 1L)
}

.writer_atomic_install <- function(
    hooks,
    temp_path,
    destination,
    entry,
    replace,
    post_install = function(path) invisible(TRUE)) {
  if (!file.exists(temp_path)) {
    .writer_abort(
      "artifact writer did not produce its temporary output", "temporary_path",
      stage = "write", record_id = entry$record_id
    )
  }
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  backup <- NULL
  if (file.exists(destination)) {
    if (!isTRUE(replace)) {
      .writer_abort(
        paste0("destination already exists: ", entry$destination_relative),
        "destination", stage = "install", record_id = entry$record_id
      )
    }
    backup <- tempfile(
      pattern = paste0(".", basename(destination), ".backup-"),
      tmpdir = dirname(destination)
    )
    if (!file.rename(destination, backup)) {
      .writer_abort(
        "could not preserve the existing artifact before replacement",
        "destination", stage = "install", record_id = entry$record_id
      )
    }
  }
  installed <- FALSE
  on.exit({
    if (!installed) {
      if (file.exists(destination)) unlink(destination, force = TRUE)
      if (!is.null(backup) && file.exists(backup)) {
        file.rename(backup, destination)
      }
    }
  }, add = TRUE)
  ok <- isTRUE(hooks$install_artifact(temp_path, destination, entry))
  if (!ok || !file.exists(destination)) {
    .writer_abort(
      "atomic artifact installation failed", "destination",
      stage = "install", record_id = entry$record_id
    )
  }
  post_install(destination)
  installed <- TRUE
  if (file.exists(temp_path)) unlink(temp_path, force = TRUE)
  if (!is.null(backup) && file.exists(backup)) unlink(backup, force = TRUE)
  invisible(destination)
}

.writer_fail_event <- function(handle, node_id, error) {
  tryCatch(
    append_lna_journal_event(handle, node_id, "failed", failure = error),
    error = function(journal_error) invisible(NULL)
  )
  invisible(NULL)
}

.writer_run_node <- function(handle, node_id, action) {
  resume <- resume_lna_journal(handle, invalid = "rebuild")
  if (node_id %in% resume$state$completed_nodes) return(invisible(FALSE))
  append_lna_journal_event(handle, node_id, "started")
  tryCatch(
    {
      action()
      append_lna_journal_event(handle, node_id, "completed")
      invisible(TRUE)
    },
    error = function(e) {
      .writer_fail_event(handle, node_id, e)
      stop(e)
    }
  )
}

.writer_write_one <- function(
    handle,
    entry,
    source_path,
    hooks,
    overwrite) {
  node_id <- paste0("write:", entry$record_id)
  resume <- resume_lna_journal(handle, invalid = "rebuild")
  if (node_id %in% resume$state$completed_nodes) return("reused")
  destination <- .writer_full_path(handle$paths$root, entry$destination_relative)
  replace <- identical(overwrite, "rebuild")
  if (file.exists(destination) && !replace) {
    .writer_abort(
      paste0("destination already exists without a reusable journal receipt: ",
        entry$destination_relative),
      "destination", stage = "preflight", record_id = entry$record_id
    )
  }
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  temp_path <- tempfile(
    pattern = paste0(".", basename(destination), ".", handle$owner_id, "-"),
    tmpdir = dirname(destination), fileext = ".lna-partial"
  )
  append_lna_journal_event(handle, node_id, "started")
  completed <- FALSE
  on.exit({
    if (!completed && file.exists(temp_path)) unlink(temp_path, force = TRUE)
  }, add = TRUE)
  tryCatch(
    {
      data <- hooks$load_source(entry, source_path)
      hooks$write_artifact(data, temp_path, entry, handle$plan)
      if (!file.exists(temp_path)) {
        .writer_abort(
          "artifact writer did not create a temporary LNA", "temporary_path",
          stage = "write", record_id = entry$record_id
        )
      }
      temp_relative <- .writer_relative_path(handle$paths$root, temp_path)
      append_lna_journal_event(
        handle, node_id, "temporary-written", output_path = temp_relative
      )
      hooks$validate_artifact(temp_path, entry, "temporary")
      .writer_atomic_install(
        hooks, temp_path, destination, entry, replace = replace,
        post_install = function(path) {
          hooks$validate_artifact(path, entry, "installed")
        }
      )
      append_lna_journal_event(handle, node_id, "completed")
      completed <- TRUE
      "written"
    },
    error = function(e) {
      .writer_fail_event(handle, node_id, e)
      stop(e)
    }
  )
}

.writer_copy_companions <- function(handle, source_root, hooks, overwrite) {
  companions <- unlist(lapply(handle$plan$artifacts, function(entry) {
    entry$expected_companions
  }), recursive = FALSE)
  if (!length(companions)) return(invisible(character()))
  destinations <- vapply(companions, function(x) x$destination_relative %||% "", character(1))
  keep <- nzchar(destinations)
  companions <- companions[keep]
  destinations <- destinations[keep]
  if (!length(companions)) return(invisible(character()))
  groups <- split(seq_along(companions), destinations)
  copied <- character()
  for (destination_relative in names(groups)) {
    items <- companions[groups[[destination_relative]]]
    checksums <- unique(vapply(items, `[[`, character(1), "source_checksum"))
    if (length(checksums) != 1L) {
      .writer_abort(
        paste0("companion destination has conflicting checksums: ", destination_relative),
        "companions", stage = "finalize"
      )
    }
    item <- items[[1L]]
    if (is.null(item$source_path)) {
      .writer_abort(
        "default companion copier requires source_path", "companions",
        stage = "finalize"
      )
    }
    source <- hooks$resolve_companion(item, source_root, handle$paths$root)
    if (!file.exists(source) || !identical(.journal_file_checksum(source), checksums)) {
      .writer_abort(
        paste0("companion is missing or has changed: ", item$source_path),
        "companions", stage = "finalize"
      )
    }
    destination <- .writer_full_path(handle$paths$root, destination_relative)
    if (file.exists(destination)) {
      if (identical(.journal_file_checksum(destination), checksums)) {
        copied <- c(copied, destination_relative)
        next
      }
      if (!identical(overwrite, "rebuild")) {
        .writer_abort(
          paste0("companion destination already exists: ", destination_relative),
          "companions", stage = "finalize"
        )
      }
    }
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    temp <- tempfile(
      pattern = paste0(".", basename(destination), ".companion-"),
      tmpdir = dirname(destination)
    )
    on.exit(if (file.exists(temp)) unlink(temp, force = TRUE), add = TRUE)
    if (!isTRUE(hooks$copy_companion(source, temp, item)) ||
        !file.exists(temp) || !identical(.journal_file_checksum(temp), checksums)) {
      .writer_abort("companion copy or checksum verification failed", "companions")
    }
    pseudo_entry <- list(
      record_id = paste0("companion:", basename(destination_relative)),
      destination_relative = destination_relative
    )
    .writer_atomic_install(
      hooks, temp, destination, pseudo_entry,
      replace = identical(overwrite, "rebuild")
    )
    copied <- c(copied, destination_relative)
  }
  invisible(copied)
}

.writer_install_manifest <- function(handle, manifest, declaration, hooks) {
  staging <- tempfile(
    pattern = paste0(".lna-finalize-", handle$owner_id, "-"),
    tmpdir = handle$paths$root
  )
  dir.create(staging, recursive = FALSE, showWarnings = FALSE)
  on.exit(unlink(staging, recursive = TRUE, force = TRUE), add = TRUE)
  hooks$write_manifest(manifest, staging)
  write_lna_dataset_declaration(declaration, staging)
  staged <- read_lna_manifest(staging)
  staged_declaration <- read_lna_dataset_declaration(staging)
  if (!identical(staged$schema, manifest$schema) ||
      !isTRUE(all.equal(staged$data, manifest$data)) ||
      !identical(staged_declaration, declaration)) {
    .writer_abort("staged dataset authority failed round-trip validation", "manifest")
  }
  declaration_source <- file.path(staging, "lna_dataset.json")
  json_source <- file.path(staging, "lna_manifest.json")
  tsv_source <- file.path(staging, "lna_manifest.tsv")
  declaration_destination <- file.path(handle$paths$root, "lna_dataset.json")
  json_destination <- file.path(handle$paths$root, "lna_manifest.json")
  tsv_destination <- file.path(handle$paths$root, "lna_manifest.tsv")
  destinations <- c(
    declaration_destination, json_destination, tsv_destination
  )
  backups <- stats::setNames(vector("list", length(destinations)), destinations)
  installed <- character()
  committed <- FALSE
  on.exit({
    if (!committed) {
      for (destination in installed[file.exists(installed)]) unlink(destination, force = TRUE)
      for (destination in names(backups)) {
        backup <- backups[[destination]]
        if (!is.null(backup) && file.exists(backup)) file.rename(backup, destination)
      }
    }
  }, add = TRUE)
  for (destination in destinations[file.exists(destinations)]) {
    backup <- tempfile(
      pattern = paste0(".", basename(destination), ".backup-"),
      tmpdir = dirname(destination)
    )
    if (!file.rename(destination, backup)) {
      .writer_abort("could not preserve existing manifest during finalization", "manifest")
    }
    backups[[destination]] <- backup
  }
  if (!file.rename(declaration_source, declaration_destination)) {
    .writer_abort("could not install dataset declaration", "manifest")
  }
  installed <- c(installed, declaration_destination)
  if (!file.rename(json_source, json_destination)) {
    .writer_abort("could not install manifest JSON sidecar", "manifest")
  }
  installed <- c(installed, json_destination)
  if (!file.rename(tsv_source, tsv_destination)) {
    .writer_abort("could not install manifest TSV authority", "manifest")
  }
  installed <- c(installed, tsv_destination)
  committed <- TRUE
  for (backup in backups) if (!is.null(backup) && file.exists(backup)) unlink(backup, force = TRUE)
  read_lna_manifest(handle$paths$root)
}

#' Write a transactional BIDS-shaped LNA exchange dataset
#'
#' @param records Source records or an existing validated dataset plan.
#' @param output_root Destination root; omitted when `records` is a plan.
#' @param source_root Root used to resolve source and companion paths.
#' @param codec_profile Codec profile forwarded to [plan_lna_dataset()].
#' @param fit_resources Shared fit resources forwarded to the planner.
#' @param metadata_policy Portable metadata and privacy policy forwarded to the
#'   planner. The default requires provenance for every resolved metadata field
#'   and excludes participants unless explicitly opted in.
#' @param overwrite `"error"` protects every unreceipted destination;
#'   `"rebuild"` authorizes transactional replacement.
#' @param hooks Optional named test/integration hooks for source loading,
#'   artifact writing/validation/installation, disk checks, companion copying,
#'   and manifest writing.
#' @param owner_id Optional journal writer identity.
#' @param recover_truncated,plan_mismatch,stale_after,takeover_stale Journal
#'   recovery and lock policies passed to [open_lna_journal()].
#' @param physical_root Optional physical transaction root. This is an advanced
#'   generation-building hook used by [update_lna_dataset()]; the serialized
#'   plan and all public paths continue to name `output_root`.
#' @return A `lna_dataset_write_result` containing plan, final manifest,
#'   written and reused record IDs, and journal resume state.
#' @export
write_lna_dataset <- function(
    records,
    output_root = NULL,
    source_root = ".",
    codec_profile = list(
      profile_id = "lna-independent-lossless-v1",
      loss_mode = "lossless",
      transforms = character(),
      transform_params = list()
    ),
    fit_resources = list(),
    metadata_policy = lna_metadata_policy(portable = TRUE),
    overwrite = c("error", "rebuild"),
    hooks = NULL,
    owner_id = NULL,
    recover_truncated = c("error", "truncate"),
    plan_mismatch = c("error", "restart"),
    stale_after = 300,
    takeover_stale = FALSE,
    physical_root = NULL) {
  overwrite <- match.arg(overwrite)
  recover_truncated <- match.arg(recover_truncated)
  plan_mismatch <- match.arg(plan_mismatch)
  source_root <- .execution_scalar_string(source_root, "source_root")
  hooks <- .writer_hooks(hooks)
  plan <- if (inherits(records, "lna_dataset_plan")) {
    if (!is.null(output_root) && !identical(output_root, records$output_root)) {
      .writer_abort("output_root disagrees with the supplied plan", "output_root")
    }
    validate_lna_dataset_plan(records)
  } else {
    output_root <- .execution_scalar_string(output_root, "output_root")
    plan_lna_dataset(
      records, output_root, codec_profile = codec_profile,
      fit_resources = fit_resources, metadata_policy = metadata_policy
    )
  }
  output_root <- plan$output_root
  archive_root <- if (is.null(physical_root)) {
    output_root
  } else {
    .execution_scalar_string(physical_root, "physical_root")
  }
  handle <- open_lna_journal(
    plan, archive_root, owner_id = owner_id, stale_after = stale_after,
    takeover_stale = takeover_stale, recover_truncated = recover_truncated,
    plan_mismatch = plan_mismatch
  )
  on.exit(close_lna_journal(handle), add = TRUE)
  journal_events <- .journal_events(handle)
  finalization_seen <- any(vapply(journal_events, function(event) {
    identical(event$node_id, "phase:finalize")
  }, logical(1)))
  reserved_outputs <- file.path(
    archive_root,
    c("lna_dataset.json", "lna_manifest.tsv", "lna_manifest.json", ".lna_plan.json")
  )
  if (any(file.exists(reserved_outputs)) && !finalization_seen &&
      identical(overwrite, "error")) {
    .writer_abort(
      "dataset authority files already exist without a journal receipt",
      "destination", stage = "preflight"
    )
  }
  sources <- stats::setNames(lapply(plan$artifacts, function(entry) {
    hooks$verify_source(entry, source_root)
  }), vapply(plan$artifacts, `[[`, character(1), "record_id"))
  hooks$check_space(plan, archive_root)
  .writer_run_node(handle, "phase:preflight", function() invisible(TRUE))
  for (resource in plan$execution$fit_resources) {
    .writer_run_node(
      handle, paste0("fit:", resource$resource_id), function() {
        hooks$fit_resource(resource, plan, archive_root)
      }
    )
  }
  .writer_run_node(handle, "phase:fit-complete", function() invisible(TRUE))
  outcomes <- vapply(plan$artifacts, function(entry) {
    .writer_write_one(
      handle, entry, sources[[entry$record_id]], hooks, overwrite
    )
  }, character(1))
  .writer_run_node(handle, "phase:write-complete", function() invisible(TRUE))
  for (entry in plan$artifacts) {
    destination <- .writer_full_path(archive_root, entry$destination_relative)
    .writer_run_node(handle, paste0("validate:", entry$record_id), function() {
      hooks$validate_artifact(destination, entry, "final")
    })
  }
  .writer_copy_companions(handle, source_root, hooks, overwrite)
  participants_receipt <- hooks$write_participants(
    plan, source_root, archive_root, overwrite
  )
  resume_before_finalize <- resume_lna_journal(handle, invalid = "error")
  if ("phase:finalize" %in% resume_before_finalize$state$completed_nodes) {
    manifest_ok <- tryCatch({
      read_lna_manifest(archive_root)
      declaration <- read_lna_dataset_declaration(archive_root)
      validate_lna_dataset_declaration(declaration, plan_hash = plan$plan_hash)
      TRUE
    }, error = function(e) FALSE)
    if (!manifest_ok) {
      latest <- .journal_latest_node_event(
        resume_before_finalize$events, "phase:finalize"
      )
      append_lna_journal_event(
        handle, "phase:finalize", "invalidated", attempt = latest$attempt,
        details = list(reason = "final manifest pair is missing or invalid")
      )
    }
  }
  final_manifest <- NULL
  .writer_run_node(handle, "phase:finalize", function() {
    manifest <- lna_manifest_from_journal(handle)
    declaration <- lna_dataset_declaration(plan, participants_receipt)
    write_lna_dataset_plan(
      plan, file.path(archive_root, ".lna_plan.json"), overwrite = TRUE
    )
    final_manifest <<- .writer_install_manifest(
      handle, manifest, declaration, hooks
    )
  })
  if (is.null(final_manifest)) final_manifest <- read_lna_manifest(archive_root)
  resume <- resume_lna_journal(handle, invalid = "error")
  outcome_ids <- vapply(plan$artifacts, `[[`, character(1), "record_id")
  structure(
    list(
      plan = plan,
      manifest = final_manifest,
      written_record_ids = outcome_ids[outcomes == "written"],
      reused_record_ids = outcome_ids[outcomes == "reused"],
      resume = resume
    ),
    class = c("lna_dataset_write_result", "list")
  )
}
