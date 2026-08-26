# Checksummed, resumable local transport for complete LNA exchange datasets.

.lna_transfer_schema_version <- "1.0.0"

.transfer_abort <- function(message, field, ...) {
  abort_lna(
    message,
    ...,
    field = field,
    .subclass = "lna_error_dataset_transfer",
    location = paste0("transfer_lna_dataset:", field)
  )
}

#' Compute the canonical SHA-256 identity of a file
#'
#' @param path Existing local file.
#' @return A string in `sha256:<64 hexadecimal digits>` form.
#' @export
lna_file_checksum <- function(path) {
  path <- .execution_scalar_string(path, "path")
  if (!file.exists(path) || dir.exists(path)) {
    .transfer_abort("path must name an existing file", "path")
  }
  .journal_file_checksum(path)
}

.transfer_relative_files <- function(root) {
  paths <- list.files(
    root, recursive = TRUE, full.names = TRUE, all.files = TRUE,
    include.dirs = FALSE, no.. = TRUE
  )
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  normalized <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  relative <- substring(normalized, nchar(root) + 2L)
  relative <- gsub("\\\\", "/", relative)
  volatile <- grepl(
    "(^|/)[.]lna[.]lock$|(^|/)[.]lna-update[.]lock($|/)|lna-partial|[.]backup-",
    relative
  )
  reserved <- basename(relative) %in% c(
    "lna_transfer.json", ".lna_transfer_journal.jsonl"
  )
  relative[!volatile & !reserved]
}

.transfer_portable_key <- function(path) {
  tolower(stringi::stri_trans_nfc(enc2utf8(path)))
}

#' Validate transfer paths for case and Unicode-normalization collisions
#'
#' @param paths Archive-relative paths.
#' @return `paths`, invisibly, or an error before transfer.
#' @export
validate_lna_transfer_paths <- function(paths) {
  if (!is.character(paths) || anyNA(paths) || any(!nzchar(paths))) {
    .transfer_abort("paths must contain non-empty strings", "paths")
  }
  normalized <- vapply(paths, .manifest_validate_path, character(1), field = "paths")
  keys <- vapply(normalized, .transfer_portable_key, character(1))
  duplicated_keys <- unique(keys[duplicated(keys)])
  if (length(duplicated_keys)) {
    collisions <- lapply(duplicated_keys, function(key) normalized[keys == key])
    .transfer_abort(
      "archive paths collide after case folding or Unicode normalization",
      "paths", collisions = collisions
    )
  }
  invisible(normalized)
}

.transfer_file_record <- function(root, relative) {
  path <- .writer_full_path(root, relative)
  list(
    path = relative,
    bytes = as.numeric(file.info(path)$size),
    checksum = .journal_file_checksum(path)
  )
}

.transfer_inventory_payload <- function(inventory) {
  payload <- unclass(inventory)
  payload$inventory_hash <- NULL
  payload$source_root <- NULL
  payload
}

.transfer_checksum <- function(x, field) {
  tryCatch(
    .source_record_checksum(x, field),
    lna_error_source_record = function(e) {
      .transfer_abort(conditionMessage(e), field)
    }
  )
}

.transfer_nonnegative_integer <- function(x, field) {
  if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
  x <- suppressWarnings(as.numeric(x))
  if (length(x) != 1L || !is.finite(x) || x < 0 || x != floor(x)) {
    .transfer_abort(paste0(field, " must be one non-negative integer"), field)
  }
  x
}

.transfer_file_records <- function(files, field = "files") {
  if (!is.list(files) || is.data.frame(files) || !length(files)) {
    .transfer_abort(paste0(field, " must be a non-empty array"), field)
  }
  required <- c("path", "bytes", "checksum")
  normalized <- lapply(seq_along(files), function(i) {
    item <- files[[i]]
    item_field <- paste0(field, "[[", i, "]]")
    if (!is.list(item) || is.data.frame(item) || is.null(names(item)) ||
        length(names(item)) != length(required) || anyDuplicated(names(item)) ||
        !setequal(names(item), required)) {
      .transfer_abort(
        paste0(item_field, " must contain exactly: ", paste(required, collapse = ", ")),
        item_field
      )
    }
    path <- tryCatch(
      .manifest_validate_path(item$path, paste0(item_field, ".path")),
      lna_error_manifest = function(e) {
        .transfer_abort(conditionMessage(e), paste0(item_field, ".path"))
      }
    )
    list(
      path = path,
      bytes = .transfer_nonnegative_integer(
        item$bytes, paste0(item_field, ".bytes")
      ),
      checksum = .transfer_checksum(
        item$checksum, paste0(item_field, ".checksum")
      )
    )
  })
  paths <- vapply(normalized, `[[`, character(1), "path")
  if (anyDuplicated(paths)) {
    .transfer_abort(paste0(field, " contains duplicate paths"), field)
  }
  validate_lna_transfer_paths(paths)
  normalized[order(paths)]
}

.transfer_summary <- function(summary, files, field = "summary") {
  required <- c("n_files", "total_bytes")
  if (!is.list(summary) || is.data.frame(summary) || is.null(names(summary)) ||
      length(names(summary)) != length(required) || anyDuplicated(names(summary)) ||
      !setequal(names(summary), required)) {
    .transfer_abort(
      paste0(field, " must contain exactly: ", paste(required, collapse = ", ")),
      field
    )
  }
  normalized <- list(
    n_files = .transfer_nonnegative_integer(summary$n_files, paste0(field, ".n_files")),
    total_bytes = .transfer_nonnegative_integer(
      summary$total_bytes, paste0(field, ".total_bytes")
    )
  )
  expected_files <- length(files)
  expected_bytes <- sum(vapply(files, `[[`, numeric(1), "bytes"))
  if (!identical(as.numeric(normalized$n_files), as.numeric(expected_files)) ||
      !identical(as.numeric(normalized$total_bytes), as.numeric(expected_bytes))) {
    .transfer_abort("transfer summary does not match its file records", field)
  }
  normalized
}

.transfer_inventory_hash <- function(schema_version, files, summary) {
  .dataset_plan_hash(list(
    schema_version = schema_version,
    files = unname(files),
    summary = summary
  ))
}

.transfer_required_paths <- function(manifest, declaration) {
  required <- c(
    "lna_dataset.json", "lna_manifest.json", "lna_manifest.tsv",
    ".lna_plan.json", ".lna_journal.jsonl",
    manifest$data$lna_path
  )
  companions <- unlist(lapply(manifest$data$associations_json, function(items) {
    if (!length(items)) return(character())
    vapply(items, function(item) item$destination_relative %||% "", character(1))
  }), use.names = FALSE)
  required <- c(required, companions[nzchar(companions)])
  shared <- declaration$fit_resources[vapply(
    declaration$fit_resources,
    function(resource) identical(resource$resource_type, "shared-basis"),
    logical(1)
  )]
  if (length(shared)) {
    required <- c(
      required,
      "bases/registry.json",
      vapply(shared, function(resource) {
        paste0("bases/", resource$parameters$basis_filename)
      }, character(1))
    )
  }
  if (!is.null(declaration$participants_receipt)) {
    required <- c(required, declaration$participants_receipt$destination)
  }
  required <- unique(required[nzchar(required)])
  validate_lna_transfer_paths(required)
  sort(required)
}

.transfer_validate_archive_closure <- function(root, manifest, declaration, records) {
  if (!all(manifest$data$status == "complete")) {
    .transfer_abort("portable manifest contains incomplete artifact rows", "manifest.status")
  }
  paths <- names(records)
  required <- .transfer_required_paths(manifest, declaration)
  missing <- setdiff(required, paths)
  if (length(missing)) {
    .transfer_abort(
      paste0("portable archive is missing required file(s): ", paste(missing, collapse = ", ")),
      "inventory.files", missing_paths = missing
    )
  }
  plan <- read_lna_dataset_plan(file.path(root, ".lna_plan.json"))
  if (!identical(plan$plan_hash, declaration$plan_hash)) {
    .transfer_abort("dataset plan and declaration hashes disagree", "plan_hash")
  }
  read_lna_journal(file.path(root, ".lna_journal.jsonl"))
  rows <- .update_manifest_rows(manifest)
  for (row in rows) {
    record_id <- .update_scalar(row, "record_id")
    relative <- .update_scalar(row, "lna_path")
    expected_checksum <- .update_scalar(row, "lna_checksum")
    expected_bytes <- .update_scalar(row, "lna_bytes")
    item <- records[[relative]]
    if (is.null(item) || is.null(expected_checksum) || is.null(expected_bytes) ||
        !identical(item$checksum, expected_checksum) ||
        !identical(as.numeric(item$bytes), as.numeric(expected_bytes))) {
      .transfer_abort(
        "manifest artifact size or checksum does not match the portable inventory",
        "manifest", record_id = record_id
      )
    }
    for (association in row$associations_json[[1L]]) {
      destination <- association$destination_relative %||% NULL
      if (is.null(destination)) next
      companion <- records[[destination]]
      if (is.null(companion) ||
          !identical(companion$checksum, association$source_checksum)) {
        .transfer_abort(
          "manifest association checksum does not match the portable inventory",
          "manifest.associations", record_id = record_id,
          path = destination
        )
      }
    }
  }
  for (resource in declaration$fit_resources) {
    if (identical(resource$resource_type, "shared-basis")) {
      .dataset_verify_shared_basis_resource(resource, plan, root)
    }
  }
  participants <- declaration$participants_receipt
  if (!is.null(participants)) {
    item <- records[[participants$destination]]
    if (is.null(item) ||
        !identical(item$checksum, participants$output_checksum) ||
        !identical(as.numeric(item$bytes), as.numeric(participants$output_bytes))) {
      .transfer_abort(
        "participants receipt does not match the portable inventory",
        "participants_receipt"
      )
    }
  }
  invisible(TRUE)
}

#' Inventory every byte-bearing file in an LNA exchange dataset
#'
#' The live manifest and declaration are validated first. Every portable file,
#' including LNA artifacts, companions, bases, and authority sidecars, receives
#' a byte size and SHA-256 digest. Operational locks and partial files are
#' excluded.
#'
#' @param source_root LNA dataset root.
#' @return An `lna_transfer_inventory`.
#' @export
lna_transfer_inventory <- function(source_root) {
  source_root <- .update_normalize_root(source_root, "source_root")
  manifest <- read_lna_manifest(source_root)
  declaration <- read_lna_dataset_declaration(source_root)
  validate_lna_dataset_declaration(declaration)
  paths <- sort(.transfer_relative_files(source_root))
  validate_lna_transfer_paths(paths)
  records <- lapply(paths, .transfer_file_record, root = source_root)
  names(records) <- paths
  .transfer_validate_archive_closure(source_root, manifest, declaration, records)
  payload <- list(
    schema_version = .lna_transfer_schema_version,
    source_root = source_root,
    files = unname(records),
    summary = list(
      n_files = length(records),
      total_bytes = sum(vapply(records, `[[`, numeric(1), "bytes"))
    )
  )
  payload$inventory_hash <- .dataset_plan_hash(.transfer_inventory_payload(payload))
  structure(payload, class = c("lna_transfer_inventory", "list"))
}

.transfer_validate_inventory <- function(inventory) {
  if (!inherits(inventory, "lna_transfer_inventory") || !is.list(inventory)) {
    .transfer_abort("inventory must come from lna_transfer_inventory()", "inventory")
  }
  required <- c("schema_version", "source_root", "files", "summary", "inventory_hash")
  if (!setequal(names(inventory), required) ||
      !identical(inventory$schema_version, .lna_transfer_schema_version)) {
    .transfer_abort("transfer inventory schema is incompatible", "inventory")
  }
  inventory$files <- .transfer_file_records(inventory$files, "inventory.files")
  inventory$summary <- .transfer_summary(
    inventory$summary, inventory$files, "inventory.summary"
  )
  inventory$inventory_hash <- .transfer_checksum(
    inventory$inventory_hash, "inventory.inventory_hash"
  )
  expected_hash <- .transfer_inventory_hash(
    inventory$schema_version, inventory$files, inventory$summary
  )
  if (!identical(expected_hash, inventory$inventory_hash)) {
    .transfer_abort("transfer inventory hash verification failed", "inventory_hash")
  }
  inventory
}

.transfer_validate_receipt <- function(receipt) {
  required <- c(
    "schema_version", "inventory_hash", "direction", "completed_at", "files", "summary"
  )
  if (!is.list(receipt) || is.data.frame(receipt) || is.null(names(receipt)) ||
      length(names(receipt)) != length(required) || anyDuplicated(names(receipt)) ||
      !setequal(names(receipt), required) ||
      !identical(receipt$schema_version, .lna_transfer_schema_version)) {
    .transfer_abort("transfer receipt schema is incompatible", "receipt")
  }
  direction <- .execution_scalar_string(receipt$direction, "receipt.direction")
  if (!direction %in% c("copy", "upload", "download")) {
    .transfer_abort("transfer receipt direction is invalid", "receipt.direction")
  }
  completed_at <- .execution_scalar_string(
    receipt$completed_at, "receipt.completed_at"
  )
  if (!grepl(
      "^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}([.][0-9]+)?Z$",
      completed_at)) {
    .transfer_abort("transfer receipt completion timestamp is invalid", "receipt.completed_at")
  }
  files <- .transfer_file_records(receipt$files, "receipt.files")
  summary <- .transfer_summary(receipt$summary, files, "receipt.summary")
  inventory_hash <- .transfer_checksum(
    receipt$inventory_hash, "receipt.inventory_hash"
  )
  expected_hash <- .transfer_inventory_hash(receipt$schema_version, files, summary)
  if (!identical(expected_hash, inventory_hash)) {
    .transfer_abort("transfer receipt inventory hash verification failed", "receipt.inventory_hash")
  }
  list(
    schema_version = receipt$schema_version,
    inventory_hash = inventory_hash,
    direction = direction,
    completed_at = completed_at,
    files = files,
    summary = summary
  )
}

.transfer_same_prefix <- function(source, partial, bytes, chunk_bytes) {
  if (bytes == 0) return(TRUE)
  source_con <- file(source, open = "rb")
  partial_con <- file(partial, open = "rb")
  on.exit(close(source_con), add = TRUE)
  on.exit(close(partial_con), add = TRUE)
  remaining <- bytes
  while (remaining > 0) {
    take <- min(chunk_bytes, remaining)
    left <- readBin(source_con, "raw", n = take)
    right <- readBin(partial_con, "raw", n = take)
    if (!identical(left, right)) return(FALSE)
    remaining <- remaining - length(left)
    if (!length(left) && remaining > 0) return(FALSE)
  }
  TRUE
}

.transfer_copy_stream <- function(
    source,
    partial,
    expected_bytes,
    chunk_bytes,
    after_chunk = function(...) invisible(TRUE)) {
  offset <- if (file.exists(partial)) as.numeric(file.info(partial)$size) else 0
  if (offset > expected_bytes ||
      !.transfer_same_prefix(source, partial, offset, chunk_bytes)) {
    unlink(partial, force = TRUE)
    offset <- 0
  }
  source_con <- file(source, open = "rb")
  on.exit(close(source_con), add = TRUE)
  seek(source_con, where = offset, origin = "start")
  destination_con <- file(partial, open = if (offset > 0) "ab" else "wb")
  destination_closed <- FALSE
  on.exit(if (!destination_closed) close(destination_con), add = TRUE)
  transferred <- offset
  peak_buffer_bytes <- 0
  while (transferred < expected_bytes) {
    buffer <- readBin(
      source_con, "raw", n = min(chunk_bytes, expected_bytes - transferred)
    )
    if (!length(buffer)) {
      .transfer_abort("sender ended before the declared byte size", "source")
    }
    peak_buffer_bytes <- max(peak_buffer_bytes, length(buffer))
    writeBin(buffer, destination_con)
    transferred <- transferred + length(buffer)
    after_chunk(
      source = source, partial = partial, transferred = transferred,
      expected_bytes = expected_bytes
    )
  }
  close(destination_con)
  destination_closed <- TRUE
  list(resumed_from = offset, peak_buffer_bytes = peak_buffer_bytes)
}

.transfer_journal_append <- function(root, inventory, relative, state, details = list()) {
  path <- file.path(root, ".lna_transfer_journal.jsonl")
  sequence <- if (file.exists(path)) length(readLines(path, warn = FALSE)) + 1L else 1L
  event <- list(
    schema_version = .lna_transfer_schema_version,
    sequence = sequence,
    inventory_hash = inventory$inventory_hash,
    path = relative,
    state = state,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
    details = details
  )
  event$event_id <- .dataset_plan_hash(event)
  cat(
    jsonlite::toJSON(event, auto_unbox = TRUE, null = "null", digits = NA),
    "\n", file = path, append = TRUE, sep = ""
  )
  invisible(event)
}

.transfer_order <- function(files) {
  authority <- c("lna_dataset.json", "lna_manifest.json", "lna_manifest.tsv")
  paths <- vapply(files, `[[`, character(1), "path")
  rank <- match(basename(paths), authority, nomatch = 0L)
  files[order(rank > 0L, rank, paths)]
}

.transfer_roots_overlap <- function(left, right) {
  identical(left, right) ||
    startsWith(left, paste0(right, "/")) ||
    startsWith(right, paste0(left, "/"))
}

.transfer_normalize_future_root <- function(path) {
  if (dir.exists(path)) {
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }
  suffix <- character()
  ancestor <- path
  while (!dir.exists(ancestor)) {
    parent <- dirname(ancestor)
    if (identical(parent, ancestor)) break
    suffix <- c(basename(ancestor), suffix)
    ancestor <- parent
  }
  ancestor <- normalizePath(ancestor, winslash = "/", mustWork = TRUE)
  if (!length(suffix)) ancestor else do.call(file.path, c(list(ancestor), as.list(suffix)))
}

.transfer_invalidate_stale_receipt <- function(root, inventory_hash) {
  path <- file.path(root, "lna_transfer.json")
  if (!file.exists(path)) return(invisible(FALSE))
  receipt <- tryCatch(
    .transfer_validate_receipt(jsonlite::read_json(path, simplifyVector = FALSE)),
    error = function(e) NULL
  )
  if (!is.null(receipt) && identical(receipt$inventory_hash, inventory_hash)) {
    return(invisible(FALSE))
  }
  status <- unlink(path, force = FALSE)
  if (!identical(status, 0L) && file.exists(path)) {
    .transfer_abort("could not invalidate stale receiver receipt", "receipt")
  }
  invisible(TRUE)
}

.transfer_write_receipt <- function(root, inventory, direction) {
  receipt <- list(
    schema_version = .lna_transfer_schema_version,
    inventory_hash = inventory$inventory_hash,
    direction = direction,
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
    files = inventory$files,
    summary = inventory$summary
  )
  jsonlite::write_json(
    receipt, file.path(root, "lna_transfer.json"), auto_unbox = TRUE,
    null = "null", digits = NA, pretty = TRUE
  )
  receipt
}

#' Verify a received LNA dataset against its transfer receipt
#'
#' @param destination_root Received dataset root.
#' @param error Whether any missing, truncated, or corrupt file is a hard error.
#' @param expected_inventory_hash Optional trusted sender inventory identity.
#'   Supplying it distinguishes an authentic sender inventory from a merely
#'   self-consistent receipt.
#' @return A verification receipt with per-file failures.
#' @export
verify_lna_transfer <- function(
    destination_root,
    error = TRUE,
    expected_inventory_hash = NULL) {
  destination_root <- .update_normalize_root(destination_root, "destination_root")
  if (!is.logical(error) || length(error) != 1L || is.na(error)) {
    .transfer_abort("error must be TRUE or FALSE", "error")
  }
  receipt_path <- file.path(destination_root, "lna_transfer.json")
  if (!file.exists(receipt_path)) {
    .transfer_abort("receiver has no lna_transfer.json receipt", "receipt")
  }
  receipt <- tryCatch(
    jsonlite::read_json(receipt_path, simplifyVector = FALSE),
    error = function(e) .transfer_abort("transfer receipt is invalid JSON", "receipt")
  )
  receipt <- .transfer_validate_receipt(receipt)
  if (!is.null(expected_inventory_hash)) {
    expected_inventory_hash <- .transfer_checksum(
      expected_inventory_hash, "expected_inventory_hash"
    )
    if (!identical(receipt$inventory_hash, expected_inventory_hash)) {
      .transfer_abort(
        "transfer receipt does not match the trusted sender inventory hash",
        "expected_inventory_hash",
        expected = expected_inventory_hash,
        actual = receipt$inventory_hash
      )
    }
  }
  failures <- list()
  for (item in receipt$files) {
    path <- .writer_full_path(destination_root, item$path)
    reason <- if (!file.exists(path)) {
      "missing"
    } else if (!identical(as.numeric(file.info(path)$size), as.numeric(item$bytes))) {
      "size-mismatch"
    } else if (!identical(.journal_file_checksum(path), item$checksum)) {
      "checksum-mismatch"
    } else NULL
    if (!is.null(reason)) failures[[item$path]] <- reason
  }
  actual_inventory <- tryCatch(
    lna_transfer_inventory(destination_root),
    error = function(e) e
  )
  if (inherits(actual_inventory, "error")) {
    failures[["<archive-closure>"]] <- paste0(
      "invalid-archive: ", conditionMessage(actual_inventory)
    )
  } else {
    actual_inventory <- .transfer_validate_inventory(actual_inventory)
    if (!identical(actual_inventory$inventory_hash, receipt$inventory_hash)) {
      expected_paths <- vapply(receipt$files, `[[`, character(1), "path")
      actual_paths <- vapply(actual_inventory$files, `[[`, character(1), "path")
      for (path in setdiff(expected_paths, actual_paths)) failures[[path]] <- "missing"
      for (path in setdiff(actual_paths, expected_paths)) failures[[path]] <- "unreceipted"
      if (!length(setdiff(c(expected_paths, actual_paths), names(failures)))) {
        failures[["<inventory>"]] <- "inventory-hash-mismatch"
      }
    }
  }
  passed <- !length(failures)
  result <- structure(list(
    passed = passed,
    inventory_hash = receipt$inventory_hash,
    n_files = length(receipt$files),
    authenticated = !is.null(expected_inventory_hash),
    failures = failures
  ), class = c("lna_transfer_verification", "list"))
  if (!passed && isTRUE(error)) {
    .transfer_abort(
      "receiver failed strict checksum verification", "verification",
      verification = result
    )
  }
  if (passed) {
    read_lna_manifest(destination_root)
    read_lna_dataset_declaration(destination_root)
  }
  result
}

#' Transfer an LNA dataset with bounded-memory resume and repair
#'
#' @param source_root Complete source archive.
#' @param destination_root Receiver directory.
#' @param direction Provenance label for copy, upload, or download paths.
#' @param chunk_bytes Maximum transfer buffer size.
#' @param resume Whether matching partial prefixes may resume.
#' @param repair Whether corrupt completed destinations may be replaced.
#' @param verify Must remain `TRUE`; strict receivers cannot disable digests.
#' @param hooks Optional `after_chunk` and `after_install` fault/integration hooks.
#' @return An `lna_dataset_transfer_result`.
#' @export
transfer_lna_dataset <- function(
    source_root,
    destination_root,
    direction = c("copy", "upload", "download"),
    chunk_bytes = 8 * 1024^2,
    resume = TRUE,
    repair = TRUE,
    verify = TRUE,
    hooks = list()) {
  direction <- match.arg(direction)
  if (!isTRUE(verify)) {
    .transfer_abort("checksum verification cannot be disabled", "verify")
  }
  if (!is.logical(resume) || length(resume) != 1L || is.na(resume) ||
      !is.logical(repair) || length(repair) != 1L || is.na(repair)) {
    .transfer_abort("resume and repair must be TRUE or FALSE", "policy")
  }
  chunk_bytes <- suppressWarnings(as.numeric(chunk_bytes))
  if (length(chunk_bytes) != 1L || !is.finite(chunk_bytes) ||
      chunk_bytes < 1 || chunk_bytes != floor(chunk_bytes)) {
    .transfer_abort("chunk_bytes must be one positive integer", "chunk_bytes")
  }
  if (!is.list(hooks) || is.data.frame(hooks)) .transfer_abort("hooks must be a list", "hooks")
  unknown_hooks <- setdiff(names(hooks), c("after_chunk", "after_install"))
  if (length(unknown_hooks)) .transfer_abort("unknown transfer hook", "hooks")
  after_chunk <- hooks$after_chunk %||% function(...) invisible(TRUE)
  after_install <- hooks$after_install %||% function(...) invisible(TRUE)
  if (!is.function(after_chunk) || !is.function(after_install)) {
    .transfer_abort("transfer hooks must be functions", "hooks")
  }
  inventory <- .transfer_validate_inventory(lna_transfer_inventory(source_root))
  source_root <- inventory$source_root
  destination_root <- .execution_scalar_string(destination_root, "destination_root")
  destination_root <- .transfer_normalize_future_root(destination_root)
  if (.transfer_roots_overlap(source_root, destination_root)) {
    .transfer_abort(
      "source and destination roots must be disjoint",
      "destination_root"
    )
  }
  if (!dir.exists(destination_root) &&
      !dir.create(destination_root, recursive = TRUE, showWarnings = FALSE)) {
    .transfer_abort("could not create destination_root", "destination_root")
  }
  destination_root <- normalizePath(destination_root, winslash = "/", mustWork = TRUE)
  .transfer_invalidate_stale_receipt(destination_root, inventory$inventory_hash)
  outcomes <- list()
  peak <- 0
  for (item in .transfer_order(inventory$files)) {
    source <- .writer_full_path(source_root, item$path)
    destination <- .writer_full_path(destination_root, item$path)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(destination)) {
      valid <- identical(as.numeric(file.info(destination)$size), as.numeric(item$bytes)) &&
        identical(.journal_file_checksum(destination), item$checksum)
      if (valid) {
        outcomes[[item$path]] <- "reused"
        .transfer_journal_append(destination_root, inventory, item$path, "reused")
        next
      }
      if (!isTRUE(repair)) {
        .transfer_abort("receiver contains a corrupt file and repair is FALSE", "repair",
          path = item$path)
      }
    }
    partial <- paste0(destination, ".lna-partial")
    if (!isTRUE(resume) && file.exists(partial)) unlink(partial, force = TRUE)
    copy <- tryCatch(
      .transfer_copy_stream(
        source, partial, item$bytes, chunk_bytes, after_chunk = after_chunk
      ),
      error = function(e) {
        .transfer_journal_append(
          destination_root, inventory, item$path, "interrupted",
          list(message = conditionMessage(e))
        )
        stop(e)
      }
    )
    peak <- max(peak, copy$peak_buffer_bytes)
    if (!identical(as.numeric(file.info(partial)$size), as.numeric(item$bytes)) ||
        !identical(.journal_file_checksum(partial), item$checksum)) {
      .transfer_journal_append(destination_root, inventory, item$path, "checksum-failed")
      .transfer_abort("transferred file failed checksum verification", "checksum",
        path = item$path)
    }
    backup <- NULL
    if (file.exists(destination)) {
      backup <- paste0(destination, ".transfer-backup")
      if (file.exists(backup)) unlink(backup, force = TRUE)
      if (!file.rename(destination, backup)) {
        .transfer_abort("could not preserve corrupt receiver file", "repair")
      }
    }
    installed <- FALSE
    if (!file.rename(partial, destination)) {
      if (!is.null(backup) && file.exists(backup)) file.rename(backup, destination)
      .transfer_abort("could not install verified receiver file", "destination")
    }
    installed <- TRUE
    if (!is.null(backup) && file.exists(backup)) unlink(backup, force = TRUE)
    after_install(path = item$path, destination = destination)
    outcomes[[item$path]] <- if (copy$resumed_from > 0) "resumed" else if (installed) "copied" else "failed"
    .transfer_journal_append(
      destination_root, inventory, item$path, "completed",
      list(resumed_from = copy$resumed_from)
    )
  }
  receipt <- .transfer_write_receipt(destination_root, inventory, direction)
  verification <- verify_lna_transfer(
    destination_root,
    error = TRUE,
    expected_inventory_hash = inventory$inventory_hash
  )
  structure(list(
    inventory = inventory,
    receipt = receipt,
    verification = verification,
    outcomes = outcomes,
    peak_buffer_bytes = peak
  ), class = c("lna_dataset_transfer_result", "list"))
}
