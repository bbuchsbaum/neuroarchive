# Append-only resume journal and single-coordinator advisory lock.

.lna_journal_schema_version <- "1.0.0"
.lna_lock_schema_version <- "1.0.0"
.lna_journal_states <- c(
  "started", "temporary-written", "completed", "failed", "invalidated",
  "lock-takeover", "journal-recovery", "plan-restart"
)

.journal_abort <- function(message, field, subclass = "lna_error_journal", ...) {
  abort_lna(
    message,
    ...,
    field = field,
    .subclass = subclass,
    location = paste0("lna_journal:", field)
  )
}

.journal_now <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
}

.journal_plan_id <- function(plan_hash) {
  paste0("urn:lna-plan:", sub("^sha256:", "", plan_hash))
}

.journal_paths <- function(path) {
  path <- .execution_scalar_string(path, "path")
  root <- if (basename(path) == ".lna_journal.jsonl") dirname(path) else path
  list(
    root = root,
    journal = file.path(root, ".lna_journal.jsonl"),
    lock = file.path(root, ".lna.lock"),
    lock_receipt = file.path(root, ".lna.lock", "owner.json"),
    lock_history = file.path(root, ".lna-lock-history"),
    journal_history = file.path(root, ".lna-journal-history"),
    recovery = file.path(root, ".lna-journal-recovery")
  )
}

.journal_file_checksum <- function(path) {
  paste0(
    "sha256:", digest::digest(file = path, algo = "sha256", serialize = FALSE)
  )
}

.journal_event_id <- function(event_without_id) {
  json <- jsonlite::toJSON(
    event_without_id, auto_unbox = TRUE, null = "null", digits = NA,
    matrix = "rowmajor", pretty = FALSE
  )
  paste0(
    "sha256:",
    digest::digest(charToRaw(enc2utf8(json)), algo = "sha256", serialize = FALSE)
  )
}

.journal_read_bytes <- function(path) {
  if (!file.exists(path)) {
    return(list(complete = raw(), fragment = raw(), file_size = 0))
  }
  size <- file.info(path)$size
  if (is.na(size) || size == 0) {
    return(list(complete = raw(), fragment = raw(), file_size = 0))
  }
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  bytes <- readBin(connection, what = "raw", n = size)
  newline <- as.raw(10L)
  if (identical(tail(bytes, 1L), newline)) {
    return(list(complete = bytes, fragment = raw(), file_size = length(bytes)))
  }
  positions <- which(bytes == newline)
  boundary <- if (length(positions)) tail(positions, 1L) else 0L
  complete <- if (boundary) bytes[seq_len(boundary)] else raw()
  fragment <- bytes[seq.int(boundary + 1L, length(bytes))]
  list(complete = complete, fragment = fragment, file_size = length(bytes))
}

.journal_lines_from_raw <- function(bytes) {
  if (!length(bytes)) return(character())
  text <- tryCatch(
    rawToChar(bytes),
    error = function(e) .journal_abort(
      "journal contains invalid byte sequences", "journal",
      subclass = "lna_error_journal_corrupt"
    )
  )
  lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
  if (length(lines) && !nzchar(tail(lines, 1L))) lines <- head(lines, -1L)
  if (any(!nzchar(lines))) {
    .journal_abort(
      "journal contains an empty complete line", "journal",
      subclass = "lna_error_journal_corrupt"
    )
  }
  lines
}

.journal_normalize_event <- function(event, line) {
  field <- paste0("line[", line, "]")
  if (!is.list(event) || is.data.frame(event)) {
    .journal_abort(
      paste0(field, " must be a JSON object"), field,
      subclass = "lna_error_journal_corrupt"
    )
  }
  required <- c(
    "journal_schema_version", "event_id", "sequence", "plan_id", "plan_hash",
    "node_id", "record_id", "source_checksum", "phase", "state", "attempt",
    "timestamp", "output_path", "output_checksum", "output_bytes", "failure",
    "writer_id", "details"
  )
  if (!setequal(names(event), required)) {
    .journal_abort(
      paste0(field, " must contain exactly the version-1 journal fields"), field,
      subclass = "lna_error_journal_corrupt"
    )
  }
  if (!identical(event$journal_schema_version, .lna_journal_schema_version)) {
    .journal_abort(
      paste0(field, " uses an unsupported journal schema version"), field,
      subclass = "lna_error_journal_corrupt"
    )
  }
  integer_scalar <- function(x, name, positive = FALSE) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != floor(x) ||
        x < if (positive) 1 else 0) {
      .journal_abort(
        paste0(field, ".", name, " must be an integer"), paste0(field, ".", name),
        subclass = "lna_error_journal_corrupt"
      )
    }
    as.double(x)
  }
  event$sequence <- as.integer(integer_scalar(
    event$sequence, "sequence", positive = TRUE
  ))
  event$attempt <- as.integer(integer_scalar(event$attempt, "attempt"))
  for (name in c("event_id", "plan_id", "plan_hash", "phase", "state", "timestamp", "writer_id")) {
    event[[name]] <- tryCatch(
      .execution_scalar_string(event[[name]], paste0(field, ".", name)),
      lna_error_execution = function(e) .journal_abort(
        conditionMessage(e), paste0(field, ".", name),
        subclass = "lna_error_journal_corrupt"
      )
    )
  }
  if (!event$state %in% .lna_journal_states) {
    .journal_abort(
      paste0(field, ".state is unsupported"), paste0(field, ".state"),
      subclass = "lna_error_journal_corrupt"
    )
  }
  if (!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}T", event$timestamp) ||
      !endsWith(event$timestamp, "Z")) {
    .journal_abort(
      paste0(field, ".timestamp is not UTC ISO-8601"), paste0(field, ".timestamp"),
      subclass = "lna_error_journal_corrupt"
    )
  }
  for (name in c("node_id", "record_id")) {
    if (!is.null(event[[name]])) {
      event[[name]] <- tryCatch(
        .execution_id(event[[name]], paste0(field, ".", name)),
        lna_error_execution = function(e) .journal_abort(
          conditionMessage(e), paste0(field, ".", name),
          subclass = "lna_error_journal_corrupt"
        )
      )
    }
  }
  for (name in c("plan_hash", "source_checksum", "output_checksum")) {
    if (!is.null(event[[name]])) {
      event[[name]] <- tryCatch(
        .source_record_checksum(event[[name]], paste0(field, ".", name)),
        lna_error_source_record = function(e) .journal_abort(
          conditionMessage(e), paste0(field, ".", name),
          subclass = "lna_error_journal_corrupt"
        )
      )
    }
  }
  if (!identical(event$plan_id, .journal_plan_id(event$plan_hash))) {
    .journal_abort(
      paste0(field, ".plan_id disagrees with plan_hash"), paste0(field, ".plan_id"),
      subclass = "lna_error_journal_corrupt"
    )
  }
  if (!is.null(event$output_path)) {
    event$output_path <- tryCatch(
      .manifest_validate_path(event$output_path, paste0(field, ".output_path")),
      lna_error_manifest = function(e) .journal_abort(
        conditionMessage(e), paste0(field, ".output_path"),
        subclass = "lna_error_journal_corrupt"
      )
    )
  }
  if (!is.null(event$output_bytes)) {
    event$output_bytes <- integer_scalar(event$output_bytes, "output_bytes")
  }
  for (name in c("failure", "details")) {
    tryCatch(
      .source_record_json_value(event[[name]], paste0(field, ".", name)),
      lna_error_source_record = function(e) .journal_abort(
        conditionMessage(e), paste0(field, ".", name),
        subclass = "lna_error_journal_corrupt"
      )
    )
  }
  supplied_id <- event$event_id
  payload <- event
  payload$event_id <- NULL
  expected_id <- .journal_event_id(payload)
  if (!identical(supplied_id, expected_id)) {
    .journal_abort(
      paste0(field, " event_id checksum mismatch"), paste0(field, ".event_id"),
      subclass = "lna_error_journal_corrupt"
    )
  }
  event
}

.journal_decode_complete <- function(bytes) {
  lines <- .journal_lines_from_raw(bytes)
  events <- lapply(seq_along(lines), function(i) {
    decoded <- tryCatch(
      jsonlite::fromJSON(lines[[i]], simplifyVector = FALSE),
      error = function(e) .journal_abort(
        paste0("journal line ", i, " is invalid JSON"), paste0("line[", i, "]"),
        subclass = "lna_error_journal_corrupt"
      )
    )
    .journal_normalize_event(decoded, i)
  })
  if (length(events)) {
    sequences <- vapply(events, `[[`, integer(1), "sequence")
    if (!identical(sequences, seq_along(events))) {
      .journal_abort(
        "journal sequence values are not contiguous", "sequence",
        subclass = "lna_error_journal_corrupt"
      )
    }
  }
  events
}

#' Read and verify an append-only LNA dataset journal
#'
#' @param path Archive root or `.lna_journal.jsonl` path.
#' @return A verified `lna_journal` object. A non-newline-terminated final
#'   fragment always fails and is never interpreted as an event.
#' @export
read_lna_journal <- function(path) {
  paths <- .journal_paths(path)
  bytes <- .journal_read_bytes(paths$journal)
  if (length(bytes$fragment)) {
    .journal_abort(
      "journal has a truncated final line; explicit recovery is required",
      "journal", subclass = "lna_error_journal_truncated",
      complete_bytes = length(bytes$complete),
      fragment_bytes = length(bytes$fragment)
    )
  }
  structure(
    list(path = paths$journal, events = .journal_decode_complete(bytes$complete)),
    class = c("lna_journal", "list")
  )
}

.journal_read_lock <- function(paths) {
  if (!dir.exists(paths$lock) || !file.exists(paths$lock_receipt)) return(NULL)
  tryCatch(
    jsonlite::read_json(paths$lock_receipt, simplifyVector = FALSE),
    error = function(e) NULL
  )
}

.journal_validate_lock_receipt <- function(receipt) {
  required <- c(
    "lock_schema_version", "owner_id", "pid", "hostname", "plan_hash",
    "acquired_at", "heartbeat_at"
  )
  if (!is.list(receipt) || !setequal(names(receipt), required) ||
      !identical(receipt$lock_schema_version, .lna_lock_schema_version)) {
    .journal_abort(
      "lock owner receipt is malformed", "lock",
      subclass = "lna_error_journal_locked"
    )
  }
  receipt$owner_id <- tryCatch(
    .execution_id(receipt$owner_id, "lock.owner_id"),
    lna_error_execution = function(e) .journal_abort(
      conditionMessage(e), "lock", subclass = "lna_error_journal_locked"
    )
  )
  receipt$plan_hash <- tryCatch(
    .source_record_checksum(receipt$plan_hash, "lock.plan_hash"),
    lna_error_source_record = function(e) .journal_abort(
      conditionMessage(e), "lock", subclass = "lna_error_journal_locked"
    )
  )
  for (name in c("hostname", "acquired_at", "heartbeat_at")) {
    if (!is.character(receipt[[name]]) || length(receipt[[name]]) != 1L ||
        is.na(receipt[[name]]) || !nzchar(receipt[[name]])) {
      .journal_abort(
        "lock owner receipt is malformed", "lock",
        subclass = "lna_error_journal_locked"
      )
    }
  }
  if (!is.numeric(receipt$pid) || length(receipt$pid) != 1L ||
      !is.finite(receipt$pid) || receipt$pid < 1 || receipt$pid != floor(receipt$pid)) {
    .journal_abort(
      "lock owner receipt is malformed", "lock",
      subclass = "lna_error_journal_locked"
    )
  }
  receipt
}

.journal_write_lock <- function(paths, receipt) {
  jsonlite::write_json(
    receipt, paths$lock_receipt, auto_unbox = TRUE, null = "null",
    digits = NA, pretty = TRUE
  )
  invisible(receipt)
}

.journal_timestamp_time <- function(x) {
  suppressWarnings(as.POSIXct(
    sub("Z$", "", x), format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC"
  ))
}

.journal_new_owner_id <- function() {
  seed <- paste(Sys.getpid(), .journal_now(), basename(tempfile("writer-")), sep = ":")
  paste0(
    "writer-", Sys.getpid(), "-",
    substr(digest::digest(seed, algo = "sha256", serialize = FALSE), 1L, 16L)
  )
}

.journal_acquire_lock <- function(
    paths,
    plan_hash,
    owner_id,
    stale_after,
    takeover_stale) {
  if (!dir.exists(paths$root) &&
      !dir.create(paths$root, recursive = TRUE, showWarnings = FALSE)) {
    .journal_abort("could not create archive root", "path")
  }
  acquired <- dir.create(paths$lock, recursive = FALSE, showWarnings = FALSE)
  takeover <- NULL
  if (!acquired) {
    existing <- .journal_read_lock(paths)
    if (is.null(existing)) {
      .journal_abort(
        "journal lock exists without a readable owner receipt", "lock",
        subclass = "lna_error_journal_locked"
      )
    }
    existing <- .journal_validate_lock_receipt(existing)
    heartbeat <- .journal_timestamp_time(existing$heartbeat_at %||% "")
    age <- if (is.na(heartbeat)) Inf else as.numeric(difftime(Sys.time(), heartbeat, units = "secs"))
    stale <- is.infinite(age) || age > stale_after
    if (!stale) {
      .journal_abort(
        paste0("archive is locked by live writer ", existing$owner_id), "lock",
        subclass = "lna_error_journal_locked", lock_receipt = existing
      )
    }
    if (!isTRUE(takeover_stale)) {
      .journal_abort(
        paste0("archive has stale lock from ", existing$owner_id,
          "; takeover_stale=TRUE is required"),
        "lock", subclass = "lna_error_journal_stale_lock",
        lock_receipt = existing
      )
    }
    dir.create(paths$lock_history, recursive = TRUE, showWarnings = FALSE)
    archive <- file.path(
      paths$lock_history,
      paste0(format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), "-", existing$owner_id)
    )
    if (!file.rename(paths$lock, archive)) {
      .journal_abort("stale lock changed during takeover", "lock")
    }
    if (!dir.create(paths$lock, recursive = FALSE, showWarnings = FALSE)) {
      .journal_abort("another writer won stale-lock takeover", "lock")
    }
    takeover <- list(previous_receipt = existing, archived_lock = basename(archive))
  }
  now <- .journal_now()
  receipt <- list(
    lock_schema_version = .lna_lock_schema_version,
    owner_id = owner_id,
    pid = Sys.getpid(),
    hostname = Sys.info()[["nodename"]] %||% "unknown",
    plan_hash = plan_hash,
    acquired_at = now,
    heartbeat_at = now
  )
  .journal_write_lock(paths, receipt)
  list(receipt = receipt, takeover = takeover)
}

.journal_release_lock_raw <- function(paths, owner_id) {
  receipt <- .journal_read_lock(paths)
  receipt <- tryCatch(
    if (is.null(receipt)) NULL else .journal_validate_lock_receipt(receipt),
    lna_error_journal_locked = function(e) NULL
  )
  if (is.null(receipt) || !identical(receipt$owner_id, owner_id)) return(FALSE)
  unlink(paths$lock, recursive = TRUE, force = FALSE)
  !dir.exists(paths$lock)
}

.journal_assert_handle <- function(handle) {
  if (!inherits(handle, "lna_journal_handle")) {
    .journal_abort("handle must be an lna_journal_handle", "handle")
  }
  receipt <- .journal_read_lock(handle$paths)
  if (!is.null(receipt)) receipt <- .journal_validate_lock_receipt(receipt)
  if (is.null(receipt) || !identical(receipt$owner_id, handle$owner_id) ||
      !identical(receipt$plan_hash, handle$plan$plan_hash)) {
    .journal_abort(
      "journal handle no longer owns the archive lock", "lock",
      subclass = "lna_error_journal_locked"
    )
  }
  invisible(receipt)
}

.journal_archive_truncated <- function(paths, bytes) {
  dir.create(paths$recovery, recursive = TRUE, showWarnings = FALSE)
  checksum <- paste0(
    "sha256:", digest::digest(bytes$fragment, algo = "sha256", serialize = FALSE)
  )
  archive_name <- paste0(
    format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), "-",
    substr(sub("^sha256:", "", checksum), 1L, 16L), ".fragment"
  )
  archive <- file.path(paths$recovery, archive_name)
  connection <- file(archive, open = "wb")
  on.exit(close(connection), add = TRUE)
  writeBin(bytes$fragment, connection)
  close(connection)
  on.exit(NULL, add = FALSE)
  journal_connection <- file(paths$journal, open = "wb")
  on.exit(close(journal_connection), add = TRUE)
  if (length(bytes$complete)) writeBin(bytes$complete, journal_connection)
  close(journal_connection)
  on.exit(NULL, add = FALSE)
  list(
    fragment_checksum = checksum,
    fragment_bytes = length(bytes$fragment),
    archive_path = .dataset_plan_join(".lna-journal-recovery", archive_name)
  )
}

.journal_existing_plan_hashes <- function(events) {
  unique(vapply(events, `[[`, character(1), "plan_hash"))
}

.journal_archive_for_restart <- function(paths, events) {
  dir.create(paths$journal_history, recursive = TRUE, showWarnings = FALSE)
  old_hash <- if (length(events)) events[[length(events)]]$plan_hash else "unknown"
  name <- paste0(
    format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), "-",
    substr(sub("^sha256:", "", old_hash), 1L, 16L), ".jsonl"
  )
  archive <- file.path(paths$journal_history, name)
  if (!file.rename(paths$journal, archive)) {
    .journal_abort("could not archive mismatched journal", "journal")
  }
  .dataset_plan_join(".lna-journal-history", name)
}

.journal_make_event <- function(
    handle,
    sequence,
    node_id,
    record_id,
    source_checksum,
    phase,
    state,
    attempt,
    output_path = NULL,
    output_checksum = NULL,
    output_bytes = NULL,
    failure = NULL,
    details = list()) {
  payload <- list(
    journal_schema_version = .lna_journal_schema_version,
    sequence = as.integer(sequence),
    plan_id = .journal_plan_id(handle$plan$plan_hash),
    plan_hash = handle$plan$plan_hash,
    node_id = node_id,
    record_id = record_id,
    source_checksum = source_checksum,
    phase = phase,
    state = state,
    attempt = as.integer(attempt),
    timestamp = .journal_now(),
    output_path = output_path,
    output_checksum = output_checksum,
    output_bytes = if (is.null(output_bytes)) NULL else as.double(output_bytes),
    failure = failure,
    writer_id = handle$owner_id,
    details = details
  )
  payload$event_id <- .journal_event_id(payload)
  payload <- payload[c(
    "journal_schema_version", "event_id", "sequence", "plan_id", "plan_hash",
    "node_id", "record_id", "source_checksum", "phase", "state", "attempt",
    "timestamp", "output_path", "output_checksum", "output_bytes", "failure",
    "writer_id", "details"
  )]
  payload
}

.journal_append_event <- function(handle, event) {
  .journal_assert_handle(handle)
  verified <- .journal_normalize_event(event, event$sequence)
  line <- jsonlite::toJSON(
    verified, auto_unbox = TRUE, null = "null", digits = NA,
    matrix = "rowmajor", pretty = FALSE
  )
  cat(line, "\n", file = handle$paths$journal, append = TRUE, sep = "")
  heartbeat_lna_journal(handle)
  invisible(verified)
}

.journal_events <- function(handle) {
  .journal_assert_handle(handle)
  journal <- read_lna_journal(handle$paths$journal)
  hashes <- .journal_existing_plan_hashes(journal$events)
  if (length(hashes) && !identical(hashes, handle$plan$plan_hash)) {
    .journal_abort("journal contains events for another plan", "plan_hash")
  }
  journal$events
}

.journal_append_admin <- function(handle, state, details) {
  events <- .journal_events(handle)
  event <- .journal_make_event(
    handle,
    sequence = length(events) + 1L,
    node_id = NULL,
    record_id = NULL,
    source_checksum = NULL,
    phase = "journal",
    state = state,
    attempt = 0L,
    details = details
  )
  .journal_append_event(handle, event)
}

#' Open a resumable LNA dataset journal under an advisory lock
#'
#' @param plan A validated `lna_dataset_plan`.
#' @param path Archive root.
#' @param owner_id Optional writer identity.
#' @param stale_after Seconds without a heartbeat before a lock is stale.
#' @param takeover_stale Explicitly take over a stale lock and audit it.
#' @param recover_truncated Either `"error"` or `"truncate"`. Truncation
#'   archives the rejected fragment and appends a recovery receipt.
#' @param plan_mismatch Either `"error"` or `"restart"`. Restart archives the
#'   old journal and starts a new audited stream.
#' @return An `lna_journal_handle` owning the archive lock.
#' @export
open_lna_journal <- function(
    plan,
    path,
    owner_id = NULL,
    stale_after = 300,
    takeover_stale = FALSE,
    recover_truncated = c("error", "truncate"),
    plan_mismatch = c("error", "restart")) {
  plan <- validate_lna_dataset_plan(plan)
  recover_truncated <- match.arg(recover_truncated)
  plan_mismatch <- match.arg(plan_mismatch)
  if (!is.numeric(stale_after) || length(stale_after) != 1L ||
      !is.finite(stale_after) || stale_after <= 0) {
    .journal_abort("stale_after must be positive seconds", "stale_after")
  }
  owner_id <- owner_id %||% .journal_new_owner_id()
  owner_id <- .execution_id(owner_id, "owner_id")
  paths <- .journal_paths(path)
  acquired <- .journal_acquire_lock(
    paths, plan$plan_hash, owner_id, stale_after, takeover_stale
  )
  released <- FALSE
  on.exit({
    if (!released) .journal_release_lock_raw(paths, owner_id)
  }, add = TRUE)
  recovery <- NULL
  bytes <- .journal_read_bytes(paths$journal)
  if (length(bytes$fragment)) {
    if (identical(recover_truncated, "error")) {
      .journal_abort(
        "journal has a truncated final line; recover_truncated='truncate' is required",
        "journal", subclass = "lna_error_journal_truncated"
      )
    }
    recovery <- .journal_archive_truncated(paths, bytes)
    bytes <- .journal_read_bytes(paths$journal)
  }
  events <- .journal_decode_complete(bytes$complete)
  hashes <- .journal_existing_plan_hashes(events)
  restart <- NULL
  if (length(hashes) && !identical(hashes, plan$plan_hash)) {
    if (identical(plan_mismatch, "error")) {
      .journal_abort(
        "existing journal belongs to a different plan; plan_mismatch='restart' is required",
        "plan_hash", existing_plan_hashes = hashes, requested_plan_hash = plan$plan_hash
      )
    }
    restart <- list(
      previous_plan_hashes = hashes,
      archived_journal = .journal_archive_for_restart(paths, events)
    )
  }
  handle <- structure(
    list(plan = plan, paths = paths, owner_id = owner_id),
    class = c("lna_journal_handle", "list")
  )
  if (!is.null(restart)) .journal_append_admin(handle, "plan-restart", restart)
  if (!is.null(recovery)) .journal_append_admin(handle, "journal-recovery", recovery)
  if (!is.null(acquired$takeover)) {
    .journal_append_admin(handle, "lock-takeover", acquired$takeover)
  }
  released <- TRUE
  handle
}

#' Maintain or release an LNA journal lock
#'
#' @param handle An open `lna_journal_handle`.
#' @return The handle invisibly for heartbeat; a logical release result for
#'   close.
#' @export
heartbeat_lna_journal <- function(handle) {
  receipt <- .journal_assert_handle(handle)
  receipt$heartbeat_at <- .journal_now()
  .journal_write_lock(handle$paths, receipt)
  invisible(handle)
}

#' @rdname heartbeat_lna_journal
#' @export
close_lna_journal <- function(handle) {
  if (!inherits(handle, "lna_journal_handle")) {
    .journal_abort("handle must be an lna_journal_handle", "handle")
  }
  invisible(.journal_release_lock_raw(handle$paths, handle$owner_id))
}

.journal_failure <- function(failure) {
  if (is.null(failure)) return(NULL)
  if (inherits(failure, "condition")) {
    return(list(
      message = conditionMessage(failure),
      class = class(failure)
    ))
  }
  if (is.character(failure) && length(failure) == 1L && !is.na(failure)) {
    return(list(message = failure, class = "character"))
  }
  tryCatch(
    .source_record_json_value(failure, "failure"),
    lna_error_source_record = function(e) .journal_abort(conditionMessage(e), "failure")
  )
  failure
}

.journal_node_map <- function(handle) {
  .execution_nodes_by_id(handle$plan)
}

.journal_latest_node_event <- function(events, node_id) {
  selected <- events[vapply(events, function(event) {
    identical(event$node_id, node_id)
  }, logical(1))]
  if (!length(selected)) return(NULL)
  selected[[length(selected)]]
}

.journal_verify_artifact_event <- function(handle, node, event) {
  if (!identical(event$record_id, node$record_id) ||
      !identical(event$source_checksum, node$source_checksum) ||
      !identical(event$output_path, node$artifact_path)) {
    return("journal identity does not match the planned artifact")
  }
  if (is.null(event$output_checksum) || is.null(event$output_bytes)) {
    return("completed artifact event lacks checksum or byte size")
  }
  path <- file.path(
    handle$paths$root,
    do.call(file.path, as.list(strsplit(event$output_path, "/", fixed = TRUE)[[1L]]))
  )
  if (!file.exists(path)) return("completed artifact is missing")
  size <- file.info(path)$size
  if (is.na(size) || !identical(as.double(size), as.double(event$output_bytes))) {
    return("completed artifact byte size changed")
  }
  checksum <- .journal_file_checksum(path)
  if (!identical(checksum, event$output_checksum)) {
    return("completed artifact checksum changed")
  }
  NULL
}

.journal_resume_core <- function(handle, invalid, audit) {
  .journal_assert_handle(handle)
  invalid <- match.arg(invalid, c("rebuild", "error"))
  events <- .journal_events(handle)
  nodes <- .journal_node_map(handle)
  completed <- character()
  invalidations <- list()
  for (node_id in handle$plan$execution$topological_order) {
    node <- nodes[[node_id]]
    event <- .journal_latest_node_event(events, node_id)
    if (is.null(event) || !identical(event$state, "completed")) next
    if (!all(node$depends_on %in% completed)) next
    reason <- NULL
    if (node$kind %in% c("artifact-write", "artifact-validation")) {
      reason <- .journal_verify_artifact_event(handle, node, event)
    }
    if (!is.null(reason)) {
      if (identical(invalid, "error")) {
        .journal_abort(
          paste0("cannot reuse ", node_id, ": ", reason),
          "resume", node_id = node_id, reason = reason
        )
      }
      invalidations[[node_id]] <- reason
      if (isTRUE(audit)) {
        append_lna_journal_event(
          handle, node_id, "invalidated", attempt = event$attempt,
          details = list(reason = reason, invalidated_event_id = event$event_id)
        )
        events <- .journal_events(handle)
      }
      next
    }
    completed <- c(completed, node_id)
  }
  state <- lna_execution_state(handle$plan, completed)
  write_ids <- paste0("write:", vapply(handle$plan$artifacts, `[[`, character(1), "record_id"))
  reusable <- sub("^write:", "", intersect(write_ids, completed))
  all_records <- vapply(handle$plan$artifacts, `[[`, character(1), "record_id")
  structure(
    list(
      state = state,
      reusable_record_ids = reusable,
      rebuild_record_ids = setdiff(all_records, reusable),
      invalidations = invalidations,
      events = events
    ),
    class = c("lna_journal_resume", "list")
  )
}

#' Resume execution from verified journal events
#'
#' @param handle An open `lna_journal_handle`.
#' @param invalid Whether invalid completed outputs are audited for rebuild or
#'   cause an error.
#' @return An `lna_journal_resume` with dependency-closed execution state and
#'   reusable/rebuild record IDs.
#' @export
resume_lna_journal <- function(handle, invalid = c("rebuild", "error")) {
  .journal_resume_core(handle, match.arg(invalid), audit = TRUE)
}

#' Append a state transition to an LNA dataset journal
#'
#' @param handle An open `lna_journal_handle`.
#' @param node_id Execution node identity.
#' @param state One of started, temporary-written, completed, failed, or
#'   invalidated.
#' @param attempt Attempt number. It is allocated automatically for started
#'   events and inferred for later states.
#' @param output_path Optional archive-relative output path.
#' @param output_checksum Optional checksum, verified for completed artifacts.
#' @param output_bytes Optional byte size, verified for completed artifacts.
#' @param failure Failure condition, message, or JSON object.
#' @param details Additional JSON-safe event details.
#' @return The appended event, invisibly.
#' @export
append_lna_journal_event <- function(
    handle,
    node_id,
    state = c("started", "temporary-written", "completed", "failed", "invalidated"),
    attempt = NULL,
    output_path = NULL,
    output_checksum = NULL,
    output_bytes = NULL,
    failure = NULL,
    details = list()) {
  .journal_assert_handle(handle)
  state <- match.arg(state)
  node_id <- .execution_id(node_id, "node_id")
  nodes <- .journal_node_map(handle)
  if (!node_id %in% names(nodes)) {
    .journal_abort(paste0("unknown execution node: ", node_id), "node_id")
  }
  node <- nodes[[node_id]]
  events <- .journal_events(handle)
  node_events <- events[vapply(events, function(event) {
    identical(event$node_id, node_id)
  }, logical(1))]
  started_attempts <- vapply(
    node_events[vapply(node_events, function(event) event$state == "started", logical(1))],
    `[[`, integer(1), "attempt"
  )
  if (identical(state, "started")) {
    resume <- .journal_resume_core(handle, "error", audit = FALSE)
    missing <- setdiff(node$depends_on, resume$state$completed_nodes)
    if (length(missing)) {
      .journal_abort(
        paste0("cannot start ", node_id, "; incomplete dependencies: ",
          paste(missing, collapse = ", ")),
        "node_id", incomplete_dependencies = missing
      )
    }
    attempt <- if (length(started_attempts)) max(started_attempts) + 1L else 1L
  } else {
    if (is.null(attempt)) {
      if (!length(started_attempts)) {
        .journal_abort("a non-started event requires a started attempt", "attempt")
      }
      attempt <- max(started_attempts)
    }
    if (!attempt %in% started_attempts && !identical(state, "invalidated")) {
      .journal_abort("event attempt has no matching started event", "attempt")
    }
  }
  attempt <- as.integer(attempt)
  latest_same_attempt <- node_events[vapply(node_events, function(event) {
    identical(event$attempt, attempt)
  }, logical(1))]
  latest_state <- if (length(latest_same_attempt)) {
    latest_same_attempt[[length(latest_same_attempt)]]$state
  } else {
    NULL
  }
  allowed <- list(
    started = c("temporary-written", "completed", "failed", "invalidated"),
    `temporary-written` = c("completed", "failed", "invalidated"),
    completed = "invalidated",
    failed = character(),
    invalidated = character()
  )
  if (!identical(state, "started") &&
      (is.null(latest_state) || !state %in% allowed[[latest_state]])) {
    .journal_abort(
      paste0("invalid journal transition from ", latest_state %||% "none", " to ", state),
      "state"
    )
  }
  if (identical(state, "failed") && is.null(failure)) {
    .journal_abort("failed events require failure details", "failure")
  }
  output_path <- output_path %||% if (
    node$kind %in% c("artifact-write", "artifact-validation") &&
      identical(state, "completed")
  ) node$artifact_path else NULL
  if (!is.null(output_path)) output_path <- .manifest_validate_path(output_path, "output_path")
  if (identical(state, "completed") &&
      node$kind %in% c("artifact-write", "artifact-validation")) {
    if (!identical(output_path, node$artifact_path)) {
      .journal_abort("completed output path disagrees with plan", "output_path")
    }
    full_path <- file.path(
      handle$paths$root,
      do.call(file.path, as.list(strsplit(output_path, "/", fixed = TRUE)[[1L]]))
    )
    if (!file.exists(full_path)) {
      .journal_abort("completed output does not exist", "output_path")
    }
    actual_checksum <- .journal_file_checksum(full_path)
    actual_bytes <- as.double(file.info(full_path)$size)
    if (!is.null(output_checksum) && !identical(output_checksum, actual_checksum)) {
      .journal_abort("supplied output checksum does not match file", "output_checksum")
    }
    if (!is.null(output_bytes) && !identical(as.double(output_bytes), actual_bytes)) {
      .journal_abort("supplied output byte size does not match file", "output_bytes")
    }
    output_checksum <- actual_checksum
    output_bytes <- actual_bytes
  }
  tryCatch(
    .source_record_json_value(details, "details"),
    lna_error_source_record = function(e) .journal_abort(conditionMessage(e), "details")
  )
  event <- .journal_make_event(
    handle,
    sequence = length(events) + 1L,
    node_id = node_id,
    record_id = node$record_id,
    source_checksum = node$source_checksum,
    phase = node$phase,
    state = state,
    attempt = attempt,
    output_path = output_path,
    output_checksum = output_checksum,
    output_bytes = output_bytes,
    failure = .journal_failure(failure),
    details = details
  )
  .journal_append_event(handle, event)
}

#' Derive the final manifest from a plan and verified journal
#'
#' Existing manifest files are never consulted. Every artifact write and
#' validation node must be complete and checksum-valid.
#'
#' @param handle An open `lna_journal_handle`.
#' @return A complete `lna_manifest` derived from the plan and journal.
#' @export
lna_manifest_from_journal <- function(handle) {
  resume <- .journal_resume_core(handle, "error", audit = FALSE)
  record_ids <- vapply(handle$plan$artifacts, `[[`, character(1), "record_id")
  required <- c(paste0("write:", record_ids), paste0("validate:", record_ids))
  missing <- setdiff(required, resume$state$completed_nodes)
  if (length(missing)) {
    .journal_abort(
      paste0("cannot derive final manifest; incomplete node(s): ",
        paste(missing, collapse = ", ")),
      "manifest"
    )
  }
  states <- stats::setNames(lapply(record_ids, function(record_id) {
    event <- .journal_latest_node_event(resume$events, paste0("write:", record_id))
    list(
      lna_checksum = event$output_checksum,
      lna_bytes = event$output_bytes,
      status = "complete"
    )
  }), record_ids)
  lna_manifest(handle$plan, artifact_state = states)
}
