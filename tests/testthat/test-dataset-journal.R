journal_test_checksum <- function(id) {
  paste0(
    "sha256:",
    digest::digest(charToRaw(id), algo = "sha256", serialize = FALSE)
  )
}

journal_test_record <- function(record_id, sub, checksum_seed = record_id) {
  lna_source_record(
    record_id = record_id,
    source_path = paste0("sub-", sub, "/func/", record_id, ".nii.gz"),
    source_checksum = journal_test_checksum(checksum_seed),
    entities = list(sub = sub, task = "rest", run = "01"),
    datatype = "func",
    suffix = "bold",
    extension = ".nii.gz",
    resolved_metadata = list(RepetitionTime = 0.8),
    geometry = list(
      dims = c(3L, 3L, 2L),
      affine = diag(c(2, 2, 2, 1)),
      units = "mm",
      voxel_ordering = "first-axis-fastest",
      n_volumes = 6L
    )
  )
}

journal_test_plan <- function(root, records = NULL) {
  if (is.null(records)) {
    records <- list(
      journal_test_record("record-a", "01"),
      journal_test_record("record-b", "02")
    )
  }
  plan_lna_dataset(records, root)
}

journal_test_artifact_path <- function(root, plan, record_id) {
  entry <- plan$artifacts[[match(
    record_id,
    vapply(plan$artifacts, `[[`, character(1), "record_id")
  )]]
  file.path(root, entry$destination_relative)
}

journal_test_write_artifact <- function(root, plan, record_id, content = record_id) {
  path <- journal_test_artifact_path(root, plan, record_id)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(content, path, useBytes = TRUE)
  path
}

journal_test_complete_node <- function(handle, node_id) {
  append_lna_journal_event(handle, node_id, "started")
  append_lna_journal_event(handle, node_id, "completed")
}

journal_test_complete_prelude <- function(handle) {
  journal_test_complete_node(handle, "phase:preflight")
  journal_test_complete_node(handle, "phase:fit-complete")
}

journal_test_complete_write <- function(handle, root, record_id, content = record_id) {
  journal_test_write_artifact(root, handle$plan, record_id, content)
  node_id <- paste0("write:", record_id)
  append_lna_journal_event(handle, node_id, "started")
  append_lna_journal_event(handle, node_id, "completed")
}

test_that("started, temporary, and completed states resume conservatively", {
  root <- withr::local_tempdir()
  plan <- journal_test_plan(root)
  handle <- open_lna_journal(plan, root, owner_id = "writer-first")
  withr::defer(close_lna_journal(handle))
  journal_test_complete_prelude(handle)

  append_lna_journal_event(handle, "write:record-a", "started")
  append_lna_journal_event(
    handle, "write:record-a", "temporary-written",
    output_path = ".lna-tmp/record-a.partial"
  )
  before_crash <- resume_lna_journal(handle)
  expect_identical(before_crash$reusable_record_ids, character())
  expect_setequal(before_crash$rebuild_record_ids, c("record-a", "record-b"))
  expect_true(close_lna_journal(handle))

  handle <- open_lna_journal(plan, root, owner_id = "writer-resume")
  journal_test_write_artifact(root, plan, "record-a", "completed-a")
  started <- append_lna_journal_event(handle, "write:record-a", "started")
  completed <- append_lna_journal_event(handle, "write:record-a", "completed")
  resumed <- resume_lna_journal(handle)

  expect_identical(started$attempt, 2L)
  expect_identical(completed$plan_id, paste0(
    "urn:lna-plan:", sub("^sha256:", "", plan$plan_hash)
  ))
  expect_identical(completed$plan_hash, plan$plan_hash)
  expect_identical(completed$record_id, "record-a")
  expect_identical(
    completed$source_checksum,
    plan$artifacts[[1L]]$source$checksum
  )
  expect_identical(completed$phase, "write")
  expect_match(completed$timestamp, "Z$")
  expect_match(completed$output_checksum, "^sha256:[0-9a-f]{64}$")
  expect_true(completed$output_bytes > 0)
  expect_identical(resumed$reusable_record_ids, "record-a")

  journal <- read_lna_journal(root)
  expect_identical(
    vapply(journal$events, `[[`, integer(1), "sequence"),
    seq_along(journal$events)
  )
  expect_true(all(vapply(journal$events, function(event) {
    identical(event$plan_hash, plan$plan_hash)
  }, logical(1))))
})

test_that("failure details are durable and do not complete a node", {
  root <- withr::local_tempdir()
  plan <- journal_test_plan(root, list(journal_test_record("record-a", "01")))
  handle <- open_lna_journal(plan, root, owner_id = "writer-failure")
  withr::defer(close_lna_journal(handle))
  append_lna_journal_event(handle, "phase:preflight", "started")
  append_lna_journal_event(
    handle, "phase:preflight", "failed",
    failure = simpleError("inventory failed"),
    details = list(stage = "inventory")
  )

  journal <- read_lna_journal(root)
  failure <- tail(journal$events, 1L)[[1L]]
  expect_identical(failure$state, "failed")
  expect_identical(failure$failure$message, "inventory failed")
  expect_true("simpleError" %in% unlist(failure$failure$class))
  expect_identical(failure$details$stage, "inventory")
  expect_identical(resume_lna_journal(handle)$state$completed_nodes, character())

  restarted <- append_lna_journal_event(handle, "phase:preflight", "started")
  expect_identical(restarted$attempt, 2L)
})

test_that("final manifest is derived from verified plan and journal", {
  root <- withr::local_tempdir()
  plan <- journal_test_plan(root)
  handle <- open_lna_journal(plan, root, owner_id = "writer-manifest")
  withr::defer(close_lna_journal(handle))
  journal_test_complete_prelude(handle)
  journal_test_complete_write(handle, root, "record-a", "artifact-a")
  journal_test_complete_write(handle, root, "record-b", "artifact-b")
  journal_test_complete_node(handle, "phase:write-complete")
  journal_test_complete_node(handle, "validate:record-a")
  journal_test_complete_node(handle, "validate:record-b")

  writeLines("not a manifest", file.path(root, "lna_manifest.tsv"))
  writeLines("not json", file.path(root, "lna_manifest.json"))
  manifest <- lna_manifest_from_journal(handle)

  expect_s3_class(manifest, "lna_manifest")
  expect_true(all(manifest$data$status == "complete"))
  expect_setequal(manifest$data$record_id, c("record-a", "record-b"))
  for (record_id in manifest$data$record_id) {
    row <- manifest$data[manifest$data$record_id == record_id, , drop = FALSE]
    expect_identical(
      row$lna_checksum,
      neuroarchive:::.journal_file_checksum(
        journal_test_artifact_path(root, plan, record_id)
      )
    )
    expect_identical(
      row$lna_bytes,
      as.double(file.info(journal_test_artifact_path(root, plan, record_id))$size)
    )
  }
  expect_false("phase:finalize" %in% resume_lna_journal(handle)$state$completed_nodes)
})

test_that("changed source or plan fails unless restart is explicit", {
  root <- withr::local_tempdir()
  old_record <- journal_test_record("record-a", "01", checksum_seed = "old")
  new_record <- journal_test_record("record-a", "01", checksum_seed = "new")
  old_plan <- journal_test_plan(root, list(old_record))
  new_plan <- journal_test_plan(root, list(new_record))
  old_handle <- open_lna_journal(old_plan, root, owner_id = "writer-old")
  journal_test_complete_node(old_handle, "phase:preflight")
  expect_true(close_lna_journal(old_handle))

  expect_error(
    open_lna_journal(new_plan, root, owner_id = "writer-new"),
    "different plan",
    class = "lna_error_journal"
  )
  new_handle <- open_lna_journal(
    new_plan, root, owner_id = "writer-new",
    plan_mismatch = "restart"
  )
  withr::defer(close_lna_journal(new_handle))
  journal <- read_lna_journal(root)

  expect_length(journal$events, 1L)
  expect_identical(journal$events[[1L]]$state, "plan-restart")
  expect_identical(journal$events[[1L]]$plan_hash, new_plan$plan_hash)
  expect_true(dir.exists(file.path(root, ".lna-journal-history")))
  expect_length(
    list.files(file.path(root, ".lna-journal-history"), pattern = "[.]jsonl$"),
    1L
  )
})

test_that("missing and corrupt completed outputs are invalidated for rebuild", {
  root <- withr::local_tempdir()
  plan <- journal_test_plan(root)
  handle <- open_lna_journal(plan, root, owner_id = "writer-output")
  withr::defer(close_lna_journal(handle))
  journal_test_complete_prelude(handle)
  path_a <- journal_test_write_artifact(root, plan, "record-a", "artifact-a")
  path_b <- journal_test_write_artifact(root, plan, "record-b", "artifact-b")
  for (record_id in c("record-a", "record-b")) {
    append_lna_journal_event(handle, paste0("write:", record_id), "started")
    append_lna_journal_event(handle, paste0("write:", record_id), "completed")
  }
  unlink(path_a)
  writeLines("artifact-B-corrupted", path_b, useBytes = TRUE)

  expect_error(
    resume_lna_journal(handle, invalid = "error"),
    "cannot reuse write:record-a",
    class = "lna_error_journal"
  )
  resumed <- resume_lna_journal(handle, invalid = "rebuild")

  expect_identical(resumed$reusable_record_ids, character())
  expect_setequal(resumed$rebuild_record_ids, c("record-a", "record-b"))
  expect_match(resumed$invalidations[["write:record-a"]], "missing")
  expect_match(resumed$invalidations[["write:record-b"]], "byte size|checksum")
  journal <- read_lna_journal(root)
  invalidated <- journal$events[vapply(journal$events, function(event) {
    event$state == "invalidated"
  }, logical(1))]
  expect_setequal(
    vapply(invalidated, `[[`, character(1), "node_id"),
    c("write:record-a", "write:record-b")
  )
})

test_that("live second writer is rejected and stale takeover is audited", {
  root <- withr::local_tempdir()
  plan <- journal_test_plan(root, list(journal_test_record("record-a", "01")))
  first <- open_lna_journal(plan, root, owner_id = "writer-live", stale_after = 60)
  withr::defer(close_lna_journal(first))

  error <- expect_error(
    open_lna_journal(plan, root, owner_id = "writer-second", stale_after = 60),
    class = "lna_error_journal_locked"
  )
  expect_match(conditionMessage(error), "live writer writer-live")
  heartbeat_lna_journal(first)

  receipt_path <- file.path(root, ".lna.lock", "owner.json")
  receipt <- jsonlite::read_json(receipt_path, simplifyVector = FALSE)
  receipt$heartbeat_at <- "2000-01-01T00:00:00.000Z"
  jsonlite::write_json(
    receipt, receipt_path, auto_unbox = TRUE, null = "null", pretty = TRUE
  )
  expect_error(
    open_lna_journal(plan, root, owner_id = "writer-takeover", stale_after = 1),
    "takeover_stale=TRUE",
    class = "lna_error_journal_stale_lock"
  )
  takeover <- open_lna_journal(
    plan, root, owner_id = "writer-takeover", stale_after = 1,
    takeover_stale = TRUE
  )
  withr::defer(close_lna_journal(takeover))

  expect_error(
    append_lna_journal_event(first, "phase:preflight", "started"),
    "no longer owns",
    class = "lna_error_journal_locked"
  )
  journal <- read_lna_journal(root)
  expect_identical(tail(journal$events, 1L)[[1L]]$state, "lock-takeover")
  expect_identical(
    tail(journal$events, 1L)[[1L]]$details$previous_receipt$owner_id,
    "writer-live"
  )
  expect_true(dir.exists(file.path(root, ".lna-lock-history")))
  expect_length(list.files(file.path(root, ".lna-lock-history")), 1L)
})

test_that("truncated tail is archived and never accepted as an event", {
  root <- withr::local_tempdir()
  plan <- journal_test_plan(root, list(journal_test_record("record-a", "01")))
  handle <- open_lna_journal(plan, root, owner_id = "writer-truncate")
  append_lna_journal_event(handle, "phase:preflight", "started")
  close_lna_journal(handle)
  journal_path <- file.path(root, ".lna_journal.jsonl")
  connection <- file(journal_path, open = "ab")
  writeBin(charToRaw('{"journal_schema_version":"1.0.0"'), connection)
  close(connection)

  error <- expect_error(
    read_lna_journal(root), class = "lna_error_journal_truncated"
  )
  expect_true(error$fragment_bytes > 0)
  expect_error(
    open_lna_journal(plan, root, owner_id = "writer-no-recovery"),
    "recover_truncated='truncate'",
    class = "lna_error_journal_truncated"
  )
  recovered <- open_lna_journal(
    plan, root, owner_id = "writer-recovery", recover_truncated = "truncate"
  )
  withr::defer(close_lna_journal(recovered))
  journal <- read_lna_journal(root)

  expect_identical(
    vapply(journal$events, `[[`, character(1), "state"),
    c("started", "journal-recovery")
  )
  expect_identical(
    vapply(journal$events, `[[`, integer(1), "sequence"),
    c(1L, 2L)
  )
  expect_true(dir.exists(file.path(root, ".lna-journal-recovery")))
  expect_length(
    list.files(file.path(root, ".lna-journal-recovery"), pattern = "[.]fragment$"),
    1L
  )
})

test_that("a corrupt complete journal line fails hard", {
  root <- withr::local_tempdir()
  plan <- journal_test_plan(root, list(journal_test_record("record-a", "01")))
  handle <- open_lna_journal(plan, root, owner_id = "writer-corrupt")
  append_lna_journal_event(handle, "phase:preflight", "started")
  close_lna_journal(handle)
  cat("not-json\n", file = file.path(root, ".lna_journal.jsonl"), append = TRUE)

  expect_error(
    read_lna_journal(root), "invalid JSON",
    class = "lna_error_journal_corrupt"
  )
  expect_error(
    open_lna_journal(
      plan, root, owner_id = "writer-cannot-truncate",
      recover_truncated = "truncate"
    ),
    "invalid JSON",
    class = "lna_error_journal_corrupt"
  )
  expect_false(dir.exists(file.path(root, ".lna.lock")))
})
