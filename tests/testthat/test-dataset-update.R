update_test_record <- function(source_root, record_id, sub, value = "v1") {
  relative <- paste0("sub-", sub, "/func/", record_id, ".source")
  path <- file.path(source_root, relative)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(paste(record_id, value, sep = ":"), path, useBytes = TRUE)
  metadata_relative <- paste0("sub-", sub, "/func/", record_id, ".json")
  metadata_path <- file.path(source_root, metadata_relative)
  jsonlite::write_json(
    list(RepetitionTime = 0.8), metadata_path, auto_unbox = TRUE, pretty = TRUE
  )
  lna_source_record(
    record_id = record_id,
    source_path = relative,
    source_checksum = neuroarchive:::.journal_file_checksum(path),
    entities = list(sub = sub, task = "rest", run = "01"),
    datatype = "func", suffix = "bold", extension = ".nii.gz",
    resolved_metadata = list(RepetitionTime = 0.8),
    associations = list(list(
      association_id = paste0("metadata-", record_id),
      role = "metadata-source", path = metadata_relative,
      checksum = neuroarchive:::.journal_file_checksum(metadata_path),
      metadata = list(
        inheritance_level = "run", precedence = 100L,
        fields = "RepetitionTime", media_type = "application/json"
      )
    )),
    geometry = list(
      dims = c(2L, 2L, 2L), affine = diag(c(2, 2, 2, 1)),
      units = "mm", voxel_ordering = "first-axis-fastest", n_volumes = 3L
    )
  )
}

update_test_hooks <- function(state = new.env(parent = emptyenv())) {
  state$writes <- state$writes %||% character()
  list(
    load_source = function(entry, source_path) {
      paste(readLines(source_path, warn = FALSE), collapse = "")
    },
    write_artifact = function(data, temp_path, entry, plan) {
      state$writes <- c(state$writes, entry$record_id)
      writeLines(paste0("LNA:", entry$record_id, ":", data), temp_path, useBytes = TRUE)
    },
    validate_artifact = function(path, entry, stage) {
      if (!startsWith(readLines(path, warn = FALSE)[[1L]], "LNA:")) {
        stop("not a test LNA")
      }
      TRUE
    },
    write_manifest = function(manifest, staging_root) {
      if (isTRUE(state$fail_manifest)) {
        state$fail_manifest <- FALSE
        stop("injected update finalization failure")
      }
      write_lna_manifest(manifest, staging_root)
    }
  )
}

update_test_archive <- function(source_root, archive_root) {
  records <- list(
    update_test_record(source_root, "record-a", "01"),
    update_test_record(source_root, "record-b", "02")
  )
  write_lna_dataset(
    records, archive_root, source_root,
    hooks = update_test_hooks(), owner_id = "update-initial"
  )
  records
}

update_statuses <- function(plan) {
  stats::setNames(
    vapply(plan$classifications, `[[`, character(1), "status"),
    vapply(plan$classifications, `[[`, character(1), "record_id")
  )
}

update_test_hard_kill <- function(plan, hooks, checkpoint) {
  testthat::skip_on_os("windows")
  bundle <- tempfile("lna-update-kill-bundle-", fileext = ".rds")
  script <- tempfile("lna-update-kill-", fileext = ".R")
  log <- tempfile("lna-update-kill-", fileext = ".log")
  saveRDS(list(plan = plan, hooks = hooks), bundle)
  writeLines(c(
    "args <- commandArgs(trailingOnly = TRUE)",
    "package_root <- args[[1L]]",
    "if (file.exists(file.path(package_root, 'R', 'dataset_update.R'))) {",
    "  pkgload::load_all(package_root, quiet = TRUE)",
    "} else {",
    "  .libPaths(c(dirname(package_root), .libPaths()))",
    "  suppressPackageStartupMessages(library('neuroarchive', character.only = TRUE))",
    "}",
    "bundle <- readRDS(args[[2L]])",
    "checkpoint <- args[[3L]]",
    "bundle$hooks$after_update_checkpoint <- function(phase, state) {",
    "  if (identical(phase, checkpoint)) tools::pskill(Sys.getpid(), signal = 9L)",
    "  invisible(TRUE)",
    "}",
    "neuroarchive::update_lna_dataset(bundle$plan, hooks = bundle$hooks)"
  ), script, useBytes = TRUE)
  package_root <- normalizePath(
    find.package("neuroarchive"), winslash = "/", mustWork = TRUE
  )
  status <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(script), shQuote(package_root), shQuote(bundle),
      shQuote(checkpoint)),
    stdout = log, stderr = log
  ))
  expect_false(identical(as.integer(status), 0L), info = checkpoint)
  invisible(log)
}

test_that("update preview classifies no-op, add, change, and missing records", {
  source_root <- withr::local_tempdir()
  archive_root <- tempfile("lna-update-classify-")
  records <- update_test_archive(source_root, archive_root)

  no_op <- plan_lna_dataset_update(
    archive_root, records, source_root, mode = "refresh"
  )
  expect_true(all(update_statuses(no_op) == "unchanged"))

  changed <- update_test_record(source_root, "record-a", "01", "v2")
  added <- update_test_record(source_root, "record-c", "03")
  preview <- plan_lna_dataset_update(
    archive_root, list(changed, added), source_root, mode = "append"
  )
  expect_identical(unname(update_statuses(preview)[c("record-a", "record-b", "record-c")]),
    c("changed", "source-missing", "added"))
  actions <- stats::setNames(
    vapply(preview$classifications, `[[`, character(1), "action"),
    vapply(preview$classifications, `[[`, character(1), "record_id")
  )
  expect_identical(unname(actions[c("record-a", "record-b", "record-c")]),
    c("retain", "retain", "add"))
})

test_that("append writes only added records and records a new revision", {
  source_root <- withr::local_tempdir()
  archive_root <- tempfile("lna-update-append-")
  records <- update_test_archive(source_root, archive_root)
  old_manifest <- read_lna_manifest(archive_root)$data
  old_checksums <- stats::setNames(old_manifest$lna_checksum, old_manifest$record_id)
  changed <- update_test_record(source_root, "record-a", "01", "v2")
  added <- update_test_record(source_root, "record-c", "03")
  preview <- plan_lna_dataset_update(
    archive_root, list(changed, records[[2L]], added), source_root, mode = "append"
  )
  state <- new.env(parent = emptyenv())
  result <- update_lna_dataset(preview, hooks = update_test_hooks(state))

  expect_identical(state$writes, "record-c")
  manifest <- read_lna_manifest(archive_root)$data
  checksums <- stats::setNames(manifest$lna_checksum, manifest$record_id)
  expect_identical(checksums[["record-a"]], old_checksums[["record-a"]])
  expect_setequal(manifest$record_id, c("record-a", "record-b", "record-c"))
  expect_identical(as.integer(result$receipt$revision), 2L)
  expect_true(dir.exists(result$previous_generation))
  expect_identical(
    read_lna_dataset_update(archive_root)$archive_plan_hash,
    read_lna_dataset_declaration(archive_root)$plan_hash
  )
  events <- readLines(file.path(archive_root, ".lna_update_journal.jsonl"), warn = FALSE)
  expect_true(any(grepl('"state":"completed"', events, fixed = TRUE)))
  repeated <- update_lna_dataset(preview, hooks = update_test_hooks())
  expect_true(repeated$recovered)
  expect_identical(repeated$receipt$update_plan_hash, preview$update_plan_hash)
})

test_that("refresh rebuilds changed and corrupt records but retains unchanged ones", {
  source_root <- withr::local_tempdir()
  archive_root <- tempfile("lna-update-refresh-")
  records <- update_test_archive(source_root, archive_root)
  records[[1L]] <- update_test_record(source_root, "record-a", "01", "v2")
  preview <- plan_lna_dataset_update(
    archive_root, records, source_root, mode = "refresh"
  )
  expect_identical(update_statuses(preview)[["record-a"]], "changed")
  state <- new.env(parent = emptyenv())
  update_lna_dataset(preview, hooks = update_test_hooks(state))
  expect_identical(state$writes, "record-a")

  manifest <- read_lna_manifest(archive_root)$data
  corrupt_path <- file.path(
    archive_root, manifest$lna_path[manifest$record_id == "record-b"]
  )
  writeLines("CORRUPT", corrupt_path, useBytes = TRUE)
  corrupt <- plan_lna_dataset_update(
    archive_root, records, source_root, mode = "refresh"
  )
  expect_identical(update_statuses(corrupt)[["record-b"]], "corrupt")
  repair_state <- new.env(parent = emptyenv())
  update_lna_dataset(corrupt, hooks = update_test_hooks(repair_state))
  expect_identical(repair_state$writes, "record-b")
})

test_that("mirror pruning requires an explicit preview hash", {
  source_root <- withr::local_tempdir()
  archive_root <- tempfile("lna-update-mirror-")
  records <- update_test_archive(source_root, archive_root)
  preview <- plan_lna_dataset_update(
    archive_root, records[1L], source_root, mode = "mirror", prune = TRUE
  )
  expect_identical(update_statuses(preview)[["record-b"]], "source-missing")
  expect_error(update_lna_dataset(preview, hooks = update_test_hooks()), "confirm_prune")
  update_lna_dataset(
    preview, confirm_prune = preview$update_plan_hash,
    hooks = update_test_hooks()
  )
  expect_identical(read_lna_manifest(archive_root)$data$record_id, "record-a")
})

test_that("profile changes are stale and require refresh with complete inventory", {
  source_root <- withr::local_tempdir()
  archive_root <- tempfile("lna-update-profile-")
  records <- update_test_archive(source_root, archive_root)
  profile <- list(
    profile_id = "lna-independent-lossless-v2",
    loss_mode = "lossless", transforms = character(), transform_params = list()
  )
  expect_error(
    plan_lna_dataset_update(
      archive_root, records, source_root, mode = "append", codec_profile = profile
    ),
    "append cannot change"
  )
  preview <- plan_lna_dataset_update(
    archive_root, records, source_root, mode = "refresh", codec_profile = profile
  )
  expect_true(all(update_statuses(preview) == "stale"))
  state <- new.env(parent = emptyenv())
  update_lna_dataset(preview, hooks = update_test_hooks(state))
  expect_setequal(state$writes, c("record-a", "record-b"))
})

test_that("failed finalization leaves the prior generation valid and resumes", {
  source_root <- withr::local_tempdir()
  archive_root <- tempfile("lna-update-resume-")
  records <- update_test_archive(source_root, archive_root)
  prior <- read_lna_dataset_declaration(archive_root)$plan_hash
  records[[1L]] <- update_test_record(source_root, "record-a", "01", "v2")
  preview <- plan_lna_dataset_update(
    archive_root, records, source_root, mode = "refresh"
  )
  state <- new.env(parent = emptyenv())
  state$fail_manifest <- TRUE
  expect_error(
    update_lna_dataset(preview, hooks = update_test_hooks(state)),
    "injected update finalization failure"
  )
  expect_identical(read_lna_dataset_declaration(archive_root)$plan_hash, prior)
  expect_no_error(read_lna_manifest(archive_root))

  state$writes <- character()
  result <- update_lna_dataset(preview, hooks = update_test_hooks(state))
  expect_identical(state$writes, character())
  expect_identical(
    read_lna_dataset_declaration(archive_root)$plan_hash,
    preview$target_plan_hash
  )
  expect_true(dir.exists(result$previous_generation))
})

test_that("a demonstrably live update owner cannot be overridden", {
  source_root <- withr::local_tempdir()
  archive_root <- tempfile("lna-update-live-lock-")
  records <- update_test_archive(source_root, archive_root)
  records[[1L]] <- update_test_record(source_root, "record-a", "01", "v2")
  preview <- plan_lna_dataset_update(
    archive_root, records, source_root, mode = "refresh"
  )
  lock <- neuroarchive:::.update_lock(
    archive_root, preview, takeover_stale = FALSE, stale_after_seconds = 0
  )
  withr::defer(neuroarchive:::.update_release_lock(lock))

  expect_error(
    update_lna_dataset(
      preview, hooks = update_test_hooks(), takeover_stale = TRUE,
      stale_after_seconds = 0
    ),
    "demonstrably live"
  )
})

test_that("hard process death is recovered at every publication boundary", {
  checkpoints <- c(
    "locked", "staged", "old-preserved", "new-installed",
    "before-completion-receipt"
  )
  for (checkpoint in checkpoints) {
    source_root <- tempfile(paste0("lna-update-kill-source-", checkpoint, "-"))
    dir.create(source_root)
    withr::defer(unlink(source_root, recursive = TRUE, force = TRUE))
    archive_root <- tempfile(paste0("lna-update-kill-archive-", checkpoint, "-"))
    records <- update_test_archive(source_root, archive_root)
    records[[1L]] <- update_test_record(
      source_root, "record-a", "01", paste0("v2-", checkpoint)
    )
    preview <- plan_lna_dataset_update(
      archive_root, records, source_root, mode = "refresh"
    )
    update_test_hard_kill(preview, update_test_hooks(), checkpoint)

    if (identical(checkpoint, "old-preserved")) {
      expect_false(dir.exists(archive_root), info = checkpoint)
      state <- neuroarchive:::.update_read_state(archive_root, required = TRUE)
      backup <- file.path(dirname(archive_root), state$latest$backup_basename)
      expect_true(
        neuroarchive:::.update_generation_status(
          backup, preview$base_plan_hash
        )$valid,
        info = checkpoint
      )
    }
    expect_error(
      update_lna_dataset(
        preview, hooks = update_test_hooks(), stale_after_seconds = 0
      ),
      "takeover_stale"
    )
    recovered <- update_lna_dataset(
      preview, hooks = update_test_hooks(), takeover_stale = TRUE,
      stale_after_seconds = 0
    )
    expect_identical(
      read_lna_dataset_declaration(archive_root)$plan_hash,
      preview$target_plan_hash,
      info = checkpoint
    )
    expect_identical(
      read_lna_dataset_update(archive_root)$update_plan_hash,
      preview$update_plan_hash,
      info = checkpoint
    )
    expect_true(dir.exists(recovered$previous_generation), info = checkpoint)
    expect_false(
      dir.exists(neuroarchive:::.update_parent_paths(archive_root)$state),
      info = checkpoint
    )
    audit <- readLines(
      neuroarchive:::.update_parent_paths(archive_root)$audit,
      warn = FALSE
    )
    expect_true(any(grepl("stale-lock-takeover", audit, fixed = TRUE)),
      info = checkpoint)
  }
})

test_that("malformed transaction state fails without touching the live root", {
  source_root <- withr::local_tempdir()
  archive_root <- tempfile("lna-update-malformed-state-")
  records <- update_test_archive(source_root, archive_root)
  records[[1L]] <- update_test_record(source_root, "record-a", "01", "v2")
  preview <- plan_lna_dataset_update(
    archive_root, records, source_root, mode = "refresh"
  )
  prior <- read_lna_dataset_declaration(archive_root)$plan_hash
  state <- neuroarchive:::.update_parent_paths(archive_root)$state
  dir.create(state)
  jsonlite::write_json(
    list(schema_version = "1.0.0"),
    file.path(state, "000001-locked.json"), auto_unbox = TRUE
  )

  expect_error(
    update_lna_dataset(preview, hooks = update_test_hooks()),
    "transaction state is malformed"
  )
  expect_identical(
    read_lna_dataset_declaration(archive_root)$plan_hash, prior
  )
})
