transfer_test_record <- function(source_root, bytes = 65537L) {
  relative <- "sub-01/func/sub-01_task-rest_run-01_bold.source"
  path <- file.path(source_root, relative)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wb")
  writeBin(as.raw(seq_len(bytes) %% 251L), con)
  close(con)
  lna_source_record(
    record_id = "record-a",
    source_path = relative,
    source_checksum = neuroarchive:::.journal_file_checksum(path),
    entities = list(sub = "01", task = "rest", run = "01"),
    datatype = "func", suffix = "bold", extension = ".nii.gz",
    resolved_metadata = list(), associations = list(),
    geometry = list(
      dims = c(2L, 2L, 2L), affine = diag(c(2, 2, 2, 1)),
      units = "mm", voxel_ordering = "first-axis-fastest", n_volumes = 3L
    )
  )
}

transfer_test_writer_hooks <- function() {
  list(
    load_source = function(entry, source_path) source_path,
    write_artifact = function(data, temp_path, entry, plan) {
      if (!isTRUE(file.copy(data, temp_path, overwrite = FALSE))) stop("copy failed")
    },
    validate_artifact = function(path, entry, stage) {
      if (!file.exists(path) || file.info(path)$size < 1) stop("invalid test artifact")
      TRUE
    }
  )
}

transfer_test_archive <- function(source_root, archive_root, bytes = 65537L) {
  record <- transfer_test_record(source_root, bytes)
  write_lna_dataset(
    record, archive_root, source_root,
    hooks = transfer_test_writer_hooks(), owner_id = "transfer-source"
  )
  invisible(record)
}

transfer_test_read_receipt <- function(root) {
  jsonlite::read_json(
    file.path(root, "lna_transfer.json"), simplifyVector = FALSE
  )
}

transfer_test_write_receipt <- function(root, receipt) {
  jsonlite::write_json(
    receipt, file.path(root, "lna_transfer.json"), auto_unbox = TRUE,
    null = "null", digits = NA, pretty = TRUE
  )
}

transfer_test_rehash_receipt <- function(receipt) {
  receipt$summary <- list(
    n_files = length(receipt$files),
    total_bytes = sum(vapply(receipt$files, `[[`, numeric(1), "bytes"))
  )
  receipt$inventory_hash <- neuroarchive:::.transfer_inventory_hash(
    receipt$schema_version, receipt$files, receipt$summary
  )
  receipt
}

test_that("transfer inventory covers every portable file with strict digests", {
  source_root <- withr::local_tempdir()
  archive_root <- tempfile("lna-transfer-source-")
  transfer_test_archive(source_root, archive_root)
  inventory <- lna_transfer_inventory(archive_root)
  paths <- vapply(inventory$files, `[[`, character(1), "path")
  expect_true(all(c(
    "lna_dataset.json", "lna_manifest.json", "lna_manifest.tsv", ".lna_plan.json"
  ) %in% paths))
  expect_true(any(endsWith(paths, ".lna.h5")))
  expect_true(all(vapply(inventory$files, function(item) {
    grepl("^sha256:[0-9a-f]{64}$", item$checksum)
  }, logical(1))))
  expect_error(transfer_lna_dataset(archive_root, tempfile(), verify = FALSE),
    "cannot be disabled")
})

test_that("upload and download paths are byte exact and bounded memory", {
  source_root <- withr::local_tempdir()
  archive_root <- tempfile("lna-transfer-source-")
  upload_root <- tempfile("lna-transfer-upload-")
  download_root <- tempfile("lna-transfer-download-")
  transfer_test_archive(source_root, archive_root, bytes = 200003L)
  uploaded <- transfer_lna_dataset(
    archive_root, upload_root, direction = "upload", chunk_bytes = 4096L
  )
  expect_true(uploaded$verification$passed)
  expect_lte(uploaded$peak_buffer_bytes, 4096L)
  downloaded <- transfer_lna_dataset(
    upload_root, download_root, direction = "download", chunk_bytes = 2048L
  )
  expect_true(downloaded$verification$passed)
  expect_lte(downloaded$peak_buffer_bytes, 2048L)
  source_inventory <- lna_transfer_inventory(archive_root)
  received_inventory <- lna_transfer_inventory(download_root)
  source_checksums <- stats::setNames(
    vapply(source_inventory$files, `[[`, character(1), "checksum"),
    vapply(source_inventory$files, `[[`, character(1), "path")
  )
  received_checksums <- stats::setNames(
    vapply(received_inventory$files, `[[`, character(1), "checksum"),
    vapply(received_inventory$files, `[[`, character(1), "path")
  )
  expect_identical(received_checksums[names(source_checksums)], source_checksums)
})

test_that("interrupted partial transfer resumes from a verified prefix", {
  source_root <- withr::local_tempdir()
  archive_root <- tempfile("lna-transfer-source-")
  destination_root <- tempfile("lna-transfer-resume-")
  transfer_test_archive(source_root, archive_root, bytes = 50001L)
  state <- new.env(parent = emptyenv())
  state$failed <- FALSE
  expect_error(
    transfer_lna_dataset(
      archive_root, destination_root, chunk_bytes = 1024L,
      hooks = list(after_chunk = function(source, partial, transferred, expected_bytes) {
        if (!state$failed && endsWith(source, ".lna.h5") && transferred >= 3072L) {
          state$failed <- TRUE
          stop("injected network interruption")
        }
      })
    ),
    "injected network interruption"
  )
  partials <- list.files(
    destination_root, pattern = "lna-partial$", recursive = TRUE, full.names = TRUE
  )
  expect_length(partials, 1L)
  partial_size <- file.info(partials[[1L]])$size
  expect_gt(partial_size, 0)

  resumed <- transfer_lna_dataset(
    archive_root, destination_root, chunk_bytes = 1024L
  )
  lna_path <- names(resumed$outcomes)[endsWith(names(resumed$outcomes), ".lna.h5")]
  expect_identical(resumed$outcomes[[lna_path]], "resumed")
  expect_true(resumed$verification$passed)
  expect_length(list.files(
    destination_root, pattern = "lna-partial$", recursive = TRUE
  ), 0L)
})

test_that("receiver detects corruption and retry repairs it", {
  source_root <- withr::local_tempdir()
  archive_root <- tempfile("lna-transfer-source-")
  destination_root <- tempfile("lna-transfer-repair-")
  transfer_test_archive(source_root, archive_root)
  transfer_lna_dataset(archive_root, destination_root, chunk_bytes = 1024L)
  manifest <- read_lna_manifest(destination_root)$data
  artifact <- file.path(destination_root, manifest$lna_path[[1L]])
  con <- file(artifact, open = "r+b")
  writeBin(as.raw(255L), con)
  close(con)
  failed <- verify_lna_transfer(destination_root, error = FALSE)
  expect_false(failed$passed)
  expect_error(verify_lna_transfer(destination_root), "strict checksum")
  expect_error(
    transfer_lna_dataset(archive_root, destination_root, repair = FALSE),
    "repair is FALSE"
  )
  repaired <- transfer_lna_dataset(archive_root, destination_root, repair = TRUE)
  expect_true(repaired$verification$passed)
  expect_identical(
    neuroarchive:::.journal_file_checksum(artifact),
    neuroarchive:::.journal_file_checksum(file.path(archive_root, manifest$lna_path[[1L]]))
  )
})

test_that("receipt hash, summary, and record structure fail closed", {
  source_root <- withr::local_tempdir()
  archive_root <- tempfile("lna-transfer-source-")
  destination_root <- tempfile("lna-transfer-receipt-")
  transfer_test_archive(source_root, archive_root)
  transfer_lna_dataset(archive_root, destination_root)
  original <- transfer_test_read_receipt(destination_root)

  forged <- original
  forged$inventory_hash <- paste0("sha256:", strrep("0", 64L))
  transfer_test_write_receipt(destination_root, forged)
  expect_error(
    verify_lna_transfer(destination_root, error = FALSE),
    "inventory hash verification"
  )

  wrong_summary <- original
  wrong_summary$summary$n_files <- wrong_summary$summary$n_files + 1L
  transfer_test_write_receipt(destination_root, wrong_summary)
  expect_error(
    verify_lna_transfer(destination_root, error = FALSE),
    "summary does not match"
  )

  malformed <- original
  malformed$files[[1L]]$unexpected <- TRUE
  transfer_test_write_receipt(destination_root, malformed)
  expect_error(
    verify_lna_transfer(destination_root, error = FALSE),
    "must contain exactly"
  )

  malformed_checksum <- original
  malformed_checksum$files[[1L]]$checksum <- "sha256:not-a-digest"
  transfer_test_write_receipt(destination_root, malformed_checksum)
  expect_error(
    verify_lna_transfer(destination_root, error = FALSE),
    "checksum"
  )

  duplicate <- original
  duplicate$files[[2L]]$path <- duplicate$files[[1L]]$path
  transfer_test_write_receipt(destination_root, duplicate)
  expect_error(
    verify_lna_transfer(destination_root, error = FALSE),
    "duplicate paths"
  )

  invalid_bytes <- original
  invalid_bytes$files[[1L]]$bytes <- -1L
  transfer_test_write_receipt(destination_root, invalid_bytes)
  expect_error(
    verify_lna_transfer(destination_root, error = FALSE),
    "non-negative integer"
  )
})

test_that("an omitted manifest artifact cannot be certified", {
  source_root <- withr::local_tempdir()
  archive_root <- tempfile("lna-transfer-source-")
  destination_root <- tempfile("lna-transfer-omitted-")
  transfer_test_archive(source_root, archive_root)
  transfer_lna_dataset(archive_root, destination_root)
  receipt <- transfer_test_read_receipt(destination_root)
  paths <- vapply(receipt$files, `[[`, character(1), "path")
  artifact_index <- which(endsWith(paths, ".lna.h5"))
  expect_length(artifact_index, 1L)
  artifact <- file.path(destination_root, paths[[artifact_index]])
  expect_identical(unlink(artifact), 0L)
  receipt$files <- receipt$files[-artifact_index]
  receipt <- transfer_test_rehash_receipt(receipt)
  transfer_test_write_receipt(destination_root, receipt)

  verification <- verify_lna_transfer(destination_root, error = FALSE)
  expect_false(verification$passed)
  expect_match(
    verification$failures[["<archive-closure>"]],
    "missing required file"
  )
  expect_error(verify_lna_transfer(destination_root), "strict checksum")
})

test_that("an omitted declared companion cannot be certified", {
  source_root <- withr::local_tempdir()
  archive_root <- tempfile("lna-transfer-source-")
  destination_root <- tempfile("lna-transfer-companion-")
  record <- transfer_test_record(source_root)
  relative <- "sub-01/func/sub-01_task-rest_run-01_events.tsv"
  source <- file.path(source_root, relative)
  writeLines(c("onset\tduration\ttrial_type", "0\t1\tgo"), source,
    useBytes = TRUE)
  record$associations <- list(list(
    association_id = "events-a", role = "events", target_record_id = NULL,
    path = relative, uri = NULL,
    checksum = neuroarchive:::.journal_file_checksum(source), metadata = list()
  ))
  record <- as_lna_source_record(as.list(record))
  write_lna_dataset(
    record, archive_root, source_root,
    hooks = transfer_test_writer_hooks(), owner_id = "transfer-companion"
  )
  transfer_lna_dataset(archive_root, destination_root)
  receipt <- transfer_test_read_receipt(destination_root)
  paths <- vapply(receipt$files, `[[`, character(1), "path")
  companion_index <- which(paths == relative)
  expect_length(companion_index, 1L)
  expect_identical(unlink(file.path(destination_root, relative)), 0L)
  receipt$files <- receipt$files[-companion_index]
  transfer_test_write_receipt(
    destination_root, transfer_test_rehash_receipt(receipt)
  )

  verification <- verify_lna_transfer(destination_root, error = FALSE)
  expect_false(verification$passed)
  expect_match(
    verification$failures[["<archive-closure>"]],
    "missing required file"
  )
})

test_that("trusted sender hash distinguishes authenticity from self-consistency", {
  source_root <- withr::local_tempdir()
  archive_root <- tempfile("lna-transfer-source-")
  destination_root <- tempfile("lna-transfer-auth-")
  transfer_test_archive(source_root, archive_root)
  dir.create(file.path(archive_root, "metadata"), showWarnings = FALSE)
  writeLines("extra portable note", file.path(archive_root, "metadata", "note.txt"))
  transferred <- transfer_lna_dataset(archive_root, destination_root)
  trusted_hash <- transferred$inventory$inventory_hash
  receipt <- transfer_test_read_receipt(destination_root)
  paths <- vapply(receipt$files, `[[`, character(1), "path")
  note_index <- which(paths == "metadata/note.txt")
  expect_length(note_index, 1L)
  expect_identical(
    unlink(file.path(destination_root, "metadata", "note.txt")), 0L
  )
  receipt$files <- receipt$files[-note_index]
  receipt <- transfer_test_rehash_receipt(receipt)
  transfer_test_write_receipt(destination_root, receipt)

  self_consistent <- verify_lna_transfer(destination_root, error = FALSE)
  expect_true(self_consistent$passed)
  expect_false(self_consistent$authenticated)
  expect_error(
    verify_lna_transfer(
      destination_root, expected_inventory_hash = trusted_hash
    ),
    "trusted sender inventory hash"
  )
})

test_that("portable closure includes associations, shared bases, and participants", {
  manifest <- list(data = data.frame(
    lna_path = "sub-01/func/sub-01_task-rest_bold.lna.h5",
    stringsAsFactors = FALSE
  ))
  manifest$data$associations_json <- I(list(list(list(
    destination_relative = "sub-01/func/sub-01_task-rest_events.tsv"
  ))))
  declaration <- list(
    fit_resources = list(list(
      resource_type = "shared-basis",
      parameters = list(
        basis_filename = paste0("sha256-", strrep("a", 64L), ".lna_basis.h5")
      )
    )),
    participants_receipt = list(destination = "participants.tsv")
  )
  required <- neuroarchive:::.transfer_required_paths(manifest, declaration)
  expect_setequal(required, c(
    "lna_dataset.json", "lna_manifest.json", "lna_manifest.tsv",
    ".lna_plan.json", ".lna_journal.jsonl",
    "sub-01/func/sub-01_task-rest_bold.lna.h5",
    "sub-01/func/sub-01_task-rest_events.tsv",
    paste0("bases/sha256-", strrep("a", 64L), ".lna_basis.h5"),
    "bases/registry.json", "participants.tsv"
  ))
})

test_that("a changed transfer invalidates the old receipt before payload replacement", {
  first_source <- withr::local_tempdir()
  second_source <- withr::local_tempdir()
  first_archive <- tempfile("lna-transfer-first-")
  second_archive <- tempfile("lna-transfer-second-")
  destination_root <- tempfile("lna-transfer-stale-receipt-")
  transfer_test_archive(first_source, first_archive, bytes = 50001L)
  transfer_test_archive(second_source, second_archive, bytes = 50003L)
  transfer_lna_dataset(first_archive, destination_root)
  expect_true(file.exists(file.path(destination_root, "lna_transfer.json")))

  expect_error(
    transfer_lna_dataset(
      second_archive,
      destination_root,
      hooks = list(after_install = function(path, destination) {
        if (endsWith(path, ".lna.h5")) stop("injected replacement interruption")
      })
    ),
    "injected replacement interruption"
  )
  expect_false(file.exists(file.path(destination_root, "lna_transfer.json")))
})

test_that("source and destination roots must be disjoint", {
  source_root <- withr::local_tempdir()
  archive_root <- tempfile("lna-transfer-source-")
  transfer_test_archive(source_root, archive_root)
  nested <- file.path(archive_root, "nested-receiver")
  expect_error(transfer_lna_dataset(archive_root, nested), "must be disjoint")
  expect_false(dir.exists(nested))
})

test_that("case and Unicode-normalization path collisions fail before transfer", {
  expect_error(
    validate_lna_transfer_paths(c("sub-01/A.txt", "sub-01/a.txt")),
    "collide"
  )
  expect_error(
    validate_lna_transfer_paths(c(
      paste0("sub-01/", "\u00e9", ".txt"),
      paste0("sub-01/", "e\u0301", ".txt")
    )),
    "collide"
  )
  expect_no_error(validate_lna_transfer_paths(c("sub-01/a.txt", "sub-02/a.txt")))
})
