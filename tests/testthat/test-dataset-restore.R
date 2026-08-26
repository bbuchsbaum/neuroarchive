restore_test_record <- function(
    source_root,
    record_id,
    sub,
    task,
    run = "01",
    scope = "raw",
    pipeline = NULL,
    with_events = TRUE,
    target_record_id = NULL) {
  source_relative <- paste0("inputs/", record_id, ".source")
  source <- file.path(source_root, source_relative)
  dir.create(dirname(source), recursive = TRUE, showWarnings = FALSE)
  writeLines(paste0("SOURCE:", record_id), source, useBytes = TRUE)
  metadata_relative <- paste0("metadata/", record_id, ".json")
  metadata <- file.path(source_root, metadata_relative)
  dir.create(dirname(metadata), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(
    list(RepetitionTime = 0.8, TaskName = task), metadata,
    auto_unbox = TRUE, pretty = TRUE
  )
  associations <- list(list(
    association_id = paste0("metadata-", record_id),
    role = "metadata-source",
    path = metadata_relative,
    checksum = neuroarchive:::.journal_file_checksum(metadata),
    metadata = list(
      inheritance_level = "run", precedence = 100L,
      fields = c("RepetitionTime", "TaskName")
    )
  ))
  if (isTRUE(with_events)) {
    events_relative <- paste0(
      "sub-", sub, "/func/sub-", sub, "_task-", task,
      "_run-", run, "_events.tsv"
    )
    events <- file.path(source_root, events_relative)
    dir.create(dirname(events), recursive = TRUE, showWarnings = FALSE)
    writeLines(
      c("onset\tduration\ttrial_type", "0\t1\tcondition-a"),
      events, useBytes = TRUE
    )
    associations[[length(associations) + 1L]] <- list(
      association_id = paste0("events-", record_id),
      role = "events",
      path = events_relative,
      checksum = neuroarchive:::.journal_file_checksum(events),
      metadata = list(media_type = "text/tab-separated-values")
    )
  }
  if (!is.null(target_record_id)) {
    associations[[length(associations) + 1L]] <- list(
      association_id = paste0("dependency-", record_id),
      role = "dependency",
      target_record_id = target_record_id,
      metadata = list(reason = "paired acquisition")
    )
  }
  lna_source_record(
    record_id = record_id,
    source_path = source_relative,
    source_checksum = neuroarchive:::.journal_file_checksum(source),
    entities = list(sub = sub, task = task, run = run),
    datatype = "func", suffix = "bold", extension = ".nii.gz",
    scope = scope, pipeline = pipeline,
    resolved_metadata = list(RepetitionTime = 0.8, TaskName = task),
    associations = associations,
    geometry = list(
      dims = c(2L, 2L, 2L), affine = diag(c(2, 2, 2, 1)),
      units = "mm", voxel_ordering = "first-axis-fastest", n_volumes = 3L
    )
  )
}

restore_pack_hooks <- function() {
  list(
    load_source = function(entry, source_path) {
      paste(readLines(source_path, warn = FALSE), collapse = "")
    },
    write_artifact = function(data, temp_path, entry, plan) {
      writeLines(paste0("LNA:", entry$record_id, ":", data), temp_path,
        useBytes = TRUE)
    },
    validate_artifact = function(path, entry, stage) {
      stopifnot(startsWith(readLines(path, warn = FALSE)[[1L]], "LNA:"))
      TRUE
    },
    install_artifact = function(temp_path, destination, entry) {
      file.rename(temp_path, destination)
    }
  )
}

restore_unpack_hooks <- function() {
  list(
    reconstruct = function(lna_path, nifti_path, row, dataset) {
      writeLines(paste0("NIFTI:", row$record_id[[1L]]), nifti_path,
        useBytes = TRUE)
    },
    validate_nifti = function(path, row) {
      stopifnot(startsWith(readLines(path, warn = FALSE)[[1L]], "NIFTI:"))
      TRUE
    }
  )
}

restore_test_archive <- function(records, source_root) {
  archive_root <- tempfile("restore-archive-")
  write_lna_dataset(
    records, archive_root, source_root, hooks = restore_pack_hooks(),
    owner_id = paste0("pack-", basename(archive_root))
  )
  archive_root
}

restore_raw_records <- function(source_root) {
  list(
    restore_test_record(source_root, "sub01-rest", "01", "rest", "01"),
    restore_test_record(source_root, "sub01-motor", "01", "motor", "01"),
    restore_test_record(source_root, "sub02-rest", "02", "rest", "01")
  )
}

test_that("a full raw archive restores canonical paths, metadata, and events", {
  source_root <- withr::local_tempdir()
  records <- restore_raw_records(source_root)
  archive <- restore_test_archive(records, source_root)
  output <- tempfile("restored-raw-")
  receipt <- restore_bids(
    archive, output, hooks = restore_unpack_hooks(), closure = "closed"
  )

  expect_identical(receipt$status, "closed")
  expect_identical(receipt$missing_dependencies, list())
  expect_true(file.exists(file.path(output, "dataset_description.json")))
  expect_true(file.exists(file.path(output, "README")))
  expect_true(file.exists(file.path(output, ".bidsignore")))
  expect_true(file.exists(file.path(output, "lna_restoration.json")))
  expect_length(list.files(output, pattern = "[.]nii[.]gz$", recursive = TRUE), 3L)

  nifti <- file.path(
    output, "sub-01/func/sub-01_task-rest_run-01_bold.nii.gz"
  )
  sidecar <- sub("[.]nii[.]gz$", ".json", nifti)
  events <- file.path(
    output, "sub-01/func/sub-01_task-rest_run-01_events.tsv"
  )
  source_events <- file.path(
    source_root, "sub-01/func/sub-01_task-rest_run-01_events.tsv"
  )
  metadata <- jsonlite::read_json(sidecar, simplifyVector = TRUE)
  description <- jsonlite::read_json(
    file.path(output, "dataset_description.json"), simplifyVector = TRUE
  )

  expect_true(file.exists(nifti))
  expect_identical(metadata$RepetitionTime, 0.8)
  expect_identical(metadata$TaskName, "rest")
  expect_identical(description$DatasetType, "raw")
  expect_identical(description$BIDSVersion, "1.11.1")
  expect_identical(
    readBin(events, "raw", n = file.info(events)$size),
    readBin(source_events, "raw", n = file.info(source_events)$size)
  )
  expect_identical(
    receipt$transformation_rules$original_file_byte_identity, "not-claimed"
  )
  expect_identical(read_lna_restoration(output)$status, "closed")
  schema <- system.file(
    "schemas", "lna-restoration.schema.json", package = "neuroarchive"
  )
  if (nzchar(schema)) {
    expect_true(jsonvalidate::json_validate(
      file.path(output, "lna_restoration.json"), schema, engine = "ajv"
    ))
  }
})

test_that("default reconstruction writes a readable NIfTI from a real LNA", {
  source_root <- withr::local_tempdir()
  record <- restore_test_record(
    source_root, "real-lna", "01", "rest", "01", with_events = FALSE
  )
  archive <- tempfile("real-restore-archive-")
  write_lna_dataset(
    record, archive, source_root,
    hooks = list(load_source = function(entry, source_path) {
      array(seq_len(24), dim = c(2L, 2L, 2L, 3L))
    }),
    owner_id = "real-restore-pack"
  )
  output <- tempfile("real-restored-")
  restore_bids(archive, output)
  nifti <- file.path(
    output, "sub-01/func/sub-01_task-rest_run-01_bold.nii.gz"
  )
  restored <- neuroim2::read_vec(nifti)

  expect_true(file.exists(nifti))
  expect_identical(dim(restored), c(2L, 2L, 2L, 3L))
  expect_equal(as.vector(as.array(restored)), seq_len(24))
  expect_identical(methods::slot(restored, "header")$xyzt_units, 10L)
  expect_equal(methods::slot(restored, "header")$pixdim[[5L]], 0.8)
})

test_that("closed subject/task subsets preserve exact selection and closure", {
  source_root <- withr::local_tempdir()
  archive <- restore_test_archive(restore_raw_records(source_root), source_root)

  task_output <- tempfile("restored-task-")
  task_receipt <- restore_bids(
    archive, task_output, selection = list(task = "rest"),
    hooks = restore_unpack_hooks()
  )
  expect_setequal(
    task_receipt$selection$restored_record_ids,
    c("sub01-rest", "sub02-rest")
  )
  expect_length(
    list.files(task_output, pattern = "[.]nii[.]gz$", recursive = TRUE), 2L
  )
  expect_false(any(grepl("task-motor", list.files(task_output, recursive = TRUE))))

  subject_output <- tempfile("restored-subject-")
  subject_receipt <- restore_bids(
    archive, subject_output, selection = list(sub = "01"),
    hooks = restore_unpack_hooks()
  )
  expect_setequal(
    subject_receipt$selection$restored_record_ids,
    c("sub01-rest", "sub01-motor")
  )
  expect_false(any(grepl("sub-02", list.files(subject_output, recursive = TRUE))))
})

test_that("closed restoration includes transitive record dependencies", {
  source_root <- withr::local_tempdir()
  records <- list(
    restore_test_record(
      source_root, "paired-a", "01", "rest", "01",
      target_record_id = "paired-b"
    ),
    restore_test_record(source_root, "paired-b", "01", "rest", "02")
  )
  archive <- restore_test_archive(records, source_root)
  output <- tempfile("restored-pair-")
  receipt <- restore_bids(
    archive, output, selection = list(run = "01"),
    hooks = restore_unpack_hooks(), closure = "closed"
  )

  expect_identical(receipt$selection$requested_record_ids, "paired-a")
  expect_setequal(receipt$selection$restored_record_ids, c("paired-a", "paired-b"))
  expect_true(file.exists(file.path(
    output, "sub-01/func/sub-01_task-rest_run-02_bold.nii.gz"
  )))
})

test_that("derivative restoration declares GeneratedBy, Sources, and DatasetLinks", {
  source_root <- withr::local_tempdir()
  record <- restore_test_record(
    source_root, "derivative", "01", "rest", "01",
    scope = "derivative", pipeline = "fmriprep-24.1"
  )
  archive <- restore_test_archive(list(record), source_root)
  output <- tempfile("restored-derivative-")
  restore_bids(archive, output, hooks = restore_unpack_hooks())

  derivative_root <- file.path(output, "derivatives", "fmriprep-24.1")
  description <- jsonlite::read_json(
    file.path(derivative_root, "dataset_description.json"),
    simplifyVector = FALSE
  )
  sidecar <- jsonlite::read_json(file.path(
    derivative_root,
    "sub-01/func/sub-01_task-rest_run-01_bold.json"
  ), simplifyVector = FALSE)

  expect_identical(description$DatasetType, "derivative")
  expect_identical(description$GeneratedBy[[1L]]$Name, "fmriprep-24.1")
  expect_true("SourceDatasets" %in% names(description))
  expect_true("DatasetLinks" %in% names(description))
  expect_match(sidecar$Sources, "^bids:neuroarchive:")
})

test_that("closed missing dependencies fail before writes and partial is marked", {
  source_root <- withr::local_tempdir()
  archive <- restore_test_archive(list(
    restore_test_record(source_root, "missing-events", "01", "rest", "01")
  ), source_root)
  events <- file.path(
    archive, "sub-01/func/sub-01_task-rest_run-01_events.tsv"
  )
  unlink(events)
  closed_output <- tempfile("closed-missing-")
  error <- expect_error(
    restore_bids(
      archive, closed_output, hooks = restore_unpack_hooks(), closure = "closed"
    ),
    class = "lna_error_restoration"
  )
  expect_match(conditionMessage(error), "missing requirement")
  expect_false(dir.exists(closed_output))
  expect_true(length(error$missing_dependencies) >= 1L)

  partial_output <- tempfile("partial-missing-")
  receipt <- restore_bids(
    archive, partial_output, hooks = restore_unpack_hooks(), closure = "partial"
  )
  readme <- paste(readLines(file.path(partial_output, "README")), collapse = "\n")
  restored <- read_lna_restoration(partial_output)

  expect_identical(receipt$status, "partial")
  expect_match(receipt$bids_validity, "not-claimed")
  expect_true(length(receipt$missing_dependencies) >= 1L)
  expect_match(readme, "explicitly partial extraction")
  expect_identical(restored$status, "partial")
})

test_that("closed fixtures pass the local structural BIDS gate", {
  source_root <- withr::local_tempdir()
  archive <- restore_test_archive(restore_raw_records(source_root), source_root)
  output <- tempfile("restored-structural-")
  restore_bids(archive, output, hooks = restore_unpack_hooks())

  description <- jsonlite::read_json(
    file.path(output, "dataset_description.json"), simplifyVector = TRUE
  )
  nifti <- list.files(output, pattern = "_bold[.]nii[.]gz$", recursive = TRUE)
  sidecars <- sub("[.]nii[.]gz$", ".json", nifti)

  expect_true(all(c("Name", "BIDSVersion") %in% names(description)))
  expect_true(file.exists(file.path(output, "README")))
  expect_true(length(nifti) > 0L)
  expect_true(all(file.exists(file.path(output, sidecars))))
  expect_true(all(grepl(
    "^sub-[A-Za-z0-9+.-]+/func/sub-[A-Za-z0-9+.-]+_task-[A-Za-z0-9+.-]+_run-[A-Za-z0-9+.-]+_bold[.]nii[.]gz$",
    nifti
  )))
})
