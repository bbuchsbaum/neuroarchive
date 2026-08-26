metadata_test_checksum <- function(path) {
  neuroarchive:::.journal_file_checksum(path)
}

metadata_test_association <- function(
    id,
    role,
    relative,
    source_root,
    metadata = list()) {
  list(
    association_id = id,
    role = role,
    path = relative,
    checksum = metadata_test_checksum(file.path(source_root, relative)),
    metadata = metadata
  )
}

metadata_test_record <- function(
    source_root,
    record_id = "record-01",
    sub = "01",
    extra_associations = list(),
    metadata_sources = NULL) {
  relative <- paste0("sub-", sub, "/func/", record_id, ".source")
  source <- file.path(source_root, relative)
  dir.create(dirname(source), recursive = TRUE, showWarnings = FALSE)
  writeLines(paste0("SOURCE:", record_id), source, useBytes = TRUE)
  if (is.null(metadata_sources)) {
    metadata_relative <- paste0("sub-", sub, "/func/", record_id, ".json")
    metadata_path <- file.path(source_root, metadata_relative)
    jsonlite::write_json(
      list(RepetitionTime = 0.8, TaskName = "rest"), metadata_path,
      auto_unbox = TRUE, pretty = TRUE
    )
    metadata_sources <- list(metadata_test_association(
      paste0("metadata-", record_id), "metadata-source", metadata_relative,
      source_root,
      list(
        inheritance_level = "run", precedence = 100L,
        fields = c("RepetitionTime", "TaskName"),
        media_type = "application/json"
      )
    ))
  }
  lna_source_record(
    record_id = record_id,
    source_path = relative,
    source_checksum = metadata_test_checksum(source),
    entities = list(sub = sub, task = "rest", run = "01"),
    datatype = "func", suffix = "bold", extension = ".nii.gz",
    resolved_metadata = list(RepetitionTime = 0.8, TaskName = "rest"),
    associations = c(metadata_sources, extra_associations),
    geometry = list(
      dims = c(2L, 2L, 2L), affine = diag(c(2, 2, 2, 1)),
      units = "mm", voxel_ordering = "first-axis-fastest", n_volumes = 3L
    )
  )
}

metadata_test_hooks <- function() {
  list(
    load_source = function(entry, source_path) {
      paste(readLines(source_path, warn = FALSE), collapse = "")
    },
    write_artifact = function(data, temp_path, entry, plan) {
      writeLines(paste0("LNA:", entry$record_id, ":", data), temp_path,
        useBytes = TRUE)
    },
    validate_artifact = function(path, entry, stage) {
      if (!startsWith(readLines(path, warn = FALSE)[[1L]], "LNA:")) {
        stop("not a test LNA")
      }
      TRUE
    },
    install_artifact = function(temp_path, destination, entry) {
      file.rename(temp_path, destination)
    }
  )
}

test_that("effective inherited metadata has deterministic field provenance", {
  source_root <- withr::local_tempdir()
  paths <- c(
    dataset = "task-rest_bold.dataset.json",
    task = "task-rest_bold.task.json",
    run = "sub-01_task-rest_run-01_bold.json"
  )
  contents <- list(
    dataset = list(Manufacturer = "ScannerCo"),
    task = list(TaskName = "rest", RepetitionTime = 1.0),
    run = list(RepetitionTime = 0.8)
  )
  levels <- c(dataset = "dataset", task = "task", run = "run")
  precedences <- c(dataset = 1L, task = 10L, run = 100L)
  fields <- list(
    dataset = "Manufacturer",
    task = c("TaskName", "RepetitionTime"),
    run = "RepetitionTime"
  )
  associations <- lapply(names(paths), function(name) {
    jsonlite::write_json(
      contents[[name]], file.path(source_root, paths[[name]]),
      auto_unbox = TRUE
    )
    metadata_test_association(
      paste0("metadata-", name), "metadata-source", paths[[name]], source_root,
      list(
        inheritance_level = levels[[name]],
        precedence = precedences[[name]],
        fields = fields[[name]]
      )
    )
  })
  source_relative <- "sub-01/func/source.nii.gz"
  dir.create(dirname(file.path(source_root, source_relative)), recursive = TRUE)
  writeLines("source", file.path(source_root, source_relative), useBytes = TRUE)
  record <- lna_source_record(
    record_id = "inherited",
    source_path = source_relative,
    source_checksum = metadata_test_checksum(file.path(source_root, source_relative)),
    entities = list(sub = "01", task = "rest", run = "01"),
    datatype = "func", suffix = "bold", extension = ".nii.gz",
    resolved_metadata = list(
      Manufacturer = "ScannerCo", TaskName = "rest", RepetitionTime = 0.8
    ),
    associations = associations,
    geometry = list(
      dims = c(2L, 2L, 2L), affine = diag(4), units = "mm",
      voxel_ordering = "first-axis-fastest", n_volumes = 3L
    )
  )
  plan <- plan_lna_dataset(
    record, "archive", metadata_policy = lna_metadata_policy(portable = TRUE)
  )
  provenance <- plan$artifacts[[1L]]$metadata_provenance

  expect_identical(provenance$status, "complete")
  expect_identical(
    provenance$field_sources$RepetitionTime$association_id, "metadata-run"
  )
  expect_identical(
    provenance$field_sources$TaskName$association_id, "metadata-task"
  )
  expect_identical(
    provenance$field_sources$Manufacturer$association_id, "metadata-dataset"
  )
  manifest_provenance <- lna_manifest(plan)$data$metadata_provenance_json[[1L]]
  expect_identical(
    manifest_provenance$field_sources$RepetitionTime$association_id,
    "metadata-run"
  )
})

test_that("portable metadata provenance and required companions fail closed", {
  source_root <- withr::local_tempdir()
  record <- metadata_test_record(source_root)
  record$associations <- list()
  record <- as_lna_source_record(as.list(record))

  error <- expect_error(
    plan_lna_dataset(
      record, "archive", metadata_policy = lna_metadata_policy(portable = TRUE)
    ),
    class = "lna_error_metadata_policy"
  )
  expect_match(conditionMessage(error), "lacks provenance")
  expect_identical(error$record_id, "record-01")

  expect_error(
    plan_lna_dataset(
      metadata_test_record(source_root), "archive",
      metadata_policy = lna_metadata_policy(portable = TRUE, events = "require")
    ),
    "required events companion is missing",
    class = "lna_error_metadata_policy"
  )
})

test_that("selected events and their JSON are byte-preserved and remain opaque", {
  source_root <- withr::local_tempdir()
  output_root <- withr::local_tempdir()
  events_relative <- "sub-01/func/sub-01_task-rest_run-01_events.tsv"
  events_json_relative <- "sub-01/func/sub-01_task-rest_run-01_events.json"
  unselected_relative <- "sub-01/func/sub-01_task-rest_run-02_events.tsv"
  dir.create(dirname(file.path(source_root, events_relative)), recursive = TRUE)
  writeLines(
    c("onset\tduration\ttrial_type", "0\t1\tgo", "2\t1\tstop"),
    file.path(source_root, events_relative), useBytes = TRUE
  )
  writeLines(
    '{"trial_type":{"Description":"opaque condition label"}}',
    file.path(source_root, events_json_relative), useBytes = TRUE
  )
  writeLines(
    c("onset\tduration\ttrial_type", "0\t1\tunselected"),
    file.path(source_root, unselected_relative), useBytes = TRUE
  )
  associations <- list(
    metadata_test_association(
      "events-tsv", "events", events_relative, source_root,
      list(media_type = "text/tab-separated-values")
    ),
    metadata_test_association(
      "events-json", "events-json", events_json_relative, source_root,
      list(media_type = "application/json")
    )
  )
  record <- metadata_test_record(
    source_root, extra_associations = associations
  )
  result <- write_lna_dataset(
    record, output_root, source_root, hooks = metadata_test_hooks(),
    owner_id = "metadata-events"
  )

  for (relative in c(events_relative, events_json_relative)) {
    source <- file.path(source_root, relative)
    destination <- file.path(output_root, relative)
    expect_true(file.exists(destination))
    expect_identical(
      readBin(destination, "raw", n = file.info(destination)$size),
      readBin(source, "raw", n = file.info(source)$size)
    )
  }
  expect_false(file.exists(file.path(output_root, unselected_relative)))
  declaration <- read_lna_dataset_declaration(output_root)
  expect_identical(
    declaration$metadata_policy$task_semantics$task_entity,
    "acquisition-label"
  )
  expect_identical(
    declaration$metadata_policy$task_semantics$within_run_conditions,
    "opaque-events-content"
  )
  expect_true(all(c("events", "events-json") %in%
    vapply(result$plan$artifacts[[1L]]$expected_companions,
      `[[`, character(1), "role")))
})

test_that("participants are excluded by default and explicitly privacy-subset", {
  source_root <- withr::local_tempdir()
  participants <- file.path(source_root, "participants.tsv")
  writeLines(
    c(
      "participant_id\tage\tprivate_note",
      "sub-01\t20\talpha",
      "sub-02\t30\tbeta",
      "sub-03\t40\tgamma"
    ),
    participants, useBytes = TRUE
  )
  records <- list(
    metadata_test_record(source_root, "record-01", "01"),
    metadata_test_record(source_root, "record-02", "02")
  )

  excluded_root <- withr::local_tempdir()
  write_lna_dataset(
    records, excluded_root, source_root, hooks = metadata_test_hooks(),
    owner_id = "participants-excluded"
  )
  expect_false(file.exists(file.path(excluded_root, "participants.tsv")))
  excluded <- read_lna_dataset_declaration(excluded_root)
  expect_identical(excluded$metadata_policy$participants$mode, "exclude")

  included_root <- withr::local_tempdir()
  policy <- lna_metadata_policy(
    portable = TRUE,
    participants = list(
      mode = "include-subset",
      source_path = "participants.tsv",
      source_checksum = metadata_test_checksum(participants),
      columns = c("participant_id", "age")
    )
  )
  write_lna_dataset(
    records, included_root, source_root, hooks = metadata_test_hooks(),
    metadata_policy = policy, owner_id = "participants-included"
  )
  subset <- utils::read.delim(
    file.path(included_root, "participants.tsv"), colClasses = "character",
    check.names = FALSE
  )
  declaration <- read_lna_dataset_declaration(included_root)

  expect_identical(subset$participant_id, c("sub-01", "sub-02"))
  expect_identical(names(subset), c("participant_id", "age"))
  expect_false("sub-03" %in% subset$participant_id)
  expect_identical(
    declaration$participants_receipt$selected_participant_ids,
    c("sub-01", "sub-02")
  )
  expect_match(declaration$participants_receipt$output_checksum,
    "^sha256:[0-9a-f]{64}$")
})

test_that("typed association modes and declaration round-trip are explicit", {
  source_root <- withr::local_tempdir()
  record <- metadata_test_record(source_root)
  policy <- lna_metadata_policy(
    portable = TRUE, sessions = "copy-if-present", scans = "exclude",
    confounds = "require", physio = "copy-if-present"
  )
  expect_identical(policy$companions$sessions, "copy-if-present")
  expect_identical(policy$companions$scans, "exclude")
  expect_identical(policy$companions$confounds, "require")
  expect_identical(policy$companions$physio, "copy-if-present")
  expect_error(
    plan_lna_dataset(record, "archive", metadata_policy = policy),
    "required confounds companion is missing"
  )

  plan <- plan_lna_dataset(
    record, "archive", metadata_policy = lna_metadata_policy(portable = TRUE)
  )
  declaration <- lna_dataset_declaration(plan)
  root <- withr::local_tempdir()
  write_lna_dataset_declaration(declaration, root)
  restored <- read_lna_dataset_declaration(root)
  schema <- system.file(
    "schemas", "lna-dataset.schema.json", package = "neuroarchive"
  )

  expect_identical(restored, declaration)
  expect_identical(restored$extension_registry[[1L]]$extension, ".lna.h5")
  expect_identical(
    restored$extension_registry[[1L]]$media_type,
    "application/vnd.neuroarchive.lna+hdf5"
  )
  expect_false(restored$extension_registry[[1L]]$bids_recognized_extension)
  legacy_declaration <- declaration
  legacy_declaration$manifest$schema_version <- "1.0.0"
  expect_identical(
    validate_lna_dataset_declaration(legacy_declaration)$manifest$schema_version,
    "1.0.0"
  )
  if (nzchar(schema)) {
    expect_true(jsonvalidate::json_validate(
      file.path(root, "lna_dataset.json"), schema, engine = "ajv"
    ))
  }
})
