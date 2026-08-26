writer_test_record <- function(source_root, record_id, sub) {
  relative <- paste0("sub-", sub, "/func/", record_id, ".source")
  path <- file.path(source_root, relative)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(paste0("SOURCE:", record_id), path, useBytes = TRUE)
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
      role = "metadata-source",
      path = metadata_relative,
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

writer_test_hooks <- function(state = new.env(parent = emptyenv())) {
  state$writes <- state$writes %||% character()
  state$validations <- state$validations %||% character()
  state$manifest_seen <- state$manifest_seen %||% FALSE
  list(
    check_space = function(plan, root) {
      if (identical(state$fail, "disk")) stop("injected disk-space failure")
      invisible(TRUE)
    },
    load_source = function(entry, source_path) {
      if (identical(state$fail, "load")) stop("injected transform input failure")
      paste(readLines(source_path, warn = FALSE), collapse = "")
    },
    write_artifact = function(data, temp_path, entry, plan) {
      state$writes <- c(state$writes, entry$record_id)
      if (identical(state$fail, "hdf5")) {
        writeLines("PARTIAL", temp_path)
        stop("injected HDF5 failure")
      }
      if (identical(state$fail, "transform")) stop("injected transform failure")
      writeLines(paste0("LNA:", entry$record_id, ":", data), temp_path, useBytes = TRUE)
    },
    validate_artifact = function(path, entry, stage) {
      state$validations <- c(state$validations, paste(entry$record_id, stage, sep = ":"))
      if (identical(state$fail, paste0("validate-", stage))) {
        stop(paste("injected validation failure", stage))
      }
      if (!startsWith(readLines(path, warn = FALSE)[[1L]], "LNA:")) {
        stop("not a test LNA")
      }
      TRUE
    },
    install_artifact = function(temp_path, destination, entry) {
      if (identical(state$fail, "rename")) return(FALSE)
      file.rename(temp_path, destination)
    },
    write_manifest = function(manifest, staging_root) {
      state$manifest_seen <- TRUE
      state$manifest_artifacts_exist <- all(file.exists(file.path(
        dirname(staging_root), manifest$data$lna_path
      )))
      state$manifest_absent_before_stage <- !file.exists(file.path(
        dirname(staging_root), "lna_manifest.tsv"
      ))
      if (identical(state$fail, "manifest")) stop("injected manifest failure")
      write_lna_manifest(manifest, staging_root)
    }
  )
}

writer_test_records <- function(source_root) {
  list(
    writer_test_record(source_root, "record-a", "01"),
    writer_test_record(source_root, "record-b", "02")
  )
}

test_that("multi-record writing is transactional and manifest-last", {
  source_root <- withr::local_tempdir()
  output_root <- withr::local_tempdir()
  state <- new.env(parent = emptyenv())
  hooks <- writer_test_hooks(state)
  result <- write_lna_dataset(
    writer_test_records(source_root), output_root, source_root,
    hooks = hooks, owner_id = "writer-success"
  )

  expect_s3_class(result, "lna_dataset_write_result")
  expect_setequal(result$written_record_ids, c("record-a", "record-b"))
  expect_identical(result$reused_record_ids, character())
  expect_true(all(file.exists(file.path(output_root, result$manifest$data$lna_path))))
  expect_true(all(result$manifest$data$status == "complete"))
  expect_true(state$manifest_seen)
  expect_true(state$manifest_artifacts_exist)
  expect_true(state$manifest_absent_before_stage)
  expect_true(file.exists(file.path(output_root, "lna_manifest.tsv")))
  expect_true(file.exists(file.path(output_root, "lna_manifest.json")))
  expect_true(file.exists(file.path(output_root, "lna_dataset.json")))
  expect_true(file.exists(file.path(output_root, ".lna_plan.json")))
  expect_false(dir.exists(file.path(output_root, ".lna.lock")))
  partials <- list.files(output_root, pattern = "lna-partial|backup-", recursive = TRUE)
  expect_length(partials, 0L)
  expect_true(all(c(
    "record-a:temporary", "record-a:installed", "record-a:final",
    "record-b:temporary", "record-b:installed", "record-b:final"
  ) %in% state$validations))
})

test_that("default artifact path writes and validates a real LNA", {
  source_root <- withr::local_tempdir()
  output_root <- withr::local_tempdir()
  record <- writer_test_record(source_root, "record-a", "01")
  result <- write_lna_dataset(
    record, output_root, source_root,
    hooks = list(load_source = function(entry, source_path) {
      array(seq_len(24), dim = c(2L, 2L, 2L, 3L))
    }),
    owner_id = "writer-real-lna"
  )
  path <- file.path(output_root, result$manifest$data$lna_path[[1L]])

  expect_true(hdf5r::is_hdf5(path))
  expect_true(validate_lna(path, strict = TRUE, checksum = TRUE))
  expect_identical(result$written_record_ids, "record-a")
  expect_identical(result$manifest$data$grid_id, record$geometry$grid_id)
})

test_that("declared path companions are copied byte-for-byte", {
  source_root <- withr::local_tempdir()
  output_root <- withr::local_tempdir()
  record <- writer_test_record(source_root, "record-a", "01")
  relative <- "sub-01/func/sub-01_task-rest_run-01_events.tsv"
  source <- file.path(source_root, relative)
  dir.create(dirname(source), recursive = TRUE, showWarnings = FALSE)
  writeLines(c("onset\tduration\ttrial_type", "0\t1\tgo"), source, useBytes = TRUE)
  record$associations <- c(record$associations, list(list(
    association_id = "events-a", role = "events", target_record_id = NULL,
    path = relative, uri = NULL,
    checksum = neuroarchive:::.journal_file_checksum(source), metadata = list()
  )))
  record <- as_lna_source_record(as.list(record))
  result <- write_lna_dataset(
    record, output_root, source_root, hooks = writer_test_hooks(),
    owner_id = "writer-companion"
  )
  copied <- file.path(output_root, relative)

  expect_true(file.exists(copied))
  expect_identical(readBin(copied, "raw", n = file.info(copied)$size),
    readBin(source, "raw", n = file.info(source)$size))
  expect_identical(
    result$manifest$data$associations_json[[1L]][[1L]]$role,
    "events"
  )
})

test_that("a completed dataset resumes from journal checksums", {
  source_root <- withr::local_tempdir()
  output_root <- withr::local_tempdir()
  records <- writer_test_records(source_root)
  first_state <- new.env(parent = emptyenv())
  first <- write_lna_dataset(
    records, output_root, source_root, hooks = writer_test_hooks(first_state),
    owner_id = "writer-first"
  )
  second_state <- new.env(parent = emptyenv())
  second <- write_lna_dataset(
    first$plan, source_root = source_root, hooks = writer_test_hooks(second_state),
    owner_id = "writer-second"
  )

  expect_identical(second$written_record_ids, character())
  expect_setequal(second$reused_record_ids, c("record-a", "record-b"))
  expect_identical(second_state$writes, character())
  expect_equal(second$manifest$data, first$manifest$data)
})

test_that("injected write-path failures leave no completed destination and resume", {
  failure_stages <- c(
    "transform", "hdf5", "validate-temporary", "rename", "validate-installed"
  )
  for (failure in failure_stages) {
    source_root <- withr::local_tempdir()
    output_root <- withr::local_tempdir()
    records <- list(writer_test_record(source_root, "record-a", "01"))
    state <- new.env(parent = emptyenv())
    state$fail <- failure
    expect_error(
      write_lna_dataset(
        records, output_root, source_root, hooks = writer_test_hooks(state),
        owner_id = paste0("writer-fail-", failure)
      ),
      "injected|atomic artifact installation failed"
    )
    plan <- plan_lna_dataset(records, output_root)
    destination <- file.path(output_root, plan$artifacts[[1L]]$destination_relative)
    expect_false(file.exists(destination), info = failure)
    expect_false(file.exists(file.path(output_root, "lna_manifest.tsv")), info = failure)
    expect_length(
      list.files(output_root, pattern = "lna-partial|backup-", recursive = TRUE),
      0L
    )
    state$fail <- NULL
    resumed <- write_lna_dataset(
      records, output_root, source_root, hooks = writer_test_hooks(state),
      owner_id = paste0("writer-resume-", failure)
    )
    expect_identical(resumed$written_record_ids, "record-a", info = failure)
    events <- read_lna_journal(output_root)$events
    attempts <- vapply(events[vapply(events, function(event) {
      identical(event$node_id, "write:record-a") && event$state == "started"
    }, logical(1))], `[[`, integer(1), "attempt")
    expect_identical(attempts, c(1L, 2L), info = failure)
  }
})

test_that("disk and manifest failures resume without rewriting completed artifacts", {
  for (failure in c("disk", "manifest")) {
    source_root <- withr::local_tempdir()
    output_root <- withr::local_tempdir()
    records <- writer_test_records(source_root)
    state <- new.env(parent = emptyenv())
    state$fail <- failure
    expect_error(write_lna_dataset(
      records, output_root, source_root, hooks = writer_test_hooks(state),
      owner_id = paste0("writer-", failure)
    ), "injected")
    writes_before <- length(state$writes %||% character())
    if (failure == "disk") expect_identical(writes_before, 0L)
    if (failure == "manifest") expect_identical(writes_before, 2L)
    state$fail <- NULL
    resumed <- write_lna_dataset(
      records, output_root, source_root, hooks = writer_test_hooks(state),
      owner_id = paste0("writer-resume-", failure)
    )
    if (failure == "manifest") {
      expect_setequal(resumed$reused_record_ids, c("record-a", "record-b"))
      expect_identical(length(state$writes), writes_before)
    }
    expect_true(file.exists(file.path(output_root, "lna_manifest.tsv")))
  }
})

test_that("pre-existing destinations require explicit transactional rebuild", {
  source_root <- withr::local_tempdir()
  output_root <- withr::local_tempdir()
  record <- writer_test_record(source_root, "record-a", "01")
  plan <- plan_lna_dataset(record, output_root)
  destination <- file.path(output_root, plan$artifacts[[1L]]$destination_relative)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  writeLines("ORIGINAL", destination)
  failing_state <- new.env(parent = emptyenv())
  failing_state$fail <- "validate-installed"

  expect_error(
    write_lna_dataset(
      record, output_root, source_root, hooks = writer_test_hooks(),
      owner_id = "writer-protected"
    ),
    "destination already exists"
  )
  expect_identical(readLines(destination), "ORIGINAL")
  expect_error(
    write_lna_dataset(
      record, output_root, source_root, hooks = writer_test_hooks(failing_state),
      overwrite = "rebuild", owner_id = "writer-rebuild-fail"
    ),
    "injected validation failure installed"
  )
  expect_identical(readLines(destination), "ORIGINAL")

  rebuilt <- write_lna_dataset(
    record, output_root, source_root, hooks = writer_test_hooks(),
    overwrite = "rebuild", owner_id = "writer-rebuild"
  )
  expect_identical(rebuilt$written_record_ids, "record-a")
  expect_true(startsWith(readLines(destination)[[1L]], "LNA:"))
})

test_that("source checksum changes and unauthorized authority files fail preflight", {
  source_root <- withr::local_tempdir()
  output_root <- withr::local_tempdir()
  record <- writer_test_record(source_root, "record-a", "01")
  source_path <- file.path(source_root, record$source$path)
  writeLines("CHANGED", source_path)
  expect_error(
    write_lna_dataset(
      record, output_root, source_root, hooks = writer_test_hooks(),
      owner_id = "writer-source-change"
    ),
    "source checksum does not match"
  )

  source_root <- withr::local_tempdir()
  output_root <- withr::local_tempdir()
  record <- writer_test_record(source_root, "record-a", "01")
  writeLines("FOREIGN", file.path(output_root, "lna_manifest.tsv"))
  expect_error(
    write_lna_dataset(
      record, output_root, source_root, hooks = writer_test_hooks(),
      owner_id = "writer-foreign-manifest"
    ),
    "authority files already exist"
  )
  expect_identical(readLines(file.path(output_root, "lna_manifest.tsv")), "FOREIGN")
})
