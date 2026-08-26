adapter_write_nifti <- function(path, affine, n_volumes = 3L, seed = 1L) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  set.seed(seed)
  values <- array(rnorm(2L * 3L * 2L * n_volumes), dim = c(2L, 3L, 2L, n_volumes))
  image <- neuroim2::NeuroVec(
    values,
    neuroim2::NeuroSpace(dim(values), trans = affine)
  )
  neuroim2::write_vec(image, path)
  invisible(values)
}

adapter_fake_api <- function(root, relative_path, event_path) {
  absolute_path <- file.path(root, relative_path)
  image_entities <- list(
    subid = "01", session = "A", task = "rest", acq = "mb",
    run = "01", echo = "2", kind = "bold", suffix = "nii.gz",
    extension = ".nii.gz", datatype = "func", type = "func"
  )
  event_entities <- list(
    subid = "01", session = "A", task = "rest", run = "01",
    kind = "events", suffix = "tsv", extension = ".tsv",
    datatype = "func", type = "func"
  )
  entity_row <- function(values) {
    tibble::as_tibble(stats::setNames(lapply(values, identity), names(values)))
  }
  api <- list(
    bids_project = function(...) NULL,
    query_files = function(x, full_path = FALSE, return = "paths", ...) {
      if (identical(return, "paths")) {
        return(if (isTRUE(full_path)) absolute_path else relative_path)
      }
      tibble::tibble(
        path = if (isTRUE(full_path)) absolute_path else relative_path,
        file = basename(relative_path),
        scope = "raw", pipeline = NA_character_, extension = ".nii.gz",
        datatype = "func", kind = "bold"
      )
    },
    bids_entities = function(paths, include_path = TRUE, coerce = TRUE) {
      values <- if (grepl("events[.]tsv$", paths)) event_entities else image_entities
      out <- entity_row(values)
      if (isTRUE(include_path)) {
        out <- tibble::add_column(out, .path = paths, .before = 1L)
      }
      out
    },
    get_metadata = function(x, file, inherit = TRUE, scope = "auto",
                            provenance = FALSE, ...) {
      result <- list(
        metadata = list(RepetitionTime = 0.8, TaskName = "rest"),
        sources = list(
          list(
            path = "task-rest_bold.json", inheritance_level = "task",
            precedence = 1L, fields = c("RepetitionTime", "TaskName")
          ),
          list(
            path = paste0(
              "sub-01/ses-A/func/",
              "sub-01_ses-A_task-rest_acq-mb_run-01_echo-2_bold.json"
            ),
            inheritance_level = "run", precedence = 2L,
            fields = "RepetitionTime"
          )
        )
      )
      if (isTRUE(provenance)) result else result$metadata
    },
    event_files = function(x, ..., full_path = TRUE) {
      if (isTRUE(full_path)) file.path(root, event_path) else event_path
    },
    read_events = function(...) tibble::tibble(),
    derivative_pipelines = function(...) tibble::tibble(),
    n_volumes = function(x, ...) stats::setNames(3L, x),
    version = "0.5.1"
  )
  api
}

test_that("bidser adapter version and public API contract fail actionably", {
  old <- expect_error(
    neuroarchive:::.bidser_adapter_validate_contract(
      "0.4.0", neuroarchive:::.lna_bidser_required_exports
    ),
    class = "lna_error_bidser_adapter"
  )
  expect_identical(old$field, "bidser_version")
  expect_match(conditionMessage(old), "bidser >= 0.5.1")

  incomplete <- expect_error(
    neuroarchive:::.bidser_adapter_validate_contract("0.5.1", "query_files"),
    class = "lna_error_bidser_adapter"
  )
  expect_identical(incomplete$field, "bidser_api")
  expect_match(conditionMessage(incomplete), "bids_entities")
})

test_that("adapter maps public discovery output to the source-record contract", {
  root <- tempfile("neuroarchive_bidser_adapter_")
  dir.create(root, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  relative <- paste0(
    "sub-01/ses-A/func/",
    "sub-01_ses-A_task-rest_acq-mb_run-01_echo-2_bold.nii.gz"
  )
  event <- paste0(
    "sub-01/ses-A/func/",
    "sub-01_ses-A_task-rest_run-01_events.tsv"
  )
  run_json <- paste0(
    "sub-01/ses-A/func/",
    "sub-01_ses-A_task-rest_acq-mb_run-01_echo-2_bold.json"
  )
  affine <- matrix(c(
    2, 0.1, 0, 7,
    0, 2.5, 0.2, 8,
    0, 0, 3, 9,
    0, 0, 0, 1
  ), nrow = 4L, byrow = TRUE)
  adapter_write_nifti(file.path(root, relative), affine)
  writeLines('{"RepetitionTime": 1.0, "TaskName": "rest"}', file.path(root, "task-rest_bold.json"))
  writeLines('{"RepetitionTime": 0.8}', file.path(root, run_json))
  dir.create(dirname(file.path(root, event)), recursive = TRUE, showWarnings = FALSE)
  writeLines("onset\tduration\ttrial_type\n0\t1\tgo", file.path(root, event))
  writeLines('{"trial_type": {"Description": "opaque condition label"}}', sub("[.]tsv$", ".json", file.path(root, event)))

  api <- adapter_fake_api(root, relative, event)
  records <- neuroarchive:::.lna_source_records_from_bidser(
    project = structure(list(), class = "fake_bidser_project"), api = api
  )

  expect_length(records, 1L)
  record <- records[[1L]]
  expect_s3_class(record, "lna_source_record")
  expect_identical(
    record$entities,
    list(acq = "mb", echo = "2", run = "01", ses = "A", sub = "01", task = "rest")
  )
  expect_identical(record$bids$suffix, "bold")
  expect_identical(record$bids$extension, ".nii.gz")
  expect_identical(record$bids$scope, "raw")
  expect_identical(record$resolved_metadata$RepetitionTime, 0.8)
  expect_identical(record$geometry$dims, c(2L, 3L, 2L))
  expect_identical(record$geometry$n_volumes, 3L)
  expect_equal(record$geometry$affine, affine, tolerance = 1e-7)
  expect_identical(record$reference_assertion$kind, "unknown-native")
  roles <- vapply(record$associations, `[[`, character(1), "role")
  expect_setequal(
    roles,
    c("metadata-source", "metadata-source", "events", "events-json")
  )

  plan <- plan_lna_dataset(
    records,
    file.path(root, "archive"),
    metadata_policy = lna_metadata_policy(portable = TRUE)
  )
  expect_identical(
    plan$artifacts[[1L]]$metadata_provenance$status, "complete"
  )
  expect_identical(
    plan$artifacts[[1L]]$metadata_provenance$field_sources$RepetitionTime$precedence,
    2L
  )

  no_events <- neuroarchive:::.lna_source_records_from_bidser(
    project = structure(list(), class = "fake_bidser_project"),
    api = api,
    include_events = FALSE
  )[[1L]]
  expect_false(any(vapply(no_events$associations, function(x) {
    startsWith(x$role, "events")
  }, logical(1))))
})

test_that("real bidser adapter covers raw, derivative, multi-echo, and events", {
  if (!requireNamespace("bidser", quietly = TRUE) ||
      utils::packageVersion("bidser") < base::package_version("0.5.1")) {
    skip("bidser >= 0.5.1 integration baseline is not installed")
  }

  root <- tempfile("neuroarchive_real_bidser_")
  dir.create(root, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines("Adapter integration fixture", file.path(root, "README"))
  writeLines(
    '{"Name":"adapter fixture","BIDSVersion":"1.11.1"}',
    file.path(root, "dataset_description.json")
  )
  writeLines(
    "participant_id\nsub-01\nsub-02",
    file.path(root, "participants.tsv")
  )
  writeLines(
    '{"RepetitionTime":0.8,"TaskName":"rest"}',
    file.path(root, "task-rest_bold.json")
  )
  affine <- diag(c(2, 2, 2.5, 1))
  raw_paths <- c(
    paste0(
      "sub-01/ses-A/func/",
      "sub-01_ses-A_task-rest_acq-mb_run-01_echo-1_bold.nii.gz"
    ),
    paste0(
      "sub-01/ses-A/func/",
      "sub-01_ses-A_task-rest_acq-mb_run-01_echo-2_bold.nii.gz"
    ),
    "sub-02/func/sub-02_task-nback_run-02_bold.nii.gz"
  )
  invisible(lapply(seq_along(raw_paths), function(i) {
    adapter_write_nifti(file.path(root, raw_paths[[i]]), affine, seed = i)
  }))
  writeLines(
    "onset\tduration\ttrial_type\n0\t1\tgo",
    file.path(
      root, "sub-01/ses-A/func/",
      "sub-01_ses-A_task-rest_run-01_events.tsv"
    )
  )

  derivative_root <- file.path(root, "derivatives", "fmriprep")
  dir.create(derivative_root, recursive = TRUE)
  writeLines(
    paste0(
      '{"Name":"fmriprep","BIDSVersion":"1.11.1",',
      '"DatasetType":"derivative","GeneratedBy":[{"Name":"fmriprep"}]}'
    ),
    file.path(derivative_root, "dataset_description.json")
  )
  writeLines(
    '{"RepetitionTime":0.8,"SpatialReference":"MNI152NLin2009cAsym"}',
    file.path(derivative_root, "task-rest_bold.json")
  )
  derivative_path <- paste0(
    "derivatives/fmriprep/sub-01/ses-A/func/",
    "sub-01_ses-A_task-rest_acq-mb_run-01_echo-1_",
    "space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz"
  )
  adapter_write_nifti(file.path(root, derivative_path), affine, seed = 10L)

  project <- bidser::bids_project(root, derivatives = "auto")
  records <- lna_source_records_from_bidser(project)
  expect_length(records, 4L)
  expect_setequal(
    vapply(records, function(record) record$entities$sub, character(1)),
    c("01", "02")
  )
  derivative <- records[vapply(records, function(record) {
    identical(record$bids$scope, "derivative")
  }, logical(1))][[1L]]
  expect_identical(derivative$bids$pipeline, "fmriprep")
  expect_identical(
    derivative$reference_assertion$identifier,
    "MNI152NLin2009cAsym"
  )
  expect_identical(
    derivative$reference_assertion$origin, "bids-space-entity"
  )
  raw_rest <- records[vapply(records, function(record) {
    identical(record$bids$scope, "raw") && identical(record$entities$task, "rest")
  }, logical(1))]
  expect_setequal(
    vapply(raw_rest, function(record) record$entities$echo, character(1)),
    c("1", "2")
  )
  expect_true(all(vapply(raw_rest, function(record) {
    "events" %in% vapply(record$associations, `[[`, character(1), "role")
  }, logical(1))))
  expect_no_error(plan_lna_dataset(
    records,
    file.path(root, "portable-archive"),
    metadata_policy = lna_metadata_policy(portable = TRUE)
  ))
})

test_that("core dataset source records remain independent of bidser", {
  expect_false("bidser" %in% names(getNamespaceImports("neuroarchive")))
  record <- lna_source_record(
    record_id = "bidser-free",
    source_path = "sub-01/func/sub-01_task-rest_bold.nii.gz",
    source_checksum = paste0("sha256:", strrep("a", 64L)),
    entities = list(sub = "01", task = "rest"),
    datatype = "func", suffix = "bold", extension = ".nii.gz",
    geometry = list(
      dims = c(2L, 2L, 2L), affine = diag(4), units = "mm",
      voxel_ordering = "first-axis-fastest", n_volumes = 2L
    )
  )
  expect_s3_class(record, "lna_source_record")
})
