query_test_checksum <- function(id) {
  paste0(
    "sha256:",
    digest::digest(charToRaw(id), algo = "sha256", serialize = FALSE)
  )
}

query_test_record <- function(record_id, entities, scope = "raw", pipeline = NULL) {
  lna_source_record(
    record_id = record_id,
    source_path = paste0("source/", record_id, ".nii.gz"),
    source_checksum = query_test_checksum(record_id),
    entities = entities,
    datatype = "func", suffix = "bold", extension = ".nii.gz",
    scope = scope, pipeline = pipeline,
    resolved_metadata = list(RepetitionTime = 0.8),
    geometry = list(
      dims = c(2L, 2L, 2L), affine = diag(c(2, 2, 2, 1)),
      units = "mm", voxel_ordering = "first-axis-fastest", n_volumes = 5L
    )
  )
}

query_test_archive <- function(records) {
  root <- tempfile("query-archive-")
  dir.create(root, recursive = TRUE)
  plan <- plan_lna_dataset(records, root)
  manifest <- lna_manifest(plan)
  write_lna_manifest(manifest, root)
  write_lna_dataset_declaration(lna_dataset_declaration(plan), root)
  for (relative in manifest$data$lna_path) {
    path <- file.path(root, relative)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines("not an HDF5 file", path, useBytes = TRUE)
  }
  list(root = root, plan = plan, manifest = manifest)
}

query_test_records <- function() {
  list(
    query_test_record(
      "rest-native", list(sub = "01", task = "rest", run = "01", space = "MNI")
    ),
    query_test_record(
      "rest-long", list(
        sub = "01", task = "rest", run = "02", echo = "01", acq = "mb8",
        space = "MNI152NLin2009cAsym", desc = "preproc"
      ), scope = "derivative", pipeline = "fmriprep-24.1"
    ),
    query_test_record(
      "rest-long-echo2", list(
        sub = "01", task = "rest", run = "02", echo = "02", acq = "mb8",
        space = "MNI152NLin2009cAsym", desc = "preproc"
      ), scope = "derivative", pipeline = "fmriprep-24.1"
    ),
    query_test_record(
      "resting-native", list(
        sub = "02", task = "resting", run = "01", space = "MNI"
      )
    )
  )
}

test_that("opening validates declaration and manifest without opening HDF5", {
  archive <- query_test_archive(query_test_records())
  dataset <- open_lna_dataset(archive$root)

  expect_s3_class(dataset, "lna_dataset_handle")
  expect_identical(dataset$declaration$profile, "NeuroarchiveBIDSExchange")
  expect_identical(nrow(dataset$manifest$data), 4L)
  expect_output(print(dataset), "<lna_dataset_handle>")
  expect_true(all(file.exists(file.path(
    archive$root, dataset$manifest$data$lna_path
  ))))
  expect_error(
    hdf5r::H5File$new(file.path(
      archive$root, dataset$manifest$data$lna_path[[1L]]
    ), mode = "r")
  )
  expect_no_error(query_lna(dataset, task = "rest"))
})

test_that("exact entity matching does not over-select prefixes", {
  dataset <- open_lna_dataset(query_test_archive(query_test_records())$root)

  rest <- query_lna(dataset, task = "rest")
  resting <- query_lna(dataset, task = "resting")
  mni <- query_lna(dataset, space = "MNI")
  mni_long <- query_lna(dataset, space = "MNI152NLin2009cAsym")

  expect_setequal(rest$record_id, c("rest-native", "rest-long", "rest-long-echo2"))
  expect_identical(resting$record_id, "resting-native")
  expect_setequal(mni$record_id, c("rest-native", "resting-native"))
  expect_setequal(mni_long$record_id, c("rest-long", "rest-long-echo2"))
  expect_true(all(c("record_id", "lna_path") %in% names(rest)))
})

test_that("multi-run, echo, acquisition, space, description, and core queries work", {
  dataset <- open_lna_dataset(query_test_archive(query_test_records())$root)
  selected <- query_lna(
    dataset,
    run = "02", echo = "02", acq = "mb8",
    space = "MNI152NLin2009cAsym", desc = "preproc",
    scope = "derivative", pipeline = "fmriprep-24.1",
    codec = "lna-independent-lossless-v1",
    status = "planned", view = "entities"
  )

  expect_identical(selected$record_id, "rest-long-echo2")
  expect_identical(selected$entity_echo, "02")
  expect_identical(selected$scope, "derivative")

  spatial <- query_lna(
    dataset, grid = dataset$manifest$data$grid_id[[1L]], view = "spatial"
  )
  expect_true(all(c("record_id", "lna_path", "grid_id", "reference_identifier") %in%
    names(spatial)))
  expect_true(nrow(spatial) >= 1L)
})

test_that("regex and glob matching are explicit and deterministic", {
  dataset <- open_lna_dataset(query_test_archive(query_test_records())$root)

  expect_length(query_lna(dataset, task = "rest*")$record_id, 0L)
  glob <- query_lna(dataset, task = "rest*", match = "glob")
  regex <- query_lna(dataset, space = "^MNI", match = "regex")

  expect_length(glob$record_id, 4L)
  expect_length(regex$record_id, 4L)
  expect_identical(glob$lna_path, sort(glob$lna_path))
  expect_error(
    query_lna(dataset, task = "(", match = "regex"),
    "invalid regex pattern",
    class = "lna_error_dataset_query"
  )
})

test_that("missing optional entities return empty views rather than over-select", {
  archive <- query_test_archive(list(query_test_record(
    "minimal", list(sub = "01", task = "rest", run = "01")
  )))
  dataset <- open_lna_dataset(archive$root)

  expect_identical(nrow(query_lna(dataset, echo = "01")), 0L)
  expect_identical(nrow(query_lna(dataset, custom = "anything")), 0L)
  entities <- lna_dataset_view(dataset, "entities")
  codec <- lna_dataset_view(dataset, "codec")
  metadata <- lna_dataset_view(dataset, "metadata")

  expect_true(all(c("record_id", "lna_path", "entity_task") %in% names(entities)))
  expect_true(all(c("record_id", "lna_path", "codec_profile_id") %in% names(codec)))
  expect_true(all(c("record_id", "lna_path", "metadata_provenance_json") %in%
    names(metadata)))
})

test_that("corrupt and incompatible authority files are refused", {
  archive <- query_test_archive(query_test_records())
  writeLines("{not-json", file.path(archive$root, "lna_manifest.json"))
  expect_error(
    open_lna_dataset(archive$root), "dataset manifest is invalid",
    class = "lna_error_dataset_query"
  )

  archive <- query_test_archive(query_test_records())
  declaration_path <- file.path(archive$root, "lna_dataset.json")
  declaration <- jsonlite::read_json(declaration_path, simplifyVector = FALSE)
  declaration$profile_version <- "2.0.0"
  jsonlite::write_json(
    declaration, declaration_path, auto_unbox = TRUE, null = "null", pretty = TRUE
  )
  expect_error(
    open_lna_dataset(archive$root), "dataset declaration is invalid",
    class = "lna_error_dataset_query"
  )
})
