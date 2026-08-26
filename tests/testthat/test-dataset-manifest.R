manifest_test_checksum <- function(id) {
  paste0(
    "sha256:",
    digest::digest(charToRaw(id), algo = "sha256", serialize = FALSE)
  )
}

manifest_test_record <- function(
    record_id,
    entities,
    scope = "raw",
    pipeline = NULL,
    datatype = "func",
    suffix = "bold",
    source_path = NULL,
    source_uri = NULL,
    metadata = list(RepetitionTime = 0.8),
    associations = list()) {
  if (is.null(source_path) && is.null(source_uri)) {
    source_path <- paste0("sub-", entities$sub, "/source/", record_id, ".nii.gz")
  }
  lna_source_record(
    record_id = record_id,
    source_path = source_path,
    source_uri = source_uri,
    source_checksum = manifest_test_checksum(record_id),
    entities = entities,
    datatype = datatype,
    suffix = suffix,
    extension = ".nii.gz",
    scope = scope,
    pipeline = pipeline,
    resolved_metadata = metadata,
    associations = associations,
    geometry = list(
      dims = c(5L, 4L, 3L),
      affine = diag(c(2, 2, 2, 1)),
      units = "mm",
      voxel_ordering = "first-axis-fastest",
      n_volumes = 12L
    )
  )
}

manifest_test_plan <- function(records) {
  plan_lna_dataset(records, "portable-archive")
}

test_that("manifest declares a fixed typed core and optional exact entities", {
  record <- manifest_test_record(
    "raw-rest", list(sub = "01", ses = "base", task = "rest", run = "01")
  )
  manifest <- lna_manifest(manifest_test_plan(record))

  expect_s3_class(manifest, "lna_manifest")
  expect_identical(manifest$schema$manifest_schema_version, "1.1.0")
  expect_identical(manifest$schema$null_encoding, "n/a")
  expect_identical(manifest$schema$row_semantics, "one-row-per-source-image")
  expect_identical(
    manifest$schema$column_order,
    names(manifest$schema$columns)
  )
  expect_true(all(c(
    "source_checksum", "source_extension", "source_bytes", "lna_checksum", "lna_bytes",
    "codec_profile_json", "shared_basis_json", "associations_json", "metadata_json",
    "metadata_provenance_json",
    "n_volumes", "repetition_time_seconds", "reference_identifier",
    "grid_id", "support_id", "status"
  ) %in% names(manifest$data)))
  expect_true(all(c(
    "entity_sub", "entity_ses", "entity_task", "entity_run"
  ) %in% names(manifest$data)))
  expect_identical(manifest$schema$entity_semantics$comparison, "exact")
  expect_identical(manifest$schema$columns$source_bytes$unit, "byte")
  expect_identical(
    manifest$schema$columns$repetition_time_seconds$unit, "second"
  )
  expect_identical(manifest$data$status, "planned")
  expect_output(print(manifest), "<lna_manifest>")
})

test_that("raw, derivative, multi-echo, and sparse entities round-trip", {
  records <- list(
    manifest_test_record(
      "raw-rest", list(sub = "01", task = "rest", run = "01")
    ),
    manifest_test_record(
      "derivative-rest",
      list(
        sub = "01", task = "rest", run = "01",
        space = "MNI152NLin2009cAsym", res = "2", desc = "preproc"
      ),
      scope = "derivative", pipeline = "fmriprep-24.1"
    ),
    manifest_test_record(
      "echo-01",
      list(sub = "02", ses = "01", task = "motor", run = "02", echo = "01")
    ),
    manifest_test_record(
      "echo-02",
      list(sub = "02", ses = "01", task = "motor", run = "02", echo = "02")
    ),
    manifest_test_record("sparse", list(sub = "03"), metadata = list())
  )
  manifest <- lna_manifest(manifest_test_plan(records))
  root <- withr::local_tempdir()

  paths <- write_lna_manifest(manifest, root)
  restored <- read_lna_manifest(paths$tsv)

  expect_s3_class(restored, "lna_manifest")
  expect_identical(restored$schema, manifest$schema)
  expect_equal(restored$data, manifest$data)
  expect_true(any(is.na(restored$data$entity_echo)))
  expect_setequal(restored$data$entity_echo[!is.na(restored$data$entity_echo)], c("01", "02"))
  expect_identical(
    restored$data$pipeline[restored$data$scope == "derivative"],
    "fmriprep-24.1"
  )
  expect_identical(
    restored$data$reference_origin,
    rep("unknown", nrow(restored$data))
  )
})

test_that("artifact sizes, checksums, codec, temporal data, and status survive", {
  record <- manifest_test_record(
    "complete", list(sub = "01", task = "rest", run = "01"),
    metadata = list(RepetitionTime = 1.25, TaskName = "rest")
  )
  plan <- plan_lna_dataset(
    record,
    "portable-archive",
    codec_profile = list(
      profile_id = "qpp-analysis-safe-v1",
      loss_mode = "lossy",
      transforms = c("parcel.mean", "qpp"),
      transform_params = list(qpp = list(bits = 12L)),
      error_contract = lna_error_contract(
        "lossy", "float32", max_abs_error = 0.1, rmse_max = 0.05
      )
    )
  )
  artifact_state <- list(complete = list(
    source_bytes = 12345,
    lna_checksum = manifest_test_checksum("lna-complete"),
    lna_bytes = 4321,
    status = "complete"
  ))
  manifest <- lna_manifest(plan, artifact_state)
  row <- manifest$data[1L, , drop = FALSE]

  expect_identical(row$source_bytes, 12345)
  expect_identical(row$lna_bytes, 4321)
  expect_identical(row$status, "complete")
  expect_identical(row$loss_mode, "lossy")
  expect_identical(row$n_volumes, 12)
  expect_identical(row$repetition_time_seconds, 1.25)
  expect_identical(row$metadata_json[[1L]]$TaskName, "rest")
  expect_identical(row$codec_profile_json[[1L]]$transform_params$qpp$bits, 12L)
})

test_that("type, nullability, unit, and complete-state violations fail", {
  record <- manifest_test_record(
    "invalid", list(sub = "01", task = "rest", run = "01")
  )
  manifest <- lna_manifest(manifest_test_plan(record))

  wrong_type <- manifest
  wrong_type$schema$columns$source_bytes$type <- "string"
  expect_error(
    validate_lna_manifest(wrong_type),
    "source_bytes has incompatible type",
    class = "lna_error_manifest"
  )

  wrong_unit <- manifest
  wrong_unit$schema$columns$repetition_time_seconds$unit <- "millisecond"
  expect_error(
    validate_lna_manifest(wrong_unit),
    "repetition_time_seconds has incompatible unit",
    class = "lna_error_manifest"
  )

  forbidden_null <- manifest
  forbidden_null$data$record_id <- NA_character_
  expect_error(
    validate_lna_manifest(forbidden_null),
    "record_id\\[1\\] may not be null",
    class = "lna_error_manifest"
  )

  bad_number <- manifest
  bad_number$data$n_volumes <- 2.5
  expect_error(
    validate_lna_manifest(bad_number),
    "n_volumes\\[1\\] must be an integer",
    class = "lna_error_manifest"
  )

  incomplete <- manifest
  incomplete$data$status <- "complete"
  expect_error(
    validate_lna_manifest(incomplete),
    "complete rows require lna_checksum and lna_bytes",
    class = "lna_error_manifest"
  )
})

test_that("relative-path and uniqueness invariants fail closed", {
  unsafe_record <- manifest_test_record(
    "unsafe", list(sub = "01", task = "rest", run = "01"),
    source_path = "../outside/sub-01_task-rest_bold.nii.gz"
  )
  expect_error(
    lna_manifest(manifest_test_plan(unsafe_record)),
    "traversal-safe archive-relative path",
    class = "lna_error_manifest"
  )

  unsafe_association <- list(
    association_id = "unsafe-events",
    role = "events",
    path = "../outside/events.tsv",
    checksum = manifest_test_checksum("unsafe-events")
  )
  unsafe_companion <- manifest_test_record(
    "unsafe-companion", list(sub = "01", task = "rest", run = "01"),
    associations = list(unsafe_association)
  )
  expect_error(
    lna_manifest(manifest_test_plan(unsafe_companion)),
    "traversal-safe archive-relative path",
    class = "lna_error_manifest"
  )

  records <- list(
    manifest_test_record("one", list(sub = "01", task = "rest", run = "01")),
    manifest_test_record("two", list(sub = "02", task = "rest", run = "01"))
  )
  manifest <- lna_manifest(manifest_test_plan(records))

  duplicate_id <- manifest
  duplicate_id$data$record_id[[2L]] <- duplicate_id$data$record_id[[1L]]
  expect_error(
    validate_lna_manifest(duplicate_id), "duplicate record_id",
    class = "lna_error_manifest"
  )

  duplicate_path <- manifest
  duplicate_path$data$lna_path[[2L]] <- duplicate_path$data$lna_path[[1L]]
  expect_error(
    validate_lna_manifest(duplicate_path), "duplicate lna_path",
    class = "lna_error_manifest"
  )

  traversal <- manifest
  traversal$data$lna_path[[1L]] <- "sub-01/../escape.lna.h5"
  expect_error(
    validate_lna_manifest(traversal), "traversal-safe archive-relative path",
    class = "lna_error_manifest"
  )
})

test_that("declared unknown optional columns are preserved", {
  records <- list(
    manifest_test_record("one", list(sub = "01", task = "rest", run = "01")),
    manifest_test_record("two", list(sub = "02", task = "rest", run = "01"))
  )
  extra <- list(future_quality_score = list(
    values = c(0.95, NA_real_),
    definition = list(
      type = "number",
      nullable = TRUE,
      unit = "proportion",
      description = "Future quality score",
      role = "future-extension",
      display_hint = "continuous"
    )
  ))
  manifest <- lna_manifest(manifest_test_plan(records), extra_columns = extra)
  first_root <- withr::local_tempdir()
  second_root <- withr::local_tempdir()

  write_lna_manifest(manifest, first_root)
  restored <- read_lna_manifest(first_root)
  write_lna_manifest(restored, second_root)
  restored_again <- read_lna_manifest(second_root)

  expect_true("future_quality_score" %in% names(restored$data))
  expect_identical(
    restored$schema$columns$future_quality_score$role,
    "future-extension"
  )
  expect_identical(
    restored$schema$columns$future_quality_score$display_hint,
    "continuous"
  )
  expect_equal(restored$data$future_quality_score, c(0.95, NA_real_))
  expect_identical(restored_again$schema, restored$schema)
  expect_equal(restored_again$data, restored$data)
})

test_that("reader validates schema before interpreting rows and rejects new majors", {
  record <- manifest_test_record(
    "versioned", list(sub = "01", task = "rest", run = "01")
  )
  root <- withr::local_tempdir()
  write_lna_manifest(lna_manifest(manifest_test_plan(record)), root)
  schema_path <- file.path(root, "lna_manifest.json")
  tsv_path <- file.path(root, "lna_manifest.tsv")
  schema <- jsonlite::read_json(schema_path, simplifyVector = FALSE)
  schema$manifest_schema_version <- "2.0.0"
  jsonlite::write_json(
    schema, schema_path, auto_unbox = TRUE, null = "null", pretty = TRUE
  )
  writeLines("this is not a valid manifest table", tsv_path)

  error <- expect_error(
    read_lna_manifest(root), class = "lna_error_manifest"
  )
  expect_match(conditionMessage(error), "unsupported manifest schema major version 2")
  expect_identical(error$field, "manifest_schema_version")
})

test_that("reader narrowly migrates the shipped 1.0.0 manifest layout", {
  record <- manifest_test_record(
    "legacy", list(sub = "01", task = "rest", run = "01")
  )
  root <- withr::local_tempdir()
  write_lna_manifest(lna_manifest(manifest_test_plan(record)), root)
  schema_path <- file.path(root, "lna_manifest.json")
  tsv_path <- file.path(root, "lna_manifest.tsv")
  schema <- jsonlite::read_json(schema_path, simplifyVector = FALSE)
  schema$manifest_schema_version <- "1.0.0"
  schema$column_order <- setdiff(
    unlist(schema$column_order, use.names = FALSE), "shared_basis_json"
  )
  schema$columns$shared_basis_json <- NULL
  jsonlite::write_json(
    schema, schema_path, auto_unbox = TRUE, null = "null", digits = NA,
    pretty = TRUE
  )
  raw <- utils::read.delim(
    tsv_path, header = TRUE, sep = "\t", quote = "", comment.char = "",
    colClasses = "character", check.names = FALSE, na.strings = character(),
    stringsAsFactors = FALSE
  )
  raw$shared_basis_json <- NULL
  utils::write.table(
    raw, tsv_path, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = "n/a", fileEncoding = "UTF-8"
  )

  restored <- read_lna_manifest(root)
  expect_identical(restored$schema$manifest_schema_version, "1.0.0")
  expect_identical(restored$data$shared_basis_json[[1L]], list(status = "none"))
  expect_identical(
    names(restored$data), restored$schema$column_order
  )

  unknown_layout <- schema
  unknown_layout$column_order <- setdiff(
    unlist(unknown_layout$column_order, use.names = FALSE), "associations_json"
  )
  unknown_layout$columns$associations_json <- NULL
  jsonlite::write_json(
    unknown_layout, schema_path, auto_unbox = TRUE, null = "null", digits = NA,
    pretty = TRUE
  )
  expect_error(
    read_lna_manifest(root), "missing core column.*associations_json",
    class = "lna_error_manifest"
  )
})

test_that("manifest row and column order is deterministic", {
  records <- list(
    manifest_test_record(
      "z-record", list(sub = "09", task = "rest", run = "02", echo = "02")
    ),
    manifest_test_record(
      "a-record", list(sub = "01", task = "motor", run = "01")
    )
  )
  forward <- lna_manifest(manifest_test_plan(records))
  reverse <- lna_manifest(manifest_test_plan(rev(records)))

  expect_identical(forward$schema$column_order, reverse$schema$column_order)
  expect_identical(forward$schema$columns, reverse$schema$columns)
  expect_equal(forward$data, reverse$data)
  expect_identical(
    forward$data$record_id,
    vapply(manifest_test_plan(records)$artifacts, `[[`, character(1), "record_id")
  )
})

test_that("manifest pair is complete and overwrite is explicit", {
  record <- manifest_test_record(
    "pair", list(sub = "01", task = "rest", run = "01")
  )
  manifest <- lna_manifest(manifest_test_plan(record))
  root <- withr::local_tempdir()
  paths <- write_lna_manifest(manifest, root)

  expect_true(all(file.exists(unlist(paths))))
  expect_error(
    write_lna_manifest(manifest, root), "already exists",
    class = "lna_error_manifest"
  )

  unlink(paths$json)
  expect_error(
    read_lna_manifest(root), "manifest pair is incomplete",
    class = "lna_error_manifest"
  )
})
