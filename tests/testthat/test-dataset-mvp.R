mvp_record <- function(
    source_root, record_id, relative_path, entities, values, affine,
    scope = "raw", pipeline = NULL, reference = NULL, with_events = FALSE) {
  source <- file.path(source_root, relative_path)
  dir.create(dirname(source), recursive = TRUE, showWarnings = FALSE)
  image <- neuroim2::NeuroVec(
    values,
    neuroim2::NeuroSpace(dim(values), trans = affine)
  )
  neuroim2::write_vec(image, source)
  metadata_relative <- sub("[.]nii[.]gz$", ".json", relative_path)
  metadata_source <- file.path(source_root, metadata_relative)
  jsonlite::write_json(
    list(RepetitionTime = 0.8, TaskName = entities$task),
    metadata_source, auto_unbox = TRUE, pretty = TRUE
  )
  associations <- list(list(
    association_id = paste0("metadata:", record_id),
    role = "metadata-source",
    path = metadata_relative,
    checksum = neuroarchive:::.journal_file_checksum(metadata_source),
    metadata = list(
      inheritance_level = "run", precedence = 100L,
      fields = c("RepetitionTime", "TaskName"),
      media_type = "application/json"
    )
  ))
  if (isTRUE(with_events)) {
    event_entities <- entities[intersect(names(entities), c("sub", "ses", "task", "run"))]
    event_stem <- paste(vapply(names(event_entities), function(name) {
      paste0(name, "-", event_entities[[name]])
    }, character(1)), collapse = "_")
    event_relative <- file.path(
      dirname(relative_path), paste0(event_stem, "_events.tsv")
    )
    event_source <- file.path(source_root, event_relative)
    writeLines(
      c("onset\tduration\ttrial_type", "0\t1\tgo"),
      event_source, useBytes = TRUE
    )
    associations[[length(associations) + 1L]] <- list(
      association_id = paste0("events:", record_id),
      role = "events", path = gsub("\\\\", "/", event_relative),
      checksum = neuroarchive:::.journal_file_checksum(event_source),
      metadata = list(media_type = "text/tab-separated-values")
    )
  }
  lna_source_record(
    record_id = record_id,
    source_path = relative_path,
    source_checksum = neuroarchive:::.journal_file_checksum(source),
    entities = entities,
    datatype = "func", suffix = "bold", extension = ".nii.gz",
    scope = scope, pipeline = pipeline,
    resolved_metadata = list(RepetitionTime = 0.8, TaskName = entities$task),
    associations = associations,
    geometry = list(
      dims = dim(values)[1:3], affine = affine, units = "mm",
      voxel_ordering = "first-axis-fastest", n_volumes = dim(values)[[4L]]
    ),
    reference_assertion = reference
  )
}

mvp_fixture <- function(source_root) {
  base_values <- array(seq_len(24), dim = c(2L, 2L, 2L, 3L))
  shifted <- diag(c(2, 2, 2, 1))
  shifted[1, 4] <- 11
  wider_values <- array(seq_len(36) / 3, dim = c(3L, 2L, 2L, 3L))
  template <- function(record_id) {
    lna_reference_assertion(
      record_id = record_id,
      identifier = "MNI152NLin2009cAsym",
      kind = "template", origin = "bids-space-entity",
      evidence = list(space_entity = "MNI152NLin2009cAsym")
    )
  }
  records <- list(
    mvp_record(
      source_root, "raw-echo-1",
      "sub-01/ses-A/func/sub-01_ses-A_task-rest_acq-mb_run-01_echo-1_bold.nii.gz",
      list(sub = "01", ses = "A", task = "rest", acq = "mb", run = "01", echo = "1"),
      base_values, diag(c(2, 2, 2, 1)), with_events = TRUE
    ),
    mvp_record(
      source_root, "raw-echo-2",
      "sub-01/ses-A/func/sub-01_ses-A_task-rest_acq-mb_run-01_echo-2_bold.nii.gz",
      list(sub = "01", ses = "A", task = "rest", acq = "mb", run = "01", echo = "2"),
      base_values + 100, diag(c(2, 2, 2, 1)), with_events = TRUE
    ),
    mvp_record(
      source_root, "raw-shifted",
      "sub-01/ses-A/func/sub-01_ses-A_task-rest_run-02_bold.nii.gz",
      list(sub = "01", ses = "A", task = "rest", run = "02"),
      base_values + 200, shifted, with_events = TRUE
    ),
    mvp_record(
      source_root, "raw-wide",
      "sub-02/func/sub-02_task-nback_run-01_bold.nii.gz",
      list(sub = "02", task = "nback", run = "01"),
      wider_values, diag(c(2, 2, 2, 1))
    )
  )
  derivative_id <- "mni-preproc"
  records[[length(records) + 1L]] <- mvp_record(
    source_root, derivative_id,
    paste0(
      "derivatives/fmriprep/sub-01/ses-A/func/",
      "sub-01_ses-A_task-rest_run-01_space-MNI152NLin2009cAsym_",
      "res-2_desc-preproc_bold.nii.gz"
    ),
    list(
      sub = "01", ses = "A", task = "rest", run = "01",
      space = "MNI152NLin2009cAsym", res = "2", desc = "preproc"
    ),
    base_values + 300, diag(c(2, 2, 2, 1)),
    scope = "derivative", pipeline = "fmriprep",
    reference = template(derivative_id)
  )
  list(
    records = records,
    values = stats::setNames(list(
      base_values, base_values + 100, base_values + 200,
      wider_values, base_values + 300
    ), vapply(records, `[[`, character(1), "record_id"))
  )
}

mvp_downgrade_manifest_1_0 <- function(archive_root) {
  schema_path <- file.path(archive_root, "lna_manifest.json")
  table_path <- file.path(archive_root, "lna_manifest.tsv")
  declaration_path <- file.path(archive_root, "lna_dataset.json")
  schema <- jsonlite::read_json(schema_path, simplifyVector = FALSE)
  schema$manifest_schema_version <- "1.0.0"
  schema$column_order <- setdiff(
    unlist(schema$column_order, use.names = FALSE), "shared_basis_json"
  )
  schema$columns$shared_basis_json <- NULL
  raw <- utils::read.delim(
    table_path, header = TRUE, sep = "\t", quote = "", comment.char = "",
    colClasses = "character", check.names = FALSE, na.strings = character(),
    stringsAsFactors = FALSE
  )
  raw$shared_basis_json <- NULL
  declaration <- jsonlite::read_json(
    declaration_path, simplifyVector = FALSE
  )
  declaration$manifest$schema_version <- "1.0.0"
  jsonlite::write_json(
    schema, schema_path, auto_unbox = TRUE, null = "null", digits = NA,
    pretty = TRUE
  )
  utils::write.table(
    raw, table_path, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = "n/a", fileEncoding = "UTF-8"
  )
  jsonlite::write_json(
    declaration, declaration_path, auto_unbox = TRUE, null = "null",
    digits = NA, pretty = TRUE
  )
  invisible(archive_root)
}

test_that("independent per-file exchange MVP packs, queries, and restores", {
  source_root <- withr::local_tempdir()
  archive_root <- tempfile("neuroarchive_mvp_archive_")
  restored_root <- tempfile("neuroarchive_mvp_restored_")
  fixture <- mvp_fixture(source_root)

  result <- write_lna_dataset(
    fixture$records, archive_root, source_root,
    metadata_policy = lna_metadata_policy(portable = TRUE),
    owner_id = "mvp-fixture"
  )
  expect_length(result$written_record_ids, 5L)
  expect_equal(length(unique(result$manifest$data$lna_path)), 5L)
  expect_true(all(result$manifest$data$status == "complete"))
  expect_true(all(file.exists(file.path(
    archive_root, result$manifest$data$lna_path
  ))))

  manifest <- result$manifest$data
  grid <- stats::setNames(manifest$grid_id, manifest$record_id)
  expect_identical(grid[["raw-echo-1"]], grid[["raw-echo-2"]])
  expect_false(identical(grid[["raw-echo-1"]], grid[["raw-shifted"]]))
  expect_false(identical(grid[["raw-echo-1"]], grid[["raw-wide"]]))
  references <- stats::setNames(manifest$reference_id, manifest$record_id)
  expect_false(identical(references[["raw-echo-1"]], references[["raw-echo-2"]]))

  dataset <- open_lna_dataset(archive_root)
  echo_two <- query_lna(
    dataset, sub = "01", ses = "A", run = "01", echo = "2"
  )
  expect_identical(echo_two$record_id, "raw-echo-2")
  derivative <- query_lna(
    dataset,
    scope = "derivative", pipeline = "fmriprep",
    space = "MNI152NLin2009cAsym", desc = "preproc"
  )
  expect_identical(derivative$record_id, "mni-preproc")

  receipt <- restore_bids(dataset, restored_root, closure = "closed")
  expect_identical(receipt$status, "closed")
  expect_length(receipt$selection$restored_record_ids, 5L)
  for (i in seq_len(nrow(manifest))) {
    row <- manifest[i, , drop = FALSE]
    nifti_relative <- sub("[.]lna[.]h5$", ".nii.gz", row$lna_path[[1L]])
    restored <- file.path(restored_root, nifti_relative)
    source <- file.path(source_root, row$source_path[[1L]])
    expect_true(file.exists(restored), info = row$record_id)
    numerical <- validate_lna_restoration(
      source, restored, lna_error_contract("lossless"),
      record_id = row$record_id
    )
    expect_true(numerical$passed, info = row$record_id)
  }
  event_files <- list.files(
    restored_root, pattern = "events[.]tsv$", recursive = TRUE,
    full.names = TRUE
  )
  expect_true(length(event_files) >= 2L)
  expect_true(all(vapply(event_files, function(path) {
    identical(
      neuroarchive:::.journal_file_checksum(path),
      neuroarchive:::.journal_file_checksum(file.path(
        source_root,
        substring(path, nchar(restored_root) + 2L)
      ))
    )
  }, logical(1))))
})

test_that("MVP fixture remains bidser-free", {
  expect_false("bidser" %in% names(getNamespaceImports("neuroarchive")))
  owned <- c("plan_lna_dataset", "validate_lna_dataset_plan")
  source <- unlist(lapply(owned, function(name) {
    deparse(body(get(name, envir = asNamespace("neuroarchive"))))
  }))
  expect_false(any(grepl("bidser::", source, fixed = TRUE)))
})

test_that("supported manifest versions query, transfer, and restore source-independently", {
  for (version in c("1.1.0", "1.0.0")) {
    source_root <- tempfile(paste0("manifest-source-", version, "-"))
    dir.create(source_root)
    archive_root <- tempfile(paste0("manifest-archive-", version, "-"))
    receiver_root <- tempfile(paste0("manifest-receiver-", version, "-"))
    restored_root <- tempfile(paste0("manifest-restored-", version, "-"))
    values <- array(seq_len(24), dim = c(2L, 2L, 2L, 3L))
    record <- mvp_record(
      source_root, paste0("version-", version),
      paste0("sub-01/func/sub-01_task-rest_run-01_bold.nii.gz"),
      list(sub = "01", task = "rest", run = "01"), values,
      diag(c(2, 2, 2, 1)), with_events = TRUE
    )
    write_lna_dataset(
      record, archive_root, source_root,
      metadata_policy = lna_metadata_policy(portable = TRUE),
      owner_id = paste0("manifest-", version)
    )
    if (identical(version, "1.0.0")) {
      mvp_downgrade_manifest_1_0(archive_root)
    }
    unlink(source_root, recursive = TRUE, force = TRUE)

    archive <- open_lna_dataset(archive_root)
    expect_identical(
      archive$manifest$schema$manifest_schema_version, version,
      info = version
    )
    expect_identical(
      archive$manifest$data$shared_basis_json[[1L]], list(status = "none"),
      info = version
    )
    expect_identical(
      query_lna(archive, sub = "01", run = "01")$record_id,
      paste0("version-", version), info = version
    )
    transferred <- transfer_lna_dataset(archive_root, receiver_root)
    expect_true(transferred$verification$passed, info = version)
    received <- open_lna_dataset(receiver_root)
    receipt <- restore_bids(received, restored_root, closure = "closed")
    expect_identical(receipt$status, "closed", info = version)
    restored_path <- file.path(
      restored_root, "sub-01/func/sub-01_task-rest_run-01_bold.nii.gz"
    )
    expect_true(file.exists(restored_path), info = version)
    restored <- neuroim2::read_vec(restored_path)
    expect_equal(
      as.vector(as.array(restored)), as.vector(values), tolerance = 0,
      info = version
    )
  }
})
