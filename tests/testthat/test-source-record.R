source_record_hash <- function(letter) {
  paste0("sha256:", strrep(letter, 64L))
}

source_record_geometry <- function(
    dims = c(64L, 64L, 36L),
    affine = diag(c(2, 2, 2, 1)),
    n_volumes = 120L) {
  list(
    dims = dims,
    affine = affine,
    units = "mm",
    voxel_ordering = "first-axis-fastest",
    n_volumes = n_volumes
  )
}

source_record_fixture <- function(
    record_id = "raw-rest-01",
    path = "sub-01/func/sub-01_task-rest_run-01_bold.nii.gz",
    hash_letter = "a",
    entities = list(sub = "01", task = "rest", run = "01"),
    scope = "raw",
    pipeline = NULL,
    associations = list(),
    geometry = source_record_geometry(),
    metadata = list(RepetitionTime = 0.8, TaskName = "rest")) {
  lna_source_record(
    record_id = record_id,
    source_path = path,
    source_checksum = source_record_hash(hash_letter),
    entities = entities,
    datatype = "func",
    suffix = "bold",
    extension = ".nii.gz",
    scope = scope,
    pipeline = pipeline,
    resolved_metadata = metadata,
    associations = associations,
    geometry = geometry
  )
}

source_record_package_root <- function() {
  path <- getwd()
  while (!file.exists(file.path(path, "DESCRIPTION"))) {
    parent <- dirname(path)
    if (identical(parent, path)) {
      path <- system.file(package = "neuroarchive")
      break
    }
    path <- parent
  }
  path
}

test_that("raw source records preserve typed entities and computed geometry", {
  record <- source_record_fixture()

  expect_s3_class(record, "lna_source_record")
  expect_identical(record$schema_version, "1.0.0")
  expect_identical(record$entities$run, "01")
  expect_identical(record$entities$sub, "01")
  expect_identical(record$bids$scope, "raw")
  expect_match(record$geometry$grid_id, "^sha256:[0-9a-f]{64}$")
  expect_identical(
    record$reference_assertion$identifier,
    "unknown-native:raw-rest-01"
  )
})

test_that("list and JSON round trips retain the canonical record", {
  record <- source_record_fixture(
    associations = list(list(
      association_id = "events-01",
      role = "events",
      path = "sub-01/func/sub-01_task-rest_run-01_events.tsv",
      checksum = source_record_hash("b"),
      metadata = list(media_type = "text/tab-separated-values")
    ))
  )

  expect_equal(as_lna_source_record(as.list(record)), record)

  json <- jsonlite::toJSON(
    as.list(record),
    auto_unbox = TRUE,
    null = "null",
    digits = NA,
    matrix = "rowmajor"
  )
  decoded <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  expect_equal(as_lna_source_record(decoded), record)

  shuffled <- as.list(record)[rev(names(as.list(record)))]
  expect_equal(validate_lna_source_records(shuffled), record)
})

test_that("raw, derivative, multi-echo, and acquisition records need no bidser", {
  raw <- source_record_fixture()
  derivative <- source_record_fixture(
    record_id = "fmriprep-rest-01",
    path = paste0(
      "derivatives/fmriprep/sub-01/func/",
      "sub-01_task-rest_run-01_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz"
    ),
    hash_letter = "b",
    entities = list(
      sub = "01", task = "rest", run = "01",
      space = "MNI152NLin2009cAsym", desc = "preproc"
    ),
    scope = "derivative",
    pipeline = "fmriprep"
  )
  echo1 <- source_record_fixture(
    record_id = "echo-01",
    path = "sub-01/func/sub-01_task-rest_run-01_echo-01_bold.nii.gz",
    hash_letter = "c",
    entities = list(sub = "01", task = "rest", run = "01", echo = "01"),
    associations = list(list(
      association_id = "echo-peer-02", role = "echo-peer",
      target_record_id = "echo-02"
    ))
  )
  echo2 <- source_record_fixture(
    record_id = "echo-02",
    path = "sub-01/func/sub-01_task-rest_run-01_echo-02_bold.nii.gz",
    hash_letter = "d",
    entities = list(sub = "01", task = "rest", run = "01", echo = "02"),
    associations = list(list(
      association_id = "echo-peer-01", role = "echo-peer",
      target_record_id = "echo-01"
    ))
  )
  acquisition <- source_record_fixture(
    record_id = "acq-mb8",
    path = "sub-02/func/sub-02_task-rest_acq-mb8_run-01_bold.nii.gz",
    hash_letter = "e",
    entities = list(sub = "02", task = "rest", acq = "mb8", run = "01")
  )

  expect_no_error(validate_lna_source_records(list(raw, derivative)))
  expect_no_error(validate_lna_source_records(list(echo1, echo2)))
  expect_identical(echo1$entities$echo, "01")
  expect_identical(acquisition$entities$acq, "mb8")
  expect_identical(derivative$bids$pipeline, "fmriprep")
  expect_false("bidser" %in% names(getNamespaceImports("neuroarchive")))

  owned <- c(
    "lna_source_record", "as_lna_source_record",
    "validate_lna_source_records"
  )
  source_text <- unlist(lapply(owned, function(name) {
    deparse(body(get(name, envir = asNamespace("neuroarchive"))))
  }))
  expect_false(any(grepl("bidser::|requireNamespace\\(.*bidser", source_text)))
})

test_that("record construction is copy-on-modify and validation is non-mutating", {
  metadata <- list(RepetitionTime = 0.8, Nested = list(Source = "sidecar"))
  record <- source_record_fixture(metadata = metadata)
  original <- as.list(record)

  metadata$Nested$Source <- "changed"
  validated <- validate_lna_source_records(record)

  expect_identical(record$resolved_metadata$Nested$Source, "sidecar")
  expect_identical(as.list(record), original)
  expect_equal(validated, record)
})

test_that("field-local diagnostics reject malformed records", {
  missing_checksum <- function() {
    lna_source_record(
      record_id = "bad-checksum",
      source_path = "sub-01/func/source.nii.gz",
      source_checksum = NULL,
      entities = list(sub = "01"),
      datatype = "func", suffix = "bold", extension = ".nii.gz",
      geometry = source_record_geometry()
    )
  }
  error <- expect_error(missing_checksum(), class = "lna_error_source_record")
  expect_identical(error$field, "source.checksum")

  numeric_entity <- function() {
    source_record_fixture(entities = list(sub = "01", run = 1L))
  }
  error <- expect_error(numeric_entity(), class = "lna_error_source_record")
  expect_identical(error$field, "entities.run")

  bad_geometry <- function() {
    source_record_fixture(geometry = source_record_geometry(dims = c(64, 64)))
  }
  error <- expect_error(bad_geometry(), class = "lna_error_source_record")
  expect_identical(error$field, "geometry.dims")

  bad_metadata <- function() {
    source_record_fixture(metadata = list(private = new.env()))
  }
  error <- expect_error(bad_metadata(), class = "lna_error_source_record")
  expect_identical(error$field, "resolved_metadata.private")

  invalid_scope <- function() {
    source_record_fixture(scope = "derived-ish")
  }
  error <- expect_error(invalid_scope(), class = "lna_error_source_record")
  expect_identical(error$field, "bids.scope")
})

test_that("batch validation rejects duplicate IDs, sources, and open targets", {
  first <- source_record_fixture()
  duplicate_id <- source_record_fixture(
    path = "sub-02/func/sub-02_task-rest_run-01_bold.nii.gz",
    hash_letter = "b",
    entities = list(sub = "02", task = "rest", run = "01")
  )
  expect_error(
    validate_lna_source_records(list(first, duplicate_id)),
    "record_id is duplicated",
    class = "lna_error_source_record"
  )

  duplicate_source <- source_record_fixture(
    record_id = "copy-of-raw-rest-01"
  )
  expect_error(
    validate_lna_source_records(list(first, duplicate_source)),
    "assigned to multiple records",
    class = "lna_error_source_record"
  )

  open_target <- source_record_fixture(
    associations = list(list(
      association_id = "missing-peer", role = "echo-peer",
      target_record_id = "echo-99"
    ))
  )
  error <- expect_error(
    validate_lna_source_records(open_target),
    class = "lna_error_source_record"
  )
  expect_identical(error$field, "associations.target_record_id")
  expect_no_error(validate_lna_source_records(
    open_target,
    require_association_closure = FALSE
  ))
})

test_that("inconsistent associations and asserted geometry fail closed", {
  expect_error(
    source_record_fixture(associations = list(list(
      association_id = "events", role = "events",
      path = "sub-01/func/sub-01_task-rest_run-01_events.tsv"
    ))),
    "requires checksum",
    class = "lna_error_source_record"
  )

  expect_error(
    source_record_fixture(associations = list(
      list(association_id = "peer", role = "echo-peer", target_record_id = "x"),
      list(association_id = "peer", role = "echo-peer", target_record_id = "y")
    )),
    "association_id is duplicated",
    class = "lna_error_source_record"
  )

  geometry <- source_record_geometry()
  geometry$grid_id <- source_record_hash("f")
  error <- expect_error(
    source_record_fixture(geometry = geometry),
    class = "lna_error_source_record"
  )
  expect_identical(error$field, "geometry.grid_id")
})
