dataset_plan_checksum <- function(id) {
  paste0(
    "sha256:",
    digest::digest(charToRaw(id), algo = "sha256", serialize = FALSE)
  )
}

dataset_plan_record <- function(
    record_id,
    entities,
    source_path = NULL,
    scope = "raw",
    pipeline = NULL,
    datatype = "func",
    suffix = "bold",
    associations = list(),
    metadata = list(RepetitionTime = 0.8),
    geometry = NULL) {
  if (is.null(source_path)) {
    source_path <- paste0("source/", record_id, ".nii.gz")
  }
  if (is.null(geometry)) {
    geometry <- list(
      dims = c(4L, 4L, 2L),
      affine = diag(c(2, 2, 2, 1)),
      units = "mm",
      voxel_ordering = "first-axis-fastest",
      n_volumes = 10L
    )
  }
  lna_source_record(
    record_id = record_id,
    source_path = source_path,
    source_checksum = dataset_plan_checksum(record_id),
    entities = entities,
    datatype = datatype,
    suffix = suffix,
    extension = ".nii.gz",
    scope = scope,
    pipeline = pipeline,
    resolved_metadata = metadata,
    associations = associations,
    geometry = geometry
  )
}

test_that("one source record maps to one BIDS-shaped LNA artifact", {
  record <- dataset_plan_record(
    "raw-rest-01",
    list(sub = "01", ses = "baseline", task = "rest", run = "01")
  )
  plan <- plan_lna_dataset(record, "/exchange")

  expect_s3_class(plan, "lna_dataset_plan")
  expect_identical(plan$archive_mode, "one-source-one-lna")
  expect_identical(plan$summary$n_records, 1L)
  expect_identical(plan$summary$n_artifacts, 1L)
  expect_identical(plan$artifacts[[1L]]$record_id, "raw-rest-01")
  expect_identical(
    plan$artifacts[[1L]]$destination_relative,
    paste0(
      "sub-01/ses-baseline/func/",
      "sub-01_ses-baseline_task-rest_run-01_bold.lna.h5"
    )
  )
  expect_identical(
    plan$artifacts[[1L]]$source$checksum,
    record$source$checksum
  )
  expect_match(plan$plan_hash, "^sha256:[0-9a-f]{64}$")
})

test_that("derivative pipeline, space, resolution, and description survive", {
  record <- dataset_plan_record(
    "fmriprep-preproc",
    list(
      sub = "01", ses = "01", task = "rest", run = "01",
      space = "MNI152NLin2009cAsym", res = "2", desc = "preproc"
    ),
    scope = "derivative",
    pipeline = "fmriprep-24.1"
  )
  plan <- plan_lna_dataset(record, "exchange")
  path <- plan$artifacts[[1L]]$destination_relative

  expect_true(startsWith(
    path, "derivatives/fmriprep-24.1/sub-01/ses-01/func/"
  ))
  expect_match(path, "space-MNI152NLin2009cAsym", fixed = TRUE)
  expect_match(path, "res-2", fixed = TRUE)
  expect_match(path, "desc-preproc", fixed = TRUE)
  expect_identical(plan$artifacts[[1L]]$bids$scope, "derivative")
  expect_identical(plan$artifacts[[1L]]$bids$pipeline, "fmriprep-24.1")
})

test_that("exact entity values prevent prefix over-selection and collisions", {
  records <- list(
    dataset_plan_record(
      "task-rest", list(sub = "01", task = "rest", run = "01")
    ),
    dataset_plan_record(
      "task-resting", list(sub = "01", task = "resting", run = "01")
    ),
    dataset_plan_record(
      "space-mni", list(
        sub = "01", task = "rest", run = "02", space = "MNI"
      )
    ),
    dataset_plan_record(
      "space-mni-long", list(
        sub = "01", task = "rest", run = "02",
        space = "MNI152NLin2009cAsym"
      )
    )
  )
  plan <- plan_lna_dataset(records, "exchange")
  paths <- vapply(
    plan$artifacts, `[[`, character(1), "destination_relative"
  )

  expect_length(unique(paths), 4L)
  expect_true(any(grepl("task-rest_", paths, fixed = TRUE)))
  expect_true(any(grepl("task-resting_", paths, fixed = TRUE)))
  expect_true(any(grepl("space-MNI_", paths, fixed = TRUE)))
  expect_true(any(grepl("space-MNI152NLin2009cAsym_", paths, fixed = TRUE)))
})

test_that("echo, acquisition, direction, reconstruction, and extra entities survive", {
  record <- dataset_plan_record(
    "multi-entity",
    list(
      sub = "07", ses = "02", task = "motor", acq = "mb8", rec = "moco",
      dir = "AP", run = "03", echo = "02", part = "mag", den = "91k",
      custom = "v1"
    )
  )
  path <- plan_lna_dataset(record, "exchange")$artifacts[[1L]]$destination_relative
  expected <- c(
    "sub-07", "ses-02", "task-motor", "acq-mb8", "rec-moco",
    "dir-AP", "run-03", "echo-02", "part-mag", "den-91k", "custom-v1",
    "bold.lna.h5"
  )
  for (component in expected) expect_match(path, component, fixed = TRUE)
})

test_that("datatype, suffix, and scope remain path-significant", {
  base_entities <- list(sub = "01", task = "rest", run = "01")
  records <- list(
    dataset_plan_record(
      "raw-bold", base_entities, datatype = "func", suffix = "bold"
    ),
    dataset_plan_record(
      "raw-events", base_entities, datatype = "func", suffix = "events"
    ),
    dataset_plan_record(
      "raw-perfusion", base_entities, datatype = "perf", suffix = "asl"
    ),
    dataset_plan_record(
      "derivative-bold", base_entities, scope = "derivative",
      pipeline = "prep", datatype = "func", suffix = "bold"
    )
  )
  paths <- vapply(
    plan_lna_dataset(records, "exchange")$artifacts,
    `[[`, character(1), "destination_relative"
  )

  expect_length(unique(paths), 4L)
  expect_true(any(endsWith(paths, "_events.lna.h5")))
  expect_true(any(grepl("/perf/", paths, fixed = TRUE)))
  expect_true(any(endsWith(paths, "_asl.lna.h5")))
  expect_true(any(startsWith(paths, "derivatives/prep/")))
})

test_that("destination collisions fail before filesystem mutation", {
  first <- dataset_plan_record(
    "collision-a",
    list(sub = "01", task = "rest", run = "01")
  )
  second <- dataset_plan_record(
    "collision-b",
    list(sub = "01", task = "rest", run = "01")
  )
  root <- tempfile("uncreated-plan-root-")

  error <- expect_error(
    plan_lna_dataset(list(first, second), root),
    class = "lna_error_dataset_plan"
  )
  expect_match(conditionMessage(error), "destination collision before writing")
  expect_setequal(error$record_ids, c("collision-a", "collision-b"))
  expect_false(dir.exists(root))
})

test_that("plan order and hash are invariant to input ordering", {
  records <- list(
    dataset_plan_record(
      "order-02", list(sub = "02", task = "rest", run = "02")
    ),
    dataset_plan_record(
      "order-01", list(sub = "01", task = "rest", run = "01")
    ),
    dataset_plan_record(
      "order-03", list(sub = "03", task = "motor", run = "01")
    )
  )
  forward <- plan_lna_dataset(records, "exchange")
  reverse <- plan_lna_dataset(rev(records), "exchange")

  expect_identical(forward$plan_hash, reverse$plan_hash)
  expect_equal(forward$artifacts, reverse$artifacts)
  expect_identical(
    vapply(forward$artifacts, `[[`, character(1), "record_id"),
    c("order-01", "order-02", "order-03")
  )
})

test_that("generated entity combinations have a no-overwrite property", {
  combinations <- expand.grid(
    task = c("rest", "resting"),
    run = c("01", "02"),
    echo = c("01", "02"),
    acq = c("mb4", "mb8"),
    dir = c("AP", "PA"),
    stringsAsFactors = FALSE
  )
  records <- lapply(seq_len(nrow(combinations)), function(i) {
    row <- combinations[i, ]
    dataset_plan_record(
      paste0("generated-", sprintf("%02d", i)),
      as.list(c(sub = "01", row))
    )
  })
  plan <- plan_lna_dataset(records, "exchange")
  paths <- vapply(
    plan$artifacts, `[[`, character(1), "destination_relative"
  )

  expect_length(paths, nrow(combinations))
  expect_length(unique(paths), nrow(combinations))
  expect_identical(plan$summary$n_artifacts, plan$summary$n_records)
})

test_that("plan includes companions, geometry, metadata, and codec profile", {
  association <- list(
    association_id = "events-01",
    role = "events",
    path = "sub-01/func/sub-01_task-rest_run-01_events.tsv",
    checksum = dataset_plan_checksum("events-01"),
    metadata = list(media_type = "text/tab-separated-values")
  )
  record <- dataset_plan_record(
    "complete-entry",
    list(sub = "01", task = "rest", run = "01"),
    associations = list(association),
    metadata = list(RepetitionTime = 0.8, TaskName = "rest")
  )
  codec <- list(
    profile_id = "qpp-analysis-safe-v1",
    loss_mode = "lossy",
    transforms = c("parcel.mean", "qpp"),
    transform_params = list(qpp = list(bits = 12L)),
    error_contract = lna_error_contract(
      "lossy", "float32", max_abs_error = 0.1, rmse_max = 0.05
    )
  )
  plan <- plan_lna_dataset(record, "exchange", codec)
  entry <- plan$artifacts[[1L]]

  expect_identical(entry$codec_profile_id, "qpp-analysis-safe-v1")
  expect_identical(plan$codec_profile$loss_mode, "lossy")
  expect_identical(entry$geometry$grid_id, record$geometry$grid_id)
  expect_identical(entry$resolved_metadata$TaskName, "rest")
  expect_length(entry$expected_companions, 1L)
  expect_identical(
    entry$expected_companions[[1L]]$destination_relative,
    association$path
  )

  table <- as.data.frame(plan)
  expect_identical(table$record_id, "complete-entry")
  expect_identical(table$codec_profile_id, "qpp-analysis-safe-v1")
  expect_output(print(plan), "<lna_dataset_plan>")
})

test_that("planner is bidser-independent and rejects unsafe layout components", {
  record <- dataset_plan_record(
    "bidser-free",
    list(sub = "01", task = "rest", run = "01")
  )
  expect_no_error(plan_lna_dataset(record, "exchange"))
  expect_false("bidser" %in% names(getNamespaceImports("neuroarchive")))

  unsafe <- dataset_plan_record(
    "unsafe-pipeline",
    list(sub = "01", task = "rest", run = "01"),
    scope = "derivative",
    pipeline = "bad/pipeline"
  )
  error <- expect_error(
    plan_lna_dataset(unsafe, "exchange"),
    class = "lna_error_dataset_plan"
  )
  expect_identical(error$field, "bids.pipeline")
})
