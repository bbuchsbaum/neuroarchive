execution_test_checksum <- function(id) {
  paste0(
    "sha256:",
    digest::digest(charToRaw(id), algo = "sha256", serialize = FALSE)
  )
}

execution_test_record <- function(
    record_id,
    sub,
    task = "rest",
    run = "01",
    checksum_seed = record_id) {
  lna_source_record(
    record_id = record_id,
    source_path = paste0("sub-", sub, "/func/", record_id, ".nii.gz"),
    source_checksum = execution_test_checksum(checksum_seed),
    entities = list(sub = sub, task = task, run = run),
    datatype = "func",
    suffix = "bold",
    extension = ".nii.gz",
    resolved_metadata = list(RepetitionTime = 0.8),
    geometry = list(
      dims = c(4L, 4L, 3L),
      affine = diag(c(2, 2, 2, 1)),
      units = "mm",
      voxel_ordering = "first-axis-fastest",
      n_volumes = 8L
    )
  )
}

execution_test_records <- function() {
  list(
    execution_test_record("record-a", "01", run = "01"),
    execution_test_record("record-b", "02", run = "02"),
    execution_test_record("record-c", "03", task = "motor", run = "01")
  )
}

execution_node_by_id <- function(plan, id) {
  nodes <- plan$execution$nodes
  nodes[[match(id, vapply(nodes, `[[`, character(1), "node_id"))]]
}

test_that("empty fit phase uses the complete phase machine", {
  plan <- plan_lna_dataset(execution_test_records()[1:2], "archive")
  execution <- plan$execution

  expect_identical(execution$schema_version, "1.0.0")
  expect_identical(
    execution$phase_order,
    c("preflight", "fit", "write", "validate", "finalize")
  )
  expect_length(execution$fit_resources, 0L)
  expect_identical(
    execution_node_by_id(plan, "phase:fit-complete")$depends_on,
    "phase:preflight"
  )
  expect_setequal(
    vapply(
      plan$execution$nodes[
        vapply(plan$execution$nodes, `[[`, character(1), "kind") == "artifact-write"
      ],
      `[[`, character(1), "node_id"
    ),
    c("write:record-a", "write:record-b")
  )
  expect_identical(tail(execution$topological_order, 1L), "phase:finalize")
  expect_no_error(validate_lna_dataset_plan(plan))
})

test_that("state transitions enforce fit and phase barriers", {
  plan <- plan_lna_dataset(execution_test_records()[1:2], "archive")
  state <- lna_execution_state(plan)

  expect_identical(
    vapply(ready_lna_execution_nodes(plan, state), `[[`, character(1), "node_id"),
    "phase:preflight"
  )
  error <- expect_error(
    advance_lna_execution(plan, state, "write:record-a"),
    class = "lna_error_execution"
  )
  expect_match(conditionMessage(error), "incomplete dependencies")
  expect_identical(error$incomplete_dependencies, "phase:fit-complete")

  state <- advance_lna_execution(plan, state, "phase:preflight")
  expect_identical(
    vapply(ready_lna_execution_nodes(plan, state), `[[`, character(1), "node_id"),
    "phase:fit-complete"
  )
  state <- advance_lna_execution(plan, state, "phase:fit-complete")
  ready_writes <- ready_lna_execution_nodes(plan, state, phase = "write")
  expect_setequal(
    vapply(ready_writes, `[[`, character(1), "node_id"),
    c("write:record-a", "write:record-b")
  )
  state <- advance_lna_execution(plan, state, "write:record-a")
  expect_false("phase:write-complete" %in% vapply(
    ready_lna_execution_nodes(plan, state), `[[`, character(1), "node_id"
  ))
  state <- advance_lna_execution(plan, state, "write:record-b")
  state <- advance_lna_execution(plan, state, "phase:write-complete")
  expect_setequal(
    vapply(ready_lna_execution_nodes(plan, state), `[[`, character(1), "node_id"),
    c("validate:record-a", "validate:record-b")
  )
  state <- advance_lna_execution(plan, state, "validate:record-a")
  state <- advance_lna_execution(plan, state, "validate:record-b")
  expect_identical(
    vapply(ready_lna_execution_nodes(plan, state), `[[`, character(1), "node_id"),
    "phase:finalize"
  )
})

test_that("shared fit dependencies fan out explicitly", {
  resources <- list(list(
    resource_id = "common-basis",
    resource_type = "shared-spatial-basis",
    input_record_ids = c("record-b", "record-a"),
    required_by_record_ids = c("record-c", "record-a", "record-b"),
    parameters = list(rank = 12L, fit_scope = "dataset")
  ))
  plan <- plan_lna_dataset(
    execution_test_records(), "archive", fit_resources = resources
  )
  fit <- execution_node_by_id(plan, "fit:common-basis")
  barrier <- execution_node_by_id(plan, "phase:fit-complete")
  writes <- plan$execution$nodes[
    vapply(plan$execution$nodes, `[[`, character(1), "kind") == "artifact-write"
  ]

  expect_identical(fit$depends_on, "phase:preflight")
  expect_setequal(
    barrier$depends_on, c("phase:preflight", "fit:common-basis")
  )
  expect_true(all(vapply(writes, function(node) {
    identical(node$depends_on, "phase:fit-complete") &&
      identical(node$required_fit_resources, "common-basis")
  }, logical(1))))
  state <- advance_lna_execution(
    plan, lna_execution_state(plan), "phase:preflight"
  )
  expect_identical(
    vapply(ready_lna_execution_nodes(plan, state), `[[`, character(1), "node_id"),
    "fit:common-basis"
  )
  expect_error(
    advance_lna_execution(plan, state, "phase:fit-complete"),
    "fit:common-basis",
    class = "lna_error_execution"
  )
})

test_that("missing dependencies and cycles fail closed", {
  plan <- plan_lna_dataset(execution_test_records()[1:2], "archive")
  missing <- plan$execution
  index <- match(
    "write:record-a",
    vapply(missing$nodes, `[[`, character(1), "node_id")
  )
  missing$nodes[[index]]$depends_on <- "fit:does-not-exist"
  expect_error(
    validate_lna_execution(missing, c("record-a", "record-b")),
    "missing dependency node",
    class = "lna_error_execution"
  )

  cyclic <- plan$execution
  index <- match(
    "phase:preflight",
    vapply(cyclic$nodes, `[[`, character(1), "node_id")
  )
  cyclic$nodes[[index]]$depends_on <- "phase:fit-complete"
  expect_error(
    validate_lna_execution(cyclic, c("record-a", "record-b")),
    "dependency cycle",
    class = "lna_error_execution"
  )
})

test_that("topological order is deterministic across source input order", {
  records <- execution_test_records()
  forward <- plan_lna_dataset(records, "archive")
  reverse <- plan_lna_dataset(rev(records), "archive")

  expect_identical(
    forward$execution$topological_order,
    reverse$execution$topological_order
  )
  expect_identical(forward$execution$nodes, reverse$execution$nodes)
  expect_identical(forward$plan_hash, reverse$plan_hash)
})

test_that("plan hash covers source, fit, codec, and layout identities", {
  base_record <- execution_test_record(
    "stable-record", "01", checksum_seed = "source-a"
  )
  changed_source <- execution_test_record(
    "stable-record", "01", checksum_seed = "source-b"
  )
  base <- plan_lna_dataset(base_record, "archive")
  source_plan <- plan_lna_dataset(changed_source, "archive")
  layout_plan <- plan_lna_dataset(base_record, "another-archive")
  codec_plan <- plan_lna_dataset(
    base_record, "archive", codec_profile = list(
      profile_id = "different-codec", loss_mode = "lossy",
      transforms = "quant", transform_params = list(bits = 8L),
      error_contract = lna_error_contract(
        "lossy", "float32", max_abs_error = 0.5, rmse_max = 0.25
      )
    )
  )
  fit_a <- plan_lna_dataset(base_record, "archive", fit_resources = list(list(
    resource_id = "basis", resource_type = "test-basis-fit",
    input_record_ids = "stable-record", parameters = list(rank = 4L)
  )))
  fit_b <- plan_lna_dataset(base_record, "archive", fit_resources = list(list(
    resource_id = "basis", resource_type = "test-basis-fit",
    input_record_ids = "stable-record", parameters = list(rank = 8L)
  )))

  hashes <- c(
    base$plan_hash, source_plan$plan_hash, layout_plan$plan_hash,
    codec_plan$plan_hash, fit_a$plan_hash, fit_b$plan_hash
  )
  expect_length(unique(hashes), length(hashes))
})

test_that("dataset plans serialize, verify, and resume safely", {
  resources <- list(list(
    resource_id = "basis",
    resource_type = "test-basis-fit",
    input_record_ids = c("record-a", "record-b"),
    parameters = list(rank = 6L)
  ))
  plan <- plan_lna_dataset(
    execution_test_records()[1:2], "archive", fit_resources = resources
  )
  path <- file.path(withr::local_tempdir(), "dataset-plan.json")
  write_lna_dataset_plan(plan, path)
  restored <- read_lna_dataset_plan(path)

  expect_s3_class(restored, "lna_dataset_plan")
  expect_identical(restored$plan_hash, plan$plan_hash)
  expect_identical(restored$execution, plan$execution)
  expect_identical(
    vapply(restored$artifacts, `[[`, character(1), "destination_relative"),
    vapply(plan$artifacts, `[[`, character(1), "destination_relative")
  )

  state <- lna_execution_state(restored, completed_nodes = "phase:preflight")
  expect_identical(state$plan_hash, plan$plan_hash)
  expect_identical(state$completed_nodes, "phase:preflight")

  decoded <- jsonlite::read_json(path, simplifyVector = FALSE)
  decoded$plan_hash <- execution_test_checksum("tampered-plan")
  jsonlite::write_json(
    decoded, path, auto_unbox = TRUE, null = "null", digits = NA,
    pretty = TRUE
  )
  expect_error(
    read_lna_dataset_plan(path), "plan_hash does not match",
    class = "lna_error_execution"
  )
})

test_that("resume state must be dependency-closed and plan-specific", {
  first <- plan_lna_dataset(execution_test_records()[1:2], "archive")
  second <- plan_lna_dataset(execution_test_records()[1:2], "other-archive")

  expect_error(
    lna_execution_state(first, completed_nodes = "write:record-a"),
    "lacks completed dependencies",
    class = "lna_error_execution"
  )
  state <- lna_execution_state(first)
  expect_error(
    ready_lna_execution_nodes(second, state),
    "different plan",
    class = "lna_error_execution"
  )
})
