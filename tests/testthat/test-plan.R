library(testthat)
library(tibble)

# Assuming Plan class definition is loaded from R/plan.R
# source("../R/plan.R") # Might be needed for interactive testing

test_that("Plan initialization works correctly", {
  plan <- Plan$new(origin_label = "run-01")

  # Check field types
  expect_true(is_tibble(plan$datasets))
  expect_true(is.list(plan$descriptors))
  expect_true(is.list(plan$payloads))
  expect_true(is.integer(plan$next_index))
  expect_true(is.character(plan$origin_label))

  # Check initial values
  expect_equal(nrow(plan$datasets), 0)
  expect_equal(length(plan$descriptors), 0)
  expect_equal(length(plan$payloads), 0)
  expect_equal(plan$next_index, 0L)
  expect_equal(plan$origin_label, "run-01")

  # Check datasets tibble structure
  expected_cols <- c(
    "path", "role", "producer", "origin", "step_index",
    "params_json", "payload_key", "write_mode", "write_mode_effective",
    "dtype"
  )
  expect_equal(names(plan$datasets), expected_cols)
  # Check column types (optional but good)
  expect_type(plan$datasets$path, "character")
  expect_type(plan$datasets$step_index, "integer")

  # Check default origin label
  plan_default <- Plan$new()
  expect_equal(plan_default$origin_label, "global")

  # Check invalid origin label
  expect_error(Plan$new(origin_label = 123))
  expect_error(Plan$new(origin_label = c("a", "b")))
})

test_that("Plan add_payload works correctly", {
  plan <- Plan$new()
  payload1 <- list(a = 1)
  payload2 <- matrix(1:4, 2)

  # Add initial payload
  plan$add_payload("payload1", payload1)
  expect_equal(length(plan$payloads), 1)
  expect_true("payload1" %in% names(plan$payloads))
  expect_identical(plan$payloads$payload1, payload1)

  # Add second payload
  plan$add_payload("payload2", payload2)
  expect_equal(length(plan$payloads), 2)
  expect_true("payload2" %in% names(plan$payloads))
  expect_identical(plan$payloads$payload2, payload2)

  # Check error on duplicate key
  expect_error(
    plan$add_payload("payload1", list(b = 2)),
    "Payload key 'payload1' already exists in plan."
  )

  # Overwrite existing payload
  expect_silent(plan$add_payload("payload1", list(b = 2), overwrite = TRUE))
  expect_identical(plan$payloads$payload1, list(b = 2))

  # Check overwrite argument type
  expect_error(plan$add_payload("x", 1, overwrite = "no"))

  # Check error on invalid key type
  expect_error(plan$add_payload(123, list()))
  expect_error(plan$add_payload(c("a", "b"), list()))
})

test_that("Plan add_dataset_def works correctly", {
  plan <- Plan$new()
  def1 <- list(
    path = "/data/raw",
    role = "input",
    producer = "initial",
    origin = "run-01",
    step_index = -1L,
    params_json = "{}",
    payload_key = "raw_data",
    write_mode = "eager"
  )

  # Add dataset definition
  plan$add_dataset_def(
    path = def1$path, role = def1$role, producer = def1$producer,
    origin = def1$origin, step_index = def1$step_index, params_json = def1$params_json,
    payload_key = def1$payload_key, write_mode = def1$write_mode,
    dtype = NA_character_
  )

  # Check results
  expect_equal(nrow(plan$datasets), 1)
  row1 <- plan$datasets[1, ]
  expect_equal(row1$path, def1$path)
  expect_equal(row1$role, def1$role)
  expect_equal(row1$producer, def1$producer)
  expect_equal(row1$origin, def1$origin)
  expect_equal(row1$step_index, def1$step_index)
  expect_equal(row1$params_json, def1$params_json)
  expect_equal(row1$payload_key, def1$payload_key)
  expect_equal(row1$write_mode, def1$write_mode)
  expect_equal(row1$write_mode_effective, NA_character_)
  expect_equal(row1$dtype, NA_character_)

  # Add another one
  plan$add_dataset_def("/basis/global", "basis", "pca", "global", 0L, '{"k": 50}', "pca_basis", "eager", dtype = NA_character_)
  expect_equal(nrow(plan$datasets), 2)

  # step_index accepts numeric integer
  plan$add_dataset_def("/data/extra", "extra", "dummy", "run-01", 1, "{}", "raw_data", "eager", dtype = NA_character_)
  expect_equal(nrow(plan$datasets), 3)
  expect_equal(plan$datasets$step_index[3], 1L)

  # invalid step_index (non integer numeric)
  expect_error(plan$add_dataset_def("/bad", "data", "dummy", "run-01", 1.5, "{}", "raw_data", "eager", dtype = NA_character_))

  # invalid write_mode
  expect_error(plan$add_dataset_def("/bad", "data", "dummy", "run-01", 0L, "{}", "raw_data", "invalid", dtype = NA_character_))

  # invalid JSON
  expect_error(plan$add_dataset_def("/bad", "data", "dummy", "run-01", 0L, "not json", "raw_data", "eager", dtype = NA_character_))

  # Check some basic type errors handled by stopifnot
  expect_error(plan$add_dataset_def(path=123, role="", producer="", origin="", step_index=0L, params_json="", payload_key="", write_mode="", dtype = NA_character_))
  expect_error(plan$add_dataset_def(path="", role="", producer="", origin="", step_index="a", params_json="", payload_key="", write_mode="", dtype = NA_character_))
})

test_that("Plan add_descriptor and get_next_filename work correctly", {
  plan <- Plan$new()
  desc1 <- list(type = "pca", k = 50)
  desc2 <- list(type = "quant", bits = 8)

  # Check initial filename
  fname1 <- plan$get_next_filename("pca")
  expect_equal(fname1, "00_pca.json")
  expect_equal(plan$next_index, 0L) # get_next_filename should not increment

  # Add first descriptor
  plan$add_descriptor(fname1, desc1)
  expect_equal(length(plan$descriptors), 1)
  expect_true(fname1 %in% names(plan$descriptors))
  expect_identical(plan$descriptors[[fname1]], desc1)
  expect_equal(plan$next_index, 1L) # add_descriptor should increment

  # Check next filename
  fname2 <- plan$get_next_filename("quant")
  expect_equal(fname2, "01_quant.json")
  expect_equal(plan$next_index, 1L) # get_next_filename should not increment

  # Add second descriptor
  plan$add_descriptor(fname2, desc2)
  expect_equal(length(plan$descriptors), 2)
  expect_true(fname2 %in% names(plan$descriptors))
  expect_identical(plan$descriptors[[fname2]], desc2)
  expect_equal(plan$next_index, 2L)

  # Check error on duplicate name
  expect_error(
    plan$add_descriptor(fname1, list()),
    "Descriptor name '00_pca.json' already exists in plan."
  )

  # Check input type errors
  expect_error(plan$add_descriptor(123, list()))
  expect_error(plan$add_descriptor("name", "not_a_list"))
  expect_error(plan$get_next_filename(123))
  expect_error(plan$get_next_filename("../bad"))
  expect_error(plan$get_next_filename("bad/type"))
  expect_error(plan$get_next_filename("bad\\\\type"))
})

test_that("Plan mark_payload_written works correctly", {
  plan <- Plan$new()
  payload1 <- list(a = 1)

  # Add payload
  plan$add_payload("payload1", payload1)
  expect_false(is.null(plan$payloads$payload1))

  # Mark as written
  plan$mark_payload_written("payload1")
  expect_true(is.null(plan$payloads$payload1))

  # Check warning on marking non-existent key
  expect_warning(
    plan$mark_payload_written("non_existent_key"),
    "Payload key 'non_existent_key' not found in plan when trying to mark as written."
  )

  # Check error on invalid key type
  expect_error(plan$mark_payload_written(123))
  expect_error(plan$mark_payload_written(c("a", "b")))
}) 