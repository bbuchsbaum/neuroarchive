fixed_normalization <- function() {
  list(
    method = "none", fit_scope = "fixed", source_id = "input-values-v1",
    calibration_identity = NULL, fit_id = NULL
  )
}

qpp_contract_descriptor <- function(policy = "fixed",
                                    normalization = fixed_normalization()) {
  asset_id <- paste0("sha256:", strrep("e", 64))
  parcel_binding <- list(
    descriptor = "00_parcel.mean.json", transform = "parcel.mean@1.1.0",
    run_id = "run-01", step_index = 0L, asset_id = asset_id,
    map_checksum = paste0("sha256:", strrep("c", 64)),
    mask_checksum = paste0("sha256:", strrep("d", 64)),
    labels_path = paste0(
      "/spatial/qpp_assets/", sub("^sha256:", "", asset_id),
      "/parcel_labels"
    )
  )
  if (identical(policy, "fixed")) {
    scope <- "fixed"
    profile <- "analysis-safe"
    extra <- list(parcel_binding = parcel_binding)
  } else if (identical(policy, "external-calibration")) {
    scope <- "external-calibration"
    profile <- "analysis-safe"
    extra <- list(
      calibration_id = "calibration-subjects-v1",
      calibration_checksum = paste0("sha256:", strrep("a", 64)),
      parcel_binding = parcel_binding
    )
  } else {
    scope <- "run"
    profile <- "fitted-run"
    extra <- list(fit_id = "run-maxabs-fit-v1")
  }
  params <- c(list(
    quantizer_policy = policy,
    profile_id = profile,
    fit_scope = scope,
    parcel_steps_checksum = paste0("sha256:", strrep("b", 64)),
    normalization_provenance = normalization
  ), extra)
  desc <- list(
    type = "qpp", version = "1.1.0", codec_version = "1.0.0",
    run_id = "run-01", params = params
  )
  desc$analysis_contract <- .qpp_analysis_contract(params)
  desc
}

parcel_contract_descriptor <- function(scope = "fixed",
                                       normalization = fixed_normalization()) {
  params <- list(
    fit_scope = scope,
    source_id = "atlas-v1",
    map_checksum = paste0("sha256:", strrep("c", 64)),
    weight_policy = "uniform",
    normalization_provenance = normalization
  )
  if (scope %in% c("fold", "run", "session", "subject")) {
    params$fit_provenance <- list(fit_id = paste0(scope, "-atlas-fit"))
  }
  if (identical(scope, "external-calibration")) {
    params$calibration_id <- "external-atlas-v1"
    params$calibration_checksum <- paste0("sha256:", strrep("d", 64))
  }
  desc <- list(type = "parcel.mean", version = "1.1.0", params = params)
  desc$analysis_contract <- .parcel_mean_analysis_contract(params)
  desc
}

generic_contract_descriptor <- function(
    type = "temporal.lift53", lossy_support = "pointwise",
    coding_dependency = "none", transformed_units = list(),
    fitted_operations = list(), normalization = fixed_normalization()) {
  contract <- list(
    schema_version = "1.0.0",
    parcellation_fit_scope = "not-applicable",
    quantizer_fit_scope = "not-applicable",
    lossy_support = lossy_support,
    coding_dependency = coding_dependency,
    independence_unit = if (lossy_support == "pointwise") "pointwise" else lossy_support,
    calibration_identity = NULL,
    normalization_provenance = normalization,
    loss_profile = list(kind = "test-profile"),
    transformed_units = transformed_units,
    fitted_operations = fitted_operations
  )
  list(type = type, version = "1.0.0", analysis_contract = contract, params = list())
}

test_that("QPP coding chunks do not constrain statistical folds", {
  desc <- qpp_contract_descriptor("fixed")
  contract <- analysis_contract(desc)
  expect_identical(contract$allowed_split_unit, "pointwise")
  expect_identical(
    contract$coding_dependencies[[1L]]$dependency,
    "chunk"
  )

  folds <- data.frame(
    observation_id = seq_len(8),
    fold = rep(c("A", "B"), 4),
    subject = "sub-01", session = "ses-01", run = "run-01"
  )
  receipt <- assert_split_compatible(desc, folds)
  expect_true(receipt$compatible)
  expect_length(receipt$checked_constraints, 0L)
  expect_identical(
    receipt$ignored_coding_dependencies[[1L]]$dependency,
    "chunk"
  )
})

test_that("run-fitted quantizers require complete runs and cannot masquerade as fixed", {
  desc <- qpp_contract_descriptor("fitted-run-maxabs")
  complete <- data.frame(
    fold = rep(c("A", "B"), each = 4),
    subject = "sub-01", session = "ses-01",
    run = rep(c("run-01", "run-02"), each = 4)
  )
  expect_true(assert_split_compatible(desc, complete)$compatible)
  expect_identical(analysis_contract(desc)$allowed_split_unit, "run")

  split <- complete
  split$fold[2] <- "B"
  expect_error(
    assert_split_compatible(desc, split),
    "qpp@1.1.0.*complete run units.*run-01",
    class = "lna_error_split_incompatible"
  )

  adversarial <- desc
  adversarial$params$fit_scope <- "fixed"
  adversarial$analysis_contract$quantizer_fit_scope <- "fixed"
  expect_error(
    analysis_contract(adversarial),
    "fitted-run-maxabs.*contradicts fit_scope 'fixed'",
    class = "lna_error_analysis_contract"
  )
})

test_that("strict prospective folds reject subject-fitted parcellations", {
  desc <- parcel_contract_descriptor("subject")
  folds <- data.frame(
    fold = rep(c("A", "B"), each = 3),
    subject = rep(c("sub-01", "sub-02"), each = 3),
    session = "ses-01", run = "run-01"
  )
  expect_true(assert_split_compatible(desc, folds)$compatible)
  expect_error(
    assert_split_compatible(desc, folds, strict_prospective = TRUE),
    "strict prospective evaluation forbids same-data fitting",
    class = "lna_error_split_incompatible"
  )
})

test_that("runwise wavelet loss rejects within-run folds and checks exact units", {
  desc <- generic_contract_descriptor(
    lossy_support = "run",
    transformed_units = list(list(scope = "run", ids = c("run-01", "run-02")))
  )
  complete <- data.frame(
    fold = rep(c("A", "B"), each = 4),
    subject = "sub-01", session = "ses-01",
    run = rep(c("run-01", "run-02"), each = 4)
  )
  expect_true(assert_split_compatible(desc, complete)$compatible)

  split <- complete
  split$fold[2] <- "B"
  expect_error(
    assert_split_compatible(desc, split),
    "temporal.lift53@1.0.0.*complete run units.*run-01",
    class = "lna_error_split_incompatible"
  )

  wrong_units <- desc
  wrong_units$analysis_contract$transformed_units[[1L]]$ids <- c(
    "run-01", "run-03"
  )
  expect_error(
    assert_split_compatible(wrong_units, complete),
    "declared transformed units are absent",
    class = "lna_error_split_incompatible"
  )
})

test_that("missing normalization and fitted provenance fail closed", {
  unknown <- qpp_contract_descriptor(
    "fixed", normalization = .analysis_normalization_record(NULL)
  )
  expect_error(
    analysis_contract(unknown),
    "unknown provenance.*normalization",
    class = "lna_error_analysis_contract"
  )

  parcel <- parcel_contract_descriptor("subject")
  parcel$params$fit_provenance <- NULL
  expect_error(
    analysis_contract(parcel),
    "require named fit_provenance",
    class = "lna_error_analysis_contract"
  )

  run_scaled <- qpp_contract_descriptor("fixed", normalization = list(
    method = "percent-signal-change", fit_scope = "run",
    source_id = "whole-run-mean", calibration_identity = NULL,
    fit_id = "run-mean-fit-v1"
  ))
  split <- data.frame(
    fold = c("A", "B"), subject = "sub-01", session = "ses-01",
    run = "run-01"
  )
  expect_error(assert_split_compatible(run_scaled, split), "run-01")
})

test_that("partial-order summaries retain incomparable boundaries", {
  hierarchical <- generic_contract_descriptor(
    lossy_support = "trial-epoch",
    transformed_units = list(list(scope = "trial-epoch", ids = c("trial-1"))),
    fitted_operations = list(
      list(operation_id = "run-fit", fit_scope = "run",
           calibration_identity = NULL, fit_id = "run-fit-v1"),
      list(operation_id = "session-fit", fit_scope = "session",
           calibration_identity = NULL, fit_id = "session-fit-v1")
    )
  )
  expect_identical(analysis_contract(hierarchical)$allowed_split_unit, "session")

  chunked <- generic_contract_descriptor(
    type = "lossy.chunk", lossy_support = "chunk",
    transformed_units = list(list(scope = "chunk", ids = "chunk-1"))
  )
  causal <- generic_contract_descriptor(
    type = "causal.filter", lossy_support = "causal-prefix",
    transformed_units = list(list(scope = "causal-prefix", ids = "prefix-1"))
  )
  expect_identical(
    analysis_contract(list(chunked = chunked, causal = causal))$allowed_split_unit,
    c("chunk", "causal-prefix")
  )
})

test_that("fold-scoped fits require the matching fold-specific fit", {
  desc <- generic_contract_descriptor(
    type = "fold.fitted",
    fitted_operations = list(list(
      operation_id = "fold-basis", fit_scope = "fold",
      calibration_identity = NULL, fit_id = "fold-specific-fit-v1"
    ))
  )
  folds <- data.frame(
    fold = c("A", "A", "B", "B"),
    fit_fold = c("A", "A", "B", "B")
  )
  expect_true(
    assert_split_compatible(desc, folds, strict_prospective = TRUE)$compatible
  )
  folds$fit_fold[4] <- "A"
  expect_error(
    assert_split_compatible(desc, folds),
    "fit_fold must identify the training fit",
    class = "lna_error_split_incompatible"
  )
})

test_that("fmrilatent provenance is conditional rather than inherently unsafe", {
  external <- list(
    operation_id = "fmrilatent-basis",
    fit_scope = "external-calibration",
    calibration_identity = list(
      id = "external-fmrilatent-basis-v1",
      checksum = paste0("sha256:", strrep("e", 64)),
      relationship = "external-disjoint"
    ),
    fit_id = NULL
  )
  safe <- generic_contract_descriptor(
    type = "fmrilatent.temporal", fitted_operations = list(external)
  )
  folds <- data.frame(
    fold = c("A", "B", "A", "B"), subject = "sub-01",
    session = "ses-01", run = "run-01"
  )
  expect_true(
    assert_split_compatible(safe, folds, strict_prospective = TRUE)$compatible
  )

  unsafe <- safe
  unsafe$analysis_contract$fitted_operations[[1L]] <- list(
    operation_id = "fmrilatent-basis", fit_scope = "subject",
    calibration_identity = NULL, fit_id = "same-subject-fit"
  )
  expect_error(
    assert_split_compatible(unsafe, folds, strict_prospective = TRUE),
    "strict prospective",
    class = "lna_error_split_incompatible"
  )
})

test_that("derived contracts survive Plan and HDF5 descriptor round trips", {
  mask <- array(c(TRUE, TRUE), dim = c(2, 1, 1))
  x <- array(c(1, 2, 3, 4, 5, 6), dim = c(2, 1, 1, 3))
  labels <- array(c(10, 20), dim = dim(mask))
  affine <- diag(c(2, 2, 2, 1))
  header <- list(
    dims = dim(mask), spacing = c(2, 2, 2), origin = c(0, 0, 0),
    transform = affine, spatial_units = "mm",
    voxel_ordering = "first-axis-fastest"
  )
  result <- core_write(
    x,
    transforms = c("parcel.mean", "qpp"),
    transform_params = list(
      "parcel.mean" = list(
        labels = labels, fit_scope = "fixed", source_id = "atlas-v1",
        signal_units = "percent-signal-change", spatial_units = "mm",
        voxel_ordering = "first-axis-fastest",
        map_grid_id = .lna_grid_id(
          dim(mask), affine, "mm", "first-axis-fastest"
        ),
        normalization_provenance = fixed_normalization()
      ),
      qpp = list(
        delta = 0.25, quantizer_policy = "fixed",
        units = "percent-signal-change"
      )
    ),
    mask = mask,
    header = header
  )
  expect_identical(analysis_contract(result$plan)$allowed_split_unit, "pointwise")

  path <- local_tempfile(fileext = ".h5")
  h5 <- neuroarchive:::open_h5(path, "w")
  materialise_plan(
    h5, result$plan, run_ids = "run-01",
    header = header
  )
  neuroarchive:::close_h5_safely(h5)
  expect_true(validate_lna(path))
  from_file <- analysis_contract(path)
  expect_identical(from_file$allowed_split_unit, "pointwise")
  expect_identical(
    unname(vapply(from_file$transforms, `[[`, character(1), "transform_id")),
    c("parcel.mean@1.1.0", "qpp@1.1.0")
  )
})

test_that("analysis contract schema is valid JSON", {
  schema <- test_inst_file("schemas", "analysis-contract.schema.json")
  expect_true(jsonlite::validate(paste(readLines(schema, warn = FALSE), collapse = "\n")))
})
