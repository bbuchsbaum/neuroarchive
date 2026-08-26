shared_dataset_record <- function(
    source_root, record_id, sub, values,
    reference_identifier = "MNI152NLin2009cAsym",
    support_id = NULL) {
  relative <- paste0("sub-", sub, "/func/sub-", sub, "_task-rest_run-01_bold.nii.gz")
  path <- file.path(source_root, relative)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  affine <- diag(c(2, 2, 2, 1))
  neuroim2::write_vec(
    neuroim2::NeuroVec(values, neuroim2::NeuroSpace(dim(values), trans = affine)),
    path
  )
  reference <- lna_reference_assertion(
    record_id = record_id,
    identifier = reference_identifier,
    kind = "template",
    origin = "bids-space-entity",
    evidence = list(space_entity = reference_identifier)
  )
  lna_source_record(
    record_id = record_id,
    source_path = relative,
    source_checksum = neuroarchive:::.journal_file_checksum(path),
    entities = list(sub = sub, task = "rest", run = "01", space = reference_identifier),
    datatype = "func", suffix = "bold", extension = ".nii.gz",
    resolved_metadata = list(), associations = list(),
    geometry = list(
      dims = dim(values)[1:3], affine = affine, units = "mm",
      voxel_ordering = "first-axis-fastest", n_volumes = dim(values)[[4L]],
      support_id = support_id
    ),
    reference_assertion = reference
  )
}

shared_dataset_profile <- function() {
  list(
    profile_id = "lna-shared-basis-pca-v1",
    loss_mode = "lossy",
    transforms = "embed.shared_basis",
    transform_params = list(),
    error_contract = lna_error_contract(
      "lossy", output_dtype = "float64",
      max_abs_error = 1e-8, rmse_max = 1e-8
    )
  )
}

test_that("shared basis grouping requires grid, support, and reference identity", {
  source_root <- withr::local_tempdir()
  values <- array(seq_len(12), dim = c(2L, 2L, 1L, 3L))
  records <- list(
    shared_dataset_record(source_root, "a", "01", values),
    shared_dataset_record(source_root, "b", "02", values + 1),
    shared_dataset_record(
      source_root, "c", "03", values + 2,
      reference_identifier = "MNI152NLin6Asym"
    )
  )
  groups <- lna_shared_basis_groups(records)
  sizes <- sort(vapply(groups, function(group) length(group$record_ids), integer(1)))
  expect_identical(unname(sizes), c(1L, 2L))
  pooled <- groups[[which(vapply(groups, function(group) length(group$record_ids) == 2L, logical(1)))]]
  expect_setequal(pooled$record_ids, c("a", "b"))
  expect_identical(pooled$support_type, "full-grid")

  native <- lapply(c("n1", "n2"), function(id) {
    record <- shared_dataset_record(source_root, id, paste0("0", substr(id, 2L, 2L)), values)
    raw <- as.list(record)
    raw$reference_assertion <- NULL
    do.call(lna_source_record, list(
      record_id = raw$record_id,
      source_path = raw$source$path,
      source_checksum = raw$source$checksum,
      entities = raw$entities,
      datatype = raw$bids$datatype, suffix = raw$bids$suffix,
      extension = raw$bids$extension,
      resolved_metadata = raw$resolved_metadata,
      associations = raw$associations,
      geometry = raw$geometry
    ))
  })
  expect_length(lna_shared_basis_groups(native), 2L)
})

test_that("masked grouping fails closed without masks and validates supplied support", {
  source_root <- withr::local_tempdir()
  values <- array(seq_len(12), dim = c(2L, 2L, 1L, 3L))
  grid <- neuroarchive:::.lna_grid_id(
    c(2L, 2L, 1L), diag(c(2, 2, 2, 1)), "mm", "first-axis-fastest"
  )
  mask <- array(c(TRUE, TRUE, FALSE, TRUE), dim = c(2L, 2L, 1L))
  support <- neuroarchive:::.lna_support_id(grid, mask)
  record <- shared_dataset_record(source_root, "masked", "01", values, support_id = support)
  expect_error(lna_shared_basis_groups(record), "explicit mask")
  groups <- lna_shared_basis_groups(record, support_masks = list(masked = mask))
  expect_identical(groups[[1L]]$support_id, support)
  expect_identical(groups[[1L]]$support_type, "masked")
  wrong <- array(c(TRUE, FALSE, FALSE, TRUE), dim = c(2L, 2L, 1L))
  expect_error(
    lna_shared_basis_groups(record, support_masks = list(masked = wrong)),
    "does not match"
  )
})

test_that("fit pass trains once and plan embeds exact basis provenance", {
  source_root <- withr::local_tempdir()
  archive_root <- tempfile("lna-shared-archive-")
  values <- array(seq_len(12) / 7, dim = c(2L, 2L, 1L, 3L))
  records <- list(
    shared_dataset_record(source_root, "a", "01", values),
    shared_dataset_record(source_root, "b", "02", values + 0.5)
  )
  groups <- lna_shared_basis_groups(records)
  resources <- fit_lna_shared_bases(
    groups, source_root, file.path(archive_root, "bases"), n_atoms = 4L
  )
  expect_length(resources, 1L)
  expect_setequal(resources[[1L]]$input_record_ids, c("a", "b"))
  basis_files <- list.files(
    file.path(archive_root, "bases"), pattern = "[.]lna_basis[.]h5$", full.names = TRUE
  )
  expect_length(basis_files, 1L)

  plan <- plan_lna_dataset(
    records, archive_root, codec_profile = shared_dataset_profile(),
    fit_resources = resources, metadata_policy = lna_metadata_policy(portable = TRUE)
  )
  expect_length(plan$execution$fit_resources, 1L)
  expect_true(all(vapply(plan$artifacts, function(entry) {
    identical(entry$shared_basis$basis_checksum, resources[[1L]]$parameters$basis_checksum)
  }, logical(1))))
  expect_true(all(vapply(plan$artifacts, function(entry) {
    identical(
      entry$codec_transform_params[["embed.shared_basis"]]$basis_id,
      resources[[1L]]$parameters$basis_id
    )
  }, logical(1))))

  result <- write_lna_dataset(plan, source_root = source_root, owner_id = "shared-basis")
  expect_true(all(vapply(result$manifest$data$shared_basis_json, function(x) {
    identical(x$status, "shared") && identical(
      x$basis_checksum, resources[[1L]]$parameters$basis_checksum
    )
  }, logical(1))))
  declaration <- read_lna_dataset_declaration(archive_root)
  expect_length(declaration$fit_resources, 1L)
  expect_identical(
    declaration$fit_resources[[1L]]$parameters$update_policy$repack,
    "reuse-only-by-exact-basis-checksum"
  )

  transfer_root <- tempfile("lna-shared-transfer-")
  transferred <- transfer_lna_dataset(archive_root, transfer_root)
  transfer_paths <- vapply(
    transferred$inventory$files, `[[`, character(1), "path"
  )
  expect_true("bases/registry.json" %in% transfer_paths)
  basis_index <- which(endsWith(transfer_paths, ".lna_basis.h5"))
  expect_length(basis_index, 1L)
  expect_identical(
    unlink(file.path(transfer_root, transfer_paths[[basis_index]])), 0L
  )
  receipt <- jsonlite::read_json(
    file.path(transfer_root, "lna_transfer.json"), simplifyVector = FALSE
  )
  receipt$files <- receipt$files[-basis_index]
  receipt$summary <- list(
    n_files = length(receipt$files),
    total_bytes = sum(vapply(receipt$files, `[[`, numeric(1), "bytes"))
  )
  receipt$inventory_hash <- neuroarchive:::.transfer_inventory_hash(
    receipt$schema_version, receipt$files, receipt$summary
  )
  jsonlite::write_json(
    receipt, file.path(transfer_root, "lna_transfer.json"), auto_unbox = TRUE,
    null = "null", digits = NA, pretty = TRUE
  )
  missing_basis <- verify_lna_transfer(transfer_root, error = FALSE)
  expect_false(missing_basis$passed)
  expect_match(
    missing_basis$failures[["<archive-closure>"]],
    "missing required file"
  )

  first <- result$manifest$data$lna_path[result$manifest$data$record_id == "a"]
  restored <- as_neurovec(file.path(archive_root, first), output_dtype = "float64")
  expect_equal(as.numeric(restored), as.numeric(values), tolerance = 1e-7)
  roi <- array(c(TRUE, FALSE, TRUE, FALSE), dim = c(2L, 2L, 1L))
  subset_handle <- read_lna(
    file.path(archive_root, first), roi_mask = roi, output_dtype = "float64"
  )
  subset_values <- neuroarchive:::as_dense_mat(
    subset_handle$pull_first(c("input_dense_mat", "dense_mat", "input"))$value
  )
  expected_subset <- t(matrix(values, nrow = 4L)[as.logical(roi), , drop = FALSE])
  expect_equal(subset_values, expected_subset, tolerance = 1e-7)
})

test_that("shared basis grid, support, reference, and source reuse mismatches reject pre-write", {
  source_root <- withr::local_tempdir()
  archive_root <- tempfile("lna-shared-mismatch-")
  values <- array(seq_len(12), dim = c(2L, 2L, 1L, 3L))
  records <- list(
    shared_dataset_record(source_root, "a", "01", values),
    shared_dataset_record(source_root, "b", "02", values + 1)
  )
  resource <- fit_lna_shared_bases(
    lna_shared_basis_groups(records), source_root,
    file.path(archive_root, "bases"), n_atoms = 4L
  )[[1L]]
  profile <- shared_dataset_profile()

  bad_grid <- resource
  bad_grid$parameters$grid_id <- paste0("sha256:", paste(rep("0", 64L), collapse = ""))
  expect_error(
    plan_lna_dataset(records, archive_root, profile, fit_resources = list(bad_grid)),
    "grid or support"
  )
  bad_reference <- resource
  bad_reference$parameters$reference_identifier <- "MNI152NLin6Asym"
  expect_error(
    plan_lna_dataset(records, archive_root, profile, fit_resources = list(bad_reference)),
    "reference assertion"
  )
  bad_source <- resource
  bad_source$parameters$source_checksums$a <- paste0(
    "sha256:", paste(rep("f", 64L), collapse = "")
  )
  expect_error(
    plan_lna_dataset(records, archive_root, profile, fit_resources = list(bad_source)),
    "source checksum"
  )
})
