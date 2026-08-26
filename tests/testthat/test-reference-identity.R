reference_hash <- function(letter) {
  paste0("sha256:", strrep(letter, 64L))
}

reference_geometry <- function(
    affine = diag(c(2, 2, 2, 1)),
    support_id = NULL) {
  list(
    dims = c(8L, 8L, 6L),
    affine = affine,
    units = "mm",
    voxel_ordering = "first-axis-fastest",
    n_volumes = 20L,
    support_id = support_id
  )
}

reference_record <- function(
    record_id,
    subject,
    hash_letter,
    assertion = NULL,
    geometry = reference_geometry()) {
  lna_source_record(
    record_id = record_id,
    source_path = paste0(
      "sub-", subject, "/func/sub-", subject,
      "_task-rest_run-01_bold.nii.gz"
    ),
    source_checksum = reference_hash(hash_letter),
    entities = list(sub = subject, task = "rest", run = "01"),
    datatype = "func",
    suffix = "bold",
    extension = ".nii.gz",
    geometry = geometry,
    reference_assertion = assertion
  )
}

bids_template_assertion <- function(record_id, identifier) {
  lna_reference_assertion(
    record_id = record_id,
    identifier = identifier,
    kind = "template",
    origin = "bids-space-entity",
    evidence = list(space_entity = identifier)
  )
}

user_reference_receipt <- function(id = "receipt-01") {
  list(
    receipt_id = id,
    asserted_by = "researcher-01",
    asserted_at = "2026-08-25T15:30:00Z",
    reason = "registration receipt and preprocessing report reviewed"
  )
}

test_that("equal-grid raw records remain reference-incompatible by default", {
  first <- reference_record("native-01", "01", "a")
  second <- reference_record("native-02", "02", "b")

  expect_identical(first$geometry$grid_id, second$geometry$grid_id)
  expect_identical(
    first$reference_assertion$identifier,
    "unknown-native:native-01"
  )
  receipt <- lna_spatial_compatibility(first, second, "grid-only")
  expect_false(receipt$compatible)
  expect_false(receipt$reference$compatible)
  expect_true(receipt$grid$compatible)
  expect_identical(receipt$reference$identity_kind, "provenance-asserted")
  expect_identical(receipt$grid$identity_kind, "computed-verified")
  expect_identical(receipt$grid$verification, "recomputed-from-geometry")
  expect_identical(receipt$reasons, "reference-assertion-mismatch")
})

test_that("unknown-native identity matches only the same record", {
  record <- reference_record("native-01", "01", "a")
  copy <- as_lna_source_record(as.list(record))

  expect_true(lna_spatial_compatibility(record, copy, "grid-only")$compatible)

  forged <- as.list(copy)
  forged$record_id <- "native-02"
  expect_error(
    as_lna_source_record(forged),
    "unknown-native identifier must equal unknown-native:native-02",
    class = "lna_error_source_record"
  )
})

test_that("matching template assertions and grids are compatible", {
  template <- "MNI152NLin2009cAsym"
  first <- reference_record(
    "mni-01", "01", "a",
    assertion = bids_template_assertion("mni-01", template),
    geometry = reference_geometry(support_id = reference_hash("c"))
  )
  second <- reference_record(
    "mni-02", "02", "b",
    assertion = bids_template_assertion("mni-02", template),
    geometry = reference_geometry(support_id = reference_hash("c"))
  )

  receipt <- lna_spatial_compatibility(first, second)
  expect_true(receipt$compatible)
  expect_true(receipt$reference$compatible)
  expect_true(receipt$grid$compatible)
  expect_true(receipt$support$compatible)
  expect_identical(receipt$policy$support, "same-support")
})

test_that("template labels remain incompatible even on the same grid", {
  old_mni <- reference_record(
    "mni6", "01", "a",
    assertion = bids_template_assertion("mni6", "MNI152NLin6Asym")
  )
  new_mni <- reference_record(
    "mni2009", "02", "b",
    assertion = bids_template_assertion("mni2009", "MNI152NLin2009cAsym")
  )

  receipt <- lna_spatial_compatibility(old_mni, new_mni, "grid-only")
  expect_false(receipt$compatible)
  expect_false(receipt$reference$compatible)
  expect_true(receipt$grid$compatible)
})

test_that("BIDS and SpatialReference origins require matching evidence", {
  expect_error(
    lna_reference_assertion(
      record_id = "bad-bids",
      identifier = "MNI152NLin6Asym",
      kind = "template",
      origin = "bids-space-entity",
      evidence = list(space_entity = "MNI152NLin2009cAsym")
    ),
    class = "lna_error_reference_assertion"
  )

  sidecar <- lna_reference_assertion(
    record_id = "sidecar-01",
    identifier = "fsaverage",
    kind = "named-reference",
    origin = "spatial-reference-sidecar",
    evidence = list(
      SpatialReference = "fsaverage",
      sidecar_path = "sub-01_task-rest_bold.json"
    )
  )
  expect_identical(sidecar$origin, "spatial-reference-sidecar")

  expect_error(
    lna_reference_assertion(
      record_id = "contradiction",
      identifier = "unknown-native:contradiction",
      kind = "unknown-native",
      origin = "bids-space-entity",
      evidence = list(space_entity = "unknown-native:contradiction")
    ),
    "known assertion origins cannot use kind='unknown-native'",
    class = "lna_error_reference_assertion"
  )
})

test_that("user overrides require receipts and preserve all prior evidence", {
  original <- reference_record("override-01", "01", "a")
  overridden <- override_lna_reference(
    original,
    identifier = "MNI152NLin2009cAsym",
    receipt = user_reference_receipt(),
    evidence = list(registration_report = "reports/sub-01.html")
  )

  expect_identical(
    original$reference_assertion$identifier,
    "unknown-native:override-01"
  )
  expect_identical(overridden$reference_assertion$status, "overridden")
  expect_identical(overridden$reference_assertion$origin, "user-receipt")
  expect_length(overridden$reference_assertion$history, 1L)
  expect_identical(
    overridden$reference_assertion$history[[1L]]$identifier,
    "unknown-native:override-01"
  )

  second <- override_lna_reference(
    overridden,
    identifier = "MNI152NLin6Asym",
    receipt = user_reference_receipt("receipt-02"),
    evidence = list(registration_report = "reports/revised.html")
  )
  expect_length(second$reference_assertion$history, 2L)
  expect_identical(
    second$reference_assertion$history[[2L]]$receipt$receipt_id,
    "receipt-01"
  )
  expect_identical(
    second$reference_assertion$history[[2L]]$evidence$registration_report,
    "reports/sub-01.html"
  )

  json <- jsonlite::toJSON(
    as.list(second), auto_unbox = TRUE, null = "null", matrix = "rowmajor"
  )
  restored <- as_lna_source_record(
    jsonlite::fromJSON(json, simplifyVector = FALSE)
  )
  expect_equal(restored, second)

  expect_error(
    override_lna_reference(
      original,
      identifier = "MNI152NLin2009cAsym",
      receipt = NULL
    ),
    "require an explicit receipt",
    class = "lna_error_reference_assertion"
  )
})

test_that("support and grid policies fail closed independently", {
  template <- "MNI152NLin2009cAsym"
  first <- reference_record(
    "support-01", "01", "a",
    bids_template_assertion("support-01", template),
    reference_geometry(support_id = reference_hash("c"))
  )
  other_support <- reference_record(
    "support-02", "02", "b",
    bids_template_assertion("support-02", template),
    reference_geometry(support_id = reference_hash("d"))
  )

  strict <- lna_spatial_compatibility(first, other_support)
  expect_false(strict$compatible)
  expect_false(strict$support$identities_match)
  expect_identical(
    strict$reasons,
    "computed-support-missing-or-mismatch"
  )
  expect_true(
    lna_spatial_compatibility(first, other_support, "grid-only")$compatible
  )

  shifted <- diag(c(2, 2, 2, 1))
  shifted[1, 4] <- 4
  other_grid <- reference_record(
    "grid-02", "02", "e",
    bids_template_assertion("grid-02", template),
    reference_geometry(
      affine = shifted,
      support_id = reference_hash("c")
    )
  )
  receipt <- lna_spatial_compatibility(first, other_grid)
  expect_false(receipt$compatible)
  expect_false(receipt$grid$compatible)
  expect_true(receipt$reference$compatible)
})
