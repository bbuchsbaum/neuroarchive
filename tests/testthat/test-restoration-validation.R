numerical_geometry <- function(dims, affine, n_volumes) {
  list(
    dims = dims,
    transform = affine,
    spatial_units = "mm",
    voxel_ordering = "first-axis-fastest",
    n_volumes = n_volumes
  )
}

numerical_lna <- function(
    x,
    affine = diag(c(2, 2, 2, 1)),
    transforms = character(),
    transform_params = list(),
    mask = NULL) {
  file <- tempfile(fileext = ".lna.h5")
  write_lna(
    x, file, transforms = transforms, transform_params = transform_params,
    mask = mask,
    run_id = "run-01",
    run_geometry = numerical_geometry(dim(x)[1:3], affine, dim(x)[4]),
    geometry_policy = "fmri", checksum = "none"
  )
  file
}

test_that("lossless integer-representable and floating fixtures decode exactly", {
  affine <- diag(c(2, 2, 2, 1))
  affine[1, 2] <- 0.25
  affine[2, 3] <- -0.125
  affine[1:3, 4] <- c(-12, 8, 24)
  fixtures <- list(
    integer_representable = array(seq_len(96), dim = c(4L, 3L, 2L, 4L)),
    floating = array(
      as.double(seq_len(96)) / 8, dim = c(4L, 3L, 2L, 4L)
    )
  )
  for (name in names(fixtures)) {
    reference <- fixtures[[name]]
    contract <- lna_error_contract(
      "lossless",
      output_dtype = if (is.integer(reference)) "integer32" else "float64"
    )
    file <- numerical_lna(reference, affine)
    restored <- as_neurovec(file, output_dtype = "float64")
    receipt <- validate_lna_restoration(
      reference, restored, contract,
      reference_affine = affine, record_id = name
    )

    expect_true(receipt$passed, info = name)
    expect_true(receipt$metrics$value_exact, info = name)
    expect_identical(receipt$metrics$max_abs_error, 0, info = name)
    expect_identical(receipt$metrics$rmse, 0, info = name)
    expect_identical(receipt$reference$dims, dim(reference), info = name)
    expect_identical(receipt$restored$dims[[4L]], dim(reference)[[4L]], info = name)
    expect_identical(
      receipt$restored$dtype,
      if (is.integer(reference)) "integer32" else "float64",
      info = name
    )
    expect_identical(receipt$metrics$affine_max_abs_error, 0, info = name)
    expect_false(receipt$serialization$source_and_nifti_bytes_compared)
    expect_output(print(receipt), "<lna_numerical_receipt>")
  }
})

test_that("lossy codec plans require dtype-specific predeclared bounds", {
  source <- tempfile(fileext = ".nii.gz")
  writeLines("source", source)
  record <- lna_source_record(
    record_id = "lossy",
    source_path = source,
    source_checksum = neuroarchive:::.journal_file_checksum(source),
    entities = list(sub = "01", task = "rest", run = "01"),
    datatype = "func", suffix = "bold", extension = ".nii.gz",
    geometry = list(
      dims = c(2L, 2L, 2L), affine = diag(4), units = "mm",
      voxel_ordering = "first-axis-fastest", n_volumes = 3L
    )
  )
  expect_error(
    plan_lna_dataset(record, "archive", codec_profile = list(
      profile_id = "lossy-unbounded", loss_mode = "lossy",
      transforms = "quant", transform_params = list(bits = 8L)
    )),
    "require an error_contract",
    class = "lna_error_dataset_plan"
  )
  contract <- lna_error_contract(
    "lossy", "float32", max_abs_error = 0.02, rmse_max = 0.01
  )
  plan <- plan_lna_dataset(record, "archive", codec_profile = list(
    profile_id = "lossy-bounded", loss_mode = "lossy",
    transforms = "quant", transform_params = list(bits = 8L),
    error_contract = contract
  ))

  expect_identical(
    plan$codec_profile$error_contract$output_dtype, "float32"
  )
  expect_identical(plan$codec_profile$error_contract$max_abs_error, 0.02)
})

test_that("quantized and adversarial values obey or fail declared bounds", {
  set.seed(42)
  reference <- array(
    c(-1e3, -1, -0.01, 0, 0.01, 1, 1e3, rnorm(120)),
    dim = c(4L, 4L, 2L, 4L)
  )
  file <- numerical_lna(
    reference,
    transforms = "quant",
    transform_params = list(quant = list(bits = 12L))
  )
  restored <- as.array(as_neurovec(file, output_dtype = "float64"))
  observed_max <- max(abs(restored - reference))
  observed_rmse <- sqrt(mean((restored - reference)^2))
  contract <- lna_error_contract(
    "lossy", "float64",
    max_abs_error = observed_max + 1e-12,
    rmse_max = observed_rmse + 1e-12
  )
  receipt <- validate_lna_restoration(
    reference, restored, contract,
    reference_affine = diag(c(2, 2, 2, 1)),
    restored_affine = diag(c(2, 2, 2, 1))
  )

  expect_true(receipt$passed)
  expect_equal(receipt$metrics$max_abs_error, observed_max)
  expect_equal(receipt$metrics$rmse, observed_rmse)

  corrupted <- restored
  corrupted[1L] <- corrupted[1L] + unname(contract$max_abs_error) + 1
  error <- expect_error(
    validate_lna_restoration(
      reference, corrupted, contract,
      reference_affine = diag(c(2, 2, 2, 1)),
      restored_affine = diag(c(2, 2, 2, 1))
    ),
    "max_abs_error",
    class = "lna_error_numerical_restoration"
  )
  expect_false(error$receipt$passed)
  expect_false(error$receipt$checks$max_abs_error)
})

test_that("dense, masked, eager, lazy, ROI, and time-subset estimands agree", {
  set.seed(7)
  reference <- array(rnorm(3L * 3L * 2L * 6L), dim = c(3L, 3L, 2L, 6L))
  mask <- array(rep(c(TRUE, FALSE), 9L), dim = c(3L, 3L, 2L))
  masked_file <- numerical_lna(reference, mask = mask)
  masked_eager <- as_neurovec(masked_file, output_dtype = "float64")
  file <- numerical_lna(
    reference, transforms = "quant",
    transform_params = list(quant = list(bits = 16L))
  )
  eager <- as_neurovec(file, output_dtype = "float64")
  lazy <- read_lna(file, lazy = TRUE, output_dtype = "float64")
  lazy_vec <- lazy$as_neurovec()
  lazy$close()
  roi <- array(FALSE, dim = dim(mask))
  roi[1, 1, 1] <- TRUE
  roi[3, 3, 1] <- TRUE
  roi_vec <- as_neurovec(
    read_lna(file, roi_mask = roi, output_dtype = "float64"), file = file
  )
  time_idx <- 2:4
  series_reader <- read_lna(file, lazy = TRUE, output_dtype = "float64")
  on.exit(series_reader$close(), add = TRUE)
  lazy_series <- series_reader$series(roi, time_idx = time_idx)
  expected_series <- matrix(as.array(eager), nrow = prod(dim(reference)[1:3]))[
    as.vector(roi), time_idx, drop = FALSE
  ]
  expected_series <- t(expected_series)
  contract <- lna_error_contract("lossless", "float64")
  receipt <- validate_lna_restoration(
    reference, masked_eager, contract, support = mask,
    reference_affine = diag(c(2, 2, 2, 1))
  )

  expect_true(receipt$passed)
  expect_equal(as.array(eager), as.array(lazy_vec))
  expect_equal(
    neuroim2::series(roi_vec, 1, 1, 1),
    neuroim2::series(eager, 1, 1, 1)
  )
  expect_equal(
    neuroim2::series(roi_vec, 3, 3, 1),
    neuroim2::series(eager, 3, 3, 1)
  )
  expect_equal(lazy_series, expected_series)
  expected_support <- neuroarchive:::.lna_support_id(
    neuroarchive:::.lna_grid_id(
      dim(reference)[1:3], diag(c(2, 2, 2, 1)), "mm",
      "first-axis-fastest"
    ),
    mask
  )
  expect_identical(receipt$support_id, expected_support)
})

test_that("rotated/sheared affine survives LNA and NIfTI value validation", {
  affine <- diag(4)
  affine[1:3, 1:3] <- matrix(
    c(0, -2, 0.25, 2, 0, 0, 0, 0, 3), nrow = 3, byrow = TRUE
  )
  affine[1:3, 4] <- c(-10, 6, 20)
  reference_values <- array(
    as.double(seq_len(72)) / 4, dim = c(3L, 3L, 2L, 4L)
  )
  reference <- neuroim2::DenseNeuroVec(
    reference_values,
    neuroim2::NeuroSpace(dim(reference_values), trans = affine)
  )
  file <- numerical_lna(reference_values, affine)
  reconstructed <- as_neurovec(file, output_dtype = "float64")
  nifti <- tempfile(fileext = ".nii.gz")
  neuroim2::write_vec(reconstructed, nifti, data_type = "DOUBLE")
  contract <- lna_error_contract(
    "lossless", "float64", affine_tolerance = 1e-6
  )
  receipt <- validate_lna_restoration(reference, nifti, contract)

  expect_true(receipt$passed)
  expect_true(receipt$metrics$value_exact)
  expect_lte(receipt$metrics$affine_max_abs_error, 1e-6)
  expect_false(identical(
    neuroarchive:::.journal_file_checksum(file),
    neuroarchive:::.journal_file_checksum(nifti)
  ))
  expect_false(receipt$serialization$source_and_nifti_bytes_compared)
})
