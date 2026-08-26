run_geometry_descriptor <- function(
    dims = c(2L, 2L, 2L),
    affine = diag(c(2, 2, 2, 1)),
    n_volumes = 3L) {
  list(
    dims = dims,
    transform = affine,
    spatial_units = "mm",
    voxel_ordering = "first-axis-fastest",
    n_volumes = n_volumes
  )
}

run_geometry_file <- function(
    runs,
    geometry,
    transforms = character(),
    policy = "fmri") {
  file <- tempfile(fileext = ".lna.h5")
  write_lna(
    runs,
    file,
    transforms = transforms,
    run_geometry = geometry,
    geometry_policy = policy,
    checksum = "none"
  )
  file
}

test_that("equal-grid runs receive explicit run-scoped geometry", {
  runs <- list(
    `run-01` = array(seq_len(24), dim = c(2, 2, 2, 3)),
    `run-02` = array(seq_len(40), dim = c(2, 2, 2, 5))
  )
  geometry <- list(
    `run-01` = run_geometry_descriptor(n_volumes = 3L),
    `run-02` = run_geometry_descriptor(n_volumes = 5L)
  )
  file <- run_geometry_file(runs, geometry)

  h5 <- hdf5r::H5File$new(file, mode = "r")
  withr::defer(try(h5$close_all(), silent = TRUE))
  expect_true(h5$exists("header/runs/run-01"))
  expect_true(h5$exists("header/runs/run-02"))
  expect_true(h5$exists("header/global"))
  expect_false(h5[["header/global"]]$attr_exists("run_id"))
  expect_false(h5[["header/global"]]$attr_exists("n_volumes"))

  first <- neuroarchive:::.read_lna_header(h5, "run-01")
  second <- neuroarchive:::.read_lna_header(h5, "run-02")
  expect_identical(first$grid_id, second$grid_id)
  expect_identical(first$n_volumes, 3L)
  expect_identical(second$n_volumes, 5L)
  expect_identical(attr(first, "lna_geometry_source"), "run-scoped")
  h5$close_all()

  decoded <- read_lna(file, run_id = "run-*", validate = FALSE)
  expect_identical(decoded[["run-01"]]$meta$header$n_volumes, 3L)
  expect_identical(decoded[["run-02"]]$meta$header$n_volumes, 5L)
  expect_true(validate_lna(file, checksum = FALSE))
})

test_that("equal dimensions with different affines never share a global grid", {
  shifted <- diag(c(2, 2, 2, 1))
  shifted[1:3, 4] <- c(8, -4, 12)
  runs <- list(
    `run-01` = array(seq_len(24), dim = c(2, 2, 2, 3)),
    `run-02` = array(seq_len(24), dim = c(2, 2, 2, 3))
  )
  geometry <- list(
    `run-01` = run_geometry_descriptor(),
    `run-02` = run_geometry_descriptor(affine = shifted)
  )
  file <- run_geometry_file(runs, geometry)

  h5 <- hdf5r::H5File$new(file, mode = "r")
  withr::defer(try(h5$close_all(), silent = TRUE))
  expect_false(h5$exists("header/global"))
  expect_error(
    neuroarchive:::.read_lna_header(h5),
    "run-scoped geometry is ambiguous",
    class = "lna_error_run_geometry"
  )
  h5$close_all()

  first <- as_neurovec(read_lna(file, run_id = "run-01"), file = file)
  second <- as_neurovec(read_lna(file, run_id = "run-02"), file = file)
  expect_equal(neuroim2::trans(neuroim2::space(first)), diag(c(2, 2, 2, 1)))
  expect_equal(neuroim2::trans(neuroim2::space(second)), shifted)
})

test_that("different dimensions round trip independently with fMRI geometry", {
  runs <- list(
    `run-01` = array(seq_len(24), dim = c(2, 2, 2, 3)),
    `run-02` = array(seq_len(54), dim = c(3, 3, 2, 3))
  )
  geometry <- list(
    `run-01` = run_geometry_descriptor(c(2L, 2L, 2L)),
    `run-02` = run_geometry_descriptor(c(3L, 3L, 2L))
  )
  file <- run_geometry_file(runs, geometry, transforms = "quant")
  decoded <- read_lna(file, run_id = "run-*", validate = FALSE)

  expect_equal(lapply(decoded, function(x) dim(x$stash$input)), lapply(runs, dim))
  expect_equal(dim(as_neurovec(decoded[["run-01"]], file = file)), dim(runs[[1L]]))
  expect_equal(dim(as_neurovec(decoded[["run-02"]], file = file)), dim(runs[[2L]]))
})

test_that("fMRI array lists require geometry for every named run", {
  runs <- list(
    `run-01` = array(1, dim = c(2, 2, 2, 3)),
    `run-02` = array(2, dim = c(2, 2, 2, 3))
  )
  file <- local_tempfile(fileext = ".lna.h5")
  error <- expect_error(
    write_lna(
      runs,
      file,
      geometry_policy = "fmri",
      checksum = "none"
    ),
    class = "lna_error_run_geometry"
  )
  expect_setequal(error$runs, c("run-01", "run-02"))
  expect_false(file.exists(file))

  error <- expect_error(
    write_lna(
      runs,
      file,
      geometry_policy = "fmri",
      run_geometry = list(`run-01` = run_geometry_descriptor()),
      checksum = "none"
    ),
    class = "lna_error_run_geometry"
  )
  expect_match(conditionMessage(error), "missing: run-02")
  expect_false(file.exists(file))
})

test_that("global spatial headers are rejected for multiple runs", {
  runs <- list(
    `run-01` = array(1, dim = c(2, 2, 2, 3)),
    `run-02` = array(2, dim = c(2, 2, 2, 3))
  )
  file <- local_tempfile(fileext = ".lna.h5")
  error <- expect_error(
    write_lna(
      runs,
      file,
      header = run_geometry_descriptor(),
      geometry_policy = "fmri",
      checksum = "none"
    ),
    class = "lna_error_run_geometry"
  )
  expect_match(conditionMessage(error), "global spatial header is ambiguous")
  expect_setequal(error$runs, names(runs))
  expect_false(file.exists(file))
})

test_that("legacy single-run global geometry remains readable", {
  run <- array(seq_len(24), dim = c(2, 2, 2, 3))
  file <- run_geometry_file(
    run,
    run_geometry_descriptor(),
    policy = "fmri"
  )
  h5 <- hdf5r::H5File$new(file, mode = "r+")
  h5$link_delete("header/runs")
  if (h5[["header"]]$attr_exists("geometry_policy")) {
    h5[["header"]]$attr_delete("geometry_policy")
  }
  h5$close_all()

  handle <- read_lna(file, validate = FALSE)
  restored <- as_neurovec(handle, file = file)
  expect_equal(dim(restored), dim(run))

  h5 <- hdf5r::H5File$new(file, mode = "r")
  withr::defer(h5$close_all())
  legacy <- neuroarchive:::.read_lna_header(h5, "run-01")
  expect_identical(attr(legacy, "lna_geometry_source"), "legacy-global")
})

test_that("validator rejects missing, malformed, and ambiguous run geometry", {
  runs <- list(
    `run-01` = array(1, dim = c(2, 2, 2, 3)),
    `run-02` = array(2, dim = c(2, 2, 2, 3))
  )
  geometry <- list(
    `run-01` = run_geometry_descriptor(),
    `run-02` = run_geometry_descriptor()
  )

  missing <- run_geometry_file(runs, geometry)
  h5 <- hdf5r::H5File$new(missing, mode = "r+")
  h5$link_delete("header/runs/run-02")
  h5$close_all()
  expect_error(
    validate_lna(missing, checksum = FALSE),
    "run geometry is missing for: run-02",
    class = "lna_error_validation"
  )

  malformed <- run_geometry_file(runs, geometry)
  h5 <- hdf5r::H5File$new(malformed, mode = "r+")
  group <- h5[["header/runs/run-02"]]
  group$attr_delete("transform")
  h5$close_all()
  expect_error(
    validate_lna(malformed, checksum = FALSE),
    "run geometry for run-02 is missing: transform",
    class = "lna_error_validation"
  )

  ambiguous <- run_geometry_file(runs, geometry)
  h5 <- hdf5r::H5File$new(ambiguous, mode = "r+")
  h5$link_delete("header/runs/run-02")
  h5$close_all()
  report <- suppressWarnings(validate_lna(
    ambiguous,
    strict = FALSE,
    checksum = FALSE
  ))
  expect_s3_class(report, "lna_validation_report")
  expect_true("run_geometry_missing" %in% report$issues$code)
})

test_that("aggregate operations reject incompatible runs before finalization", {
  runs <- list(
    `run-01` = array(1, dim = c(2, 2, 2, 3)),
    `run-02` = array(2, dim = c(3, 2, 2, 3))
  )
  geometry <- list(
    `run-01` = run_geometry_descriptor(c(2L, 2L, 2L)),
    `run-02` = run_geometry_descriptor(c(3L, 2L, 2L))
  )
  file <- local_tempfile(fileext = ".lna.h5")
  error <- expect_error(
    write_lna(
      runs,
      file,
      transforms = "aggregate_runs",
      run_geometry = geometry,
      geometry_policy = "fmri",
      checksum = "none"
    ),
    class = "lna_error_run_geometry"
  )
  expect_match(conditionMessage(error), "incompatible runs: run-02")
  expect_identical(error$runs, "run-02")
  expect_false(file.exists(file))
})

test_that("multi-run NeuroObj inputs derive each run's own geometry", {
  first_affine <- diag(c(2, 2, 2, 1))
  second_affine <- first_affine
  second_affine[2, 4] <- 10
  runs <- list(
    `run-01` = neuroim2::DenseNeuroVec(
      array(seq_len(24), dim = c(2, 2, 2, 3)),
      neuroim2::NeuroSpace(c(2, 2, 2, 3), trans = first_affine)
    ),
    `run-02` = neuroim2::DenseNeuroVec(
      array(seq_len(24), dim = c(2, 2, 2, 3)),
      neuroim2::NeuroSpace(c(2, 2, 2, 3), trans = second_affine)
    )
  )
  file <- local_tempfile(fileext = ".lna.h5")
  result <- write_lna(
    runs,
    file,
    geometry_policy = "fmri",
    checksum = "none"
  )

  expect_equal(result$run_geometry[["run-01"]]$transform, first_affine)
  expect_equal(result$run_geometry[["run-02"]]$transform, second_affine)
  decoded <- read_lna(file, run_id = "run-*", validate = FALSE)
  expect_equal(
    neuroim2::trans(neuroim2::space(as_neurovec(decoded[["run-01"]], file = file))),
    first_affine
  )
  expect_equal(
    neuroim2::trans(neuroim2::space(as_neurovec(decoded[["run-02"]], file = file))),
    second_affine
  )
})

test_that("geometry is checked against its own run rather than the first run", {
  runs <- list(
    `run-01` = array(1, dim = c(2, 2, 2, 3)),
    `run-02` = array(2, dim = c(3, 2, 2, 3))
  )
  geometry <- list(
    `run-01` = run_geometry_descriptor(c(2L, 2L, 2L)),
    `run-02` = run_geometry_descriptor(c(2L, 2L, 2L))
  )
  file <- local_tempfile(fileext = ".lna.h5")
  error <- expect_error(
    write_lna(
      runs,
      file,
      run_geometry = geometry,
      geometry_policy = "fmri",
      checksum = "none"
    ),
    class = "lna_error_run_geometry"
  )
  expect_identical(error$field, "run-02.dims")
  expect_identical(error$runs, "run-02")
  expect_false(file.exists(file))

  malformed <- geometry
  malformed[["run-02"]]$transform <- 1:15
  error <- expect_error(
    write_lna(
      runs,
      file,
      run_geometry = malformed,
      geometry_policy = "fmri",
      checksum = "none"
    ),
    class = "lna_error_run_geometry"
  )
  expect_identical(error$field, "run-02.transform")
  expect_false(file.exists(file))
})
