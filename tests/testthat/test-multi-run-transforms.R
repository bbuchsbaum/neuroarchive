library(testthat)
library(neuroarchive)

test_that("raw multi-run storage preserves independent shapes", {
  runs <- list(
    `run-01` = array(seq_len(16), dim = c(2, 2, 1, 4)),
    `run-02` = array(seq_len(60), dim = c(3, 2, 2, 5))
  )
  file <- local_tempfile(fileext = ".h5")

  write_lna(runs, file, transforms = character(), checksum = "none")
  decoded <- read_lna(file, run_id = "run-*", validate = FALSE)

  expect_equal(names(decoded), names(runs))
  for (run_id in names(runs)) {
    expect_equal(decoded[[run_id]]$stash$input, runs[[run_id]])
  }
})

test_that("quant encodes multiple runs independently", {
  runs <- stats::setNames(
    lapply(seq_len(6), function(i) {
      array(
        seq_len(4L * (i + 1L)) * (i + 0.25) + 100 * i,
        dim = c(2, 2, 1, i + 1L)
      )
    }),
    sprintf("run-%02d", seq_len(6))
  )

  multi_file <- local_tempfile(fileext = ".h5")
  result <- write_lna(
    runs,
    multi_file,
    transforms = "quant",
    transform_params = list(quant = list(bits = 12)),
    checksum = "none"
  )
  expect_true(validate_lna(multi_file, checksum = FALSE))

  expect_equal(
    unname(vapply(result$plan$descriptors, `[[`, character(1), "run_id")),
    names(runs)
  )

  h5 <- hdf5r::H5File$new(multi_file, mode = "r")
  withr::defer(try(h5$close_all(), silent = TRUE))
  for (run_id in names(runs)) {
    expect_true(h5$exists(file.path("scans", run_id, "quantized")))
    expect_true(h5$exists(file.path("scans", run_id, "quant_scale")))
    expect_true(h5$exists(file.path("scans", run_id, "quant_offset")))
  }
  h5$close_all()

  decoded <- read_lna(
    multi_file,
    run_id = "run-*",
    validate = FALSE,
    lazy = FALSE
  )
  expect_equal(names(decoded), names(runs))

  lazy_run <- read_lna(
    multi_file,
    run_id = "run-06",
    validate = FALSE,
    lazy = TRUE
  )
  withr::defer(lazy_run$close())
  expect_equal(
    lazy_run$data()$stash$input,
    decoded[["run-06"]]$stash$input,
    tolerance = 1e-12
  )
  lazy_run$close()

  for (run_id in names(runs)) {
    single_file <- local_tempfile(fileext = ".h5")
    write_lna(
      runs[[run_id]],
      single_file,
      transforms = "quant",
      transform_params = list(quant = list(bits = 12)),
      checksum = "none"
    )
    single <- read_lna(single_file, validate = FALSE)

    expect_equal(
      dim(decoded[[run_id]]$stash$input),
      dim(single$stash$input)
    )
    expect_equal(
      length(decoded[[run_id]]$stash$input),
      length(runs[[run_id]])
    )
    expect_equal(
      decoded[[run_id]]$stash$input,
      single$stash$input,
      tolerance = 1e-12
    )
  }
})

test_that("voxel-scoped quantisation keeps blockwise state inside each run", {
  runs <- list(
    `run-01` = array(seq_len(32), dim = c(2, 2, 2, 4)),
    `run-02` = array(seq_len(60) * 3, dim = c(3, 2, 2, 5))
  )
  file <- local_tempfile(fileext = ".h5")

  write_lna(
    runs,
    file,
    transforms = "quant",
    transform_params = list(
      quant = list(bits = 12, scale_scope = "voxel")
    ),
    checksum = "none"
  )

  decoded <- read_lna(file, run_id = "run-*", validate = FALSE)
  expect_equal(names(decoded), names(runs))
  expect_equal(
    lapply(decoded, function(x) dim(x$stash$input)),
    lapply(runs, dim)
  )

  for (run_id in names(runs)) {
    single_file <- local_tempfile(fileext = ".h5")
    write_lna(
      runs[[run_id]],
      single_file,
      transforms = "quant",
      transform_params = list(
        quant = list(bits = 12, scale_scope = "voxel")
      ),
      checksum = "none"
    )
    single <- read_lna(single_file, validate = FALSE)
    expect_equal(
      decoded[[run_id]]$stash$input,
      single$stash$input,
      tolerance = 1e-12
    )
  }
})

test_that("run-local basis assets are linked to the matching embed step", {
  set.seed(42)
  runs <- list(
    `run-01` = array(rnorm(3 * 3 * 1 * 12), dim = c(3, 3, 1, 12)),
    `run-02` = array(rnorm(3 * 3 * 1 * 15) + 20, dim = c(3, 3, 1, 15))
  )
  file <- local_tempfile(fileext = ".h5")

  result <- write_lna(
    runs,
    file,
    transforms = c("basis", "embed", "quant"),
    transform_params = list(
      basis = list(k = 3),
      quant = list(bits = 16)
    ),
    checksum = "none"
  )
  expect_true(validate_lna(file, checksum = FALSE))

  desc <- result$plan$descriptors
  expect_equal(
    unname(vapply(desc, `[[`, character(1), "run_id")),
    rep(names(runs), each = 3L)
  )

  for (run_id in names(runs)) {
    run_desc <- desc[vapply(desc, function(x) identical(x$run_id, run_id), logical(1))]
    basis_desc <- run_desc[vapply(run_desc, function(x) identical(x$type, "basis"), logical(1))][[1]]
    embed_desc <- run_desc[vapply(run_desc, function(x) identical(x$type, "embed"), logical(1))][[1]]
    basis_path <- basis_desc$datasets[[1]]$path

    expect_equal(embed_desc$params$basis_path, basis_path)
    expect_true(startsWith(basis_path, "/basis/"))
  }

  basis_paths <- vapply(
    desc[vapply(desc, function(x) identical(x$type, "basis"), logical(1))],
    function(x) x$datasets[[1]]$path,
    character(1)
  )
  expect_length(unique(basis_paths), 2L)

  decoded <- read_lna(file, run_id = "run-*", validate = FALSE)
  expect_equal(names(decoded), names(runs))
  for (run_id in names(runs)) {
    actual <- as_dense_mat(decoded[[run_id]]$stash$input)
    expected <- as_dense_mat(runs[[run_id]])

    expect_equal(dim(actual), dim(expected))
    expect_true(all(is.finite(actual)))
    expect_equal(colMeans(actual), colMeans(expected), tolerance = 1e-3)
    expect_gt(stats::cor(as.vector(actual), as.vector(expected)), 0.6)
  }
})

test_that("multi-run identities fail closed on collisions", {
  run <- array(seq_len(16), dim = c(2, 2, 1, 4))
  file <- local_tempfile(fileext = ".h5")

  expect_error(
    write_lna(
      list(run, run),
      file,
      transforms = "quant",
      run_id = c("run-01", "run-01"),
      checksum = "none"
    ),
    "run_id values must be unique",
    class = "lna_error_validation"
  )

})
