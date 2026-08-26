test_that("poisson_disk_sample_neuroim2 does not touch R RNG state", {
  skip_on_cran()
  dims <- c(10L, 10L, 6L)
  arr <- array(TRUE, dims)
  mask <- neuroim2::LogicalNeuroVol(arr, neuroim2::NeuroSpace(dims))

  set.seed(12345)
  seed_before <- .Random.seed
  # Call with a fixed seed; internal C++ sampler uses its own RNG
  out <- poisson_disk_sample_neuroim2(mask, radius_mm = 3, seed = 42)
  expect_true(is.matrix(out) && ncol(out) == 3)
  expect_identical(.Random.seed, seed_before)
})

