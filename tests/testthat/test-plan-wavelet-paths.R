library(testthat)

# Ensure Plan$add_dataset_def handles wavelet and aux_meta paths

test_that("Plan handles wavelet and aux_meta dataset paths", {
  plan <- Plan$new()
  json <- "{}"
  plan$add_dataset_def("/wavelet/level_0/detail_coefficients", "wavelet_coefficients", "spat.haar_octwave", "run-01", 0L, json, "d1", "eager")
  plan$add_dataset_def("/aux_meta/haar_octwave/valid_blocks_L-1", "aux_meta", "spat.haar_octwave", "run-01", 0L, json, "b1", "eager")
  expect_equal(nrow(plan$datasets), 2)
  expect_true("/wavelet/level_0/detail_coefficients" %in% plan$datasets$path)
  expect_true("/aux_meta/haar_octwave/valid_blocks_L-1" %in% plan$datasets$path)
})
