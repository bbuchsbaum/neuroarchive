test_that("wendland_c6 equals legacy wendland_c4 on small mask", {
  skip_on_cran()
  # Small 3D mask
  dims <- c(8L, 8L, 4L)
  arr <- array(TRUE, dims)
  mask <- neuroim2::LogicalNeuroVol(arr, neuroim2::NeuroSpace(dims))

  params_c6 <- list(sigma0 = 3, levels = 1L, radius_factor = 2.5,
                    kernel_type = "wendland_c6", seed = 1L)
  params_c4 <- list(sigma0 = 3, levels = 1L, radius_factor = 2.5,
                    kernel_type = "wendland_c4", seed = 1L)

  B6 <- hrbf_generate_basis(params_c6, mask)
  B4 <- hrbf_generate_basis(params_c4, mask)

  expect_equal(dim(B6), dim(B4))
  # exact equality is plausible at this size; allow small tolerance
  expect_true(max(abs(Matrix::as.matrix(B6 - B4))) < 1e-8)
})

