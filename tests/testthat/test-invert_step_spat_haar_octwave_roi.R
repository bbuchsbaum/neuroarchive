library(testthat)
library(hdf5r)

# Tests for ROI streaming logic in invert_step.spat.haar_octwave (HWT-S2-4)

test_that("invert_step.spat.haar_octwave loads ROI subset using valid blocks", {
  mask <- array(TRUE, dim = c(2,2,2))
  roi <- array(FALSE, dim = c(2,2,2))
  roi[1,1,1] <- TRUE
  root <- matrix(1, nrow = 2, ncol = 1)
  detail <- matrix(seq_len(16), nrow = 2)
  valid <- neuroarchive:::get_valid_finest_blocks(mask)

  tmp <- local_tempfile(fileext = ".h5")
  h5 <- neuroarchive:::open_h5(tmp, mode = "w")
  h5[["/wavelet/level_ROOT/coefficients"]] <- root
  h5[["/wavelet/level_0/detail_coefficients"]] <- detail
  h5[["/aux_meta/haar_octwave/valid_blocks_L-1"]] <- valid
  h5$close_all()

  h5 <- neuroarchive:::open_h5(tmp, mode = "r")
  handle <- DataHandle$new(h5 = h5,
                           subset = list(roi_mask = roi),
                           mask_info = list(mask = mask))
  desc <- list(
    type = "spat.haar_octwave",
    params = list(levels = 1, valid_finest_blocks_path = "/aux_meta/haar_octwave/valid_blocks_L-1"),
    datasets = list(
      list(path = "/wavelet/level_ROOT/coefficients"),
      list(path = "/wavelet/level_0/detail_coefficients"),
      list(path = "/aux_meta/haar_octwave/valid_blocks_L-1")
    ),
    outputs = "wavelet_coefficients",
    inputs = "input_dense_mat"
  )

  called_subset <- FALSE
  local_mocked_bindings(
    h5_read_subset = function(h5g, path, index) {
      called_subset <<- TRUE
      h5g[[path]]$read(args = index)
    },
    .env = asNamespace("neuroarchive")
  )

  out <- neuroarchive:::invert_step.spat.haar_octwave("spat.haar_octwave", desc, handle)
  expect_true(called_subset)
  expect_s3_class(out, "DataHandle")

  neuroarchive:::close_h5_safely(h5)
})
