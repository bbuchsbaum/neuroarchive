library(testthat)

# Test get_roi_detail_indices on a simple ROI within a two-box mask

test_that("get_roi_detail_indices returns expected indices", {
  mask <- make_disjoint_boxes_mask()
  roi <- array(FALSE, dim = dim(mask))
  roi[1:2, 1:2, 1:2] <- TRUE

  out <- neuroarchive:::get_roi_detail_indices(roi, mask, levels = 2)

  expected_level1 <- 1:8
  expected_level2 <- 1:2

  expect_equal(out[[1]], expected_level1)
  expect_equal(out[[2]], expected_level2)
})
