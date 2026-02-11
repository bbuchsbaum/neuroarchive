test_that("label_components treats NA as background", {
  m <- array(FALSE, c(4,4,4))
  m[2,2,2] <- TRUE
  m[3,3,3] <- NA
  res <- neuroarchive:::label_components(m)
  expect_equal(res$count, 1L)
  lab <- res$labels
  expect_equal(lab[2,2,2], lab[2,2,2])
  expect_equal(lab[3,3,3], 0L)
})

