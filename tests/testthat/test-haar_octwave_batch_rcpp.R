test_that("Haar octwave batch Rcpp matches scalar Rcpp path", {
  has_call <- function(sym) {
    dll <- getLoadedDLLs()[["neuroarchive"]]
    if (is.null(dll)) return(FALSE)
    routines <- getDLLRegisteredRoutines(dll)
    !is.null(routines) && sym %in% names(routines$`.Call`)
  }

  skip_if_not(has_call("_neuroarchive_forward_lift_matrix_rcpp"))
  skip_if_not(has_call("_neuroarchive_inverse_lift_matrix_rcpp"))
  skip_if_not(has_call("_neuroarchive_forward_lift_rcpp"))
  skip_if_not(has_call("_neuroarchive_inverse_lift_rcpp"))

  withr::local_options(list(
    lna.hwt.use_rcpp = TRUE,
    lna.hwt.use_rcpp_batch = TRUE
  ))

  set.seed(1)
  mask <- array(TRUE, dim = c(6, 5, 4))
  mask[1, 1, 1] <- FALSE
  T <- 7L
  X <- matrix(rnorm(sum(mask) * T), nrow = T)
  levels <- 3L

  coeff_batch <- neuroarchive:::perform_haar_lift_analysis(X, mask, levels)

  withr::local_options(list(lna.hwt.use_rcpp_batch = FALSE))
  coeff_scalar <- neuroarchive:::perform_haar_lift_analysis(X, mask, levels)

  expect_equal(dim(coeff_batch$root), dim(coeff_scalar$root))
  expect_equal(length(coeff_batch$detail), length(coeff_scalar$detail))
  for (i in seq_along(coeff_batch$detail)) {
    expect_equal(dim(coeff_batch$detail[[i]]), dim(coeff_scalar$detail[[i]]))
  }

  expect_equal(coeff_batch$root, coeff_scalar$root, tolerance = 1e-12)
  for (i in seq_along(coeff_batch$detail)) {
    expect_equal(coeff_batch$detail[[i]], coeff_scalar$detail[[i]], tolerance = 1e-12)
  }

  withr::local_options(list(lna.hwt.use_rcpp_batch = TRUE))
  reco_batch <- neuroarchive:::perform_haar_lift_synthesis(coeff_batch, mask, levels)
  withr::local_options(list(lna.hwt.use_rcpp_batch = FALSE))
  reco_scalar <- neuroarchive:::perform_haar_lift_synthesis(coeff_batch, mask, levels)

  expect_equal(dim(reco_batch), dim(reco_scalar))
  expect_equal(reco_batch, reco_scalar, tolerance = 1e-10)
})

