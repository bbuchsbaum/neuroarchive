test_that("quantize_voxel_slab_rcpp matches quantize_voxel_block_rcpp on extracted blocks", {
  has_call <- function(sym) {
    dll <- getLoadedDLLs()[["neuroarchive"]]
    if (is.null(dll)) return(FALSE)
    routines <- getDLLRegisteredRoutines(dll)
    !is.null(routines) && sym %in% names(routines$`.Call`)
  }

  skip_if_not(has_call("_neuroarchive_quantize_voxel_slab_rcpp"))
  skip_if_not(has_call("_neuroarchive_quantize_voxel_block_rcpp"))

  set.seed(1)
  arr <- array(rnorm(10 * 9 * 6 * 12), dim = c(10, 9, 6, 12))

  # pick an interior slab
  start <- c(3L, 2L, 4L)
  slab_dims <- c(4L, 5L, 2L)
  xi <- start[1]:(start[1] + slab_dims[1] - 1L)
  yi <- start[2]:(start[2] + slab_dims[2] - 1L)
  zi <- start[3]:(start[3] + slab_dims[3] - 1L)

  bits <- 8L
  method <- "range"
  center <- TRUE

  slab_res <- neuroarchive:::quantize_voxel_slab_rcpp(arr, start, slab_dims, bits, method, center)
  block <- arr[xi, yi, zi, , drop = FALSE]
  block_res <- neuroarchive:::quantize_voxel_block_rcpp(block, bits, method, center)

  expect_equal(dim(slab_res$q), dim(block_res$q))
  expect_equal(dim(slab_res$scale), dim(block_res$scale))
  expect_equal(dim(slab_res$offset), dim(block_res$offset))

  expect_equal(slab_res$q, block_res$q)
  expect_equal(slab_res$scale, block_res$scale, tolerance = 1e-12)
  expect_equal(slab_res$offset, block_res$offset, tolerance = 1e-12)
  expect_equal(slab_res$n_clipped_total, block_res$n_clipped_total)
})

test_that("quantize_voxel_slab_fill_rcpp reuses caller buffers and matches slab output", {
  has_call <- function(sym) {
    dll <- getLoadedDLLs()[["neuroarchive"]]
    if (is.null(dll)) return(FALSE)
    routines <- getDLLRegisteredRoutines(dll)
    !is.null(routines) && sym %in% names(routines$`.Call`)
  }

  skip_if_not(has_call("_neuroarchive_quantize_voxel_slab_rcpp"))
  skip_if_not(has_call("_neuroarchive_quantize_voxel_slab_fill_rcpp"))

  set.seed(2)
  arr <- array(rnorm(8 * 7 * 6 * 9), dim = c(8, 7, 6, 9))
  slab_dims <- c(3L, 2L, 2L)
  nt <- dim(arr)[4]

  q_buf <- array(integer(prod(slab_dims) * nt), dim = c(slab_dims, nt))
  scale_buf <- array(numeric(prod(slab_dims)), dim = slab_dims)
  offset_buf <- array(numeric(prod(slab_dims)), dim = slab_dims)

  bits <- 12L
  method <- "sd"
  center <- FALSE
  starts <- list(c(1L, 1L, 1L), c(4L, 3L, 5L))

  for (start in starts) {
    ref <- neuroarchive:::quantize_voxel_slab_rcpp(arr, start, slab_dims, bits, method, center)
    nclip <- neuroarchive:::quantize_voxel_slab_fill_rcpp(
      arr, start, slab_dims, bits, method, center,
      q_buf, scale_buf, offset_buf
    )

    expect_equal(q_buf, ref$q)
    expect_equal(scale_buf, ref$scale, tolerance = 1e-12)
    expect_equal(offset_buf, ref$offset, tolerance = 1e-12)
    expect_equal(nclip, ref$n_clipped_total)
  }
})
