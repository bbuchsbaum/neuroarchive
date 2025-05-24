library(neuroarchive)
library(withr)

FakeSpace <- function(dim, spacing_v, origin_v=c(0,0,0)) {
  structure(list(dim=dim, spacing=spacing_v, origin=origin_v), class="FakeSpace")
}
space.LogicalNeuroVol <- function(x, ...) attr(x,"space")
spacing.FakeSpace <- function(x, ...) x$spacing
origin.FakeSpace <- function(x, ...) x$origin
as.array.LogicalNeuroVol <- function(x, ...) x$arr


test_that("basis.empirical_hrbf_compressed roundtrip", {
  mask <- array(TRUE, dim=c(1,1,1))
  vol <- structure(list(arr=mask), class="LogicalNeuroVol")
  attr(vol,"space") <- FakeSpace(c(1,1,1), c(1,1,1))

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.LogicalNeuroVol", space.LogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("origin.FakeSpace", origin.FakeSpace, envir=.GlobalEnv)
  assign("as.array.LogicalNeuroVol", as.array.LogicalNeuroVol, envir=.GlobalEnv)
  defer({
    rm(FakeSpace, space.LogicalNeuroVol, spacing.FakeSpace,
       origin.FakeSpace, as.array.LogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  B <- matrix(1, nrow=1, ncol=1)  # Simplified to 1x1
  plan <- Plan$new()
  tmp <- tempfile(fileext=".h5")
  on.exit(unlink(tmp), add = TRUE)
  h5 <- H5File$new(tmp, mode="w")
  h5$create_group("/transforms")

  dict_desc <- list(type="spat.hrbf", params=list(sigma0=1, levels=0, radius_factor=10,
                                                  kernel_type="gaussian", seed=1))
  neuroarchive:::write_json_descriptor(h5[["/" ]], "dict.json", dict_desc)

  desc <- list(type="basis.empirical_hrbf_compressed",
               params=list(hrbf_dictionary_descriptor_path="/dict.json",
                           svd_rank=1L, omp_quant_bits=8L, omp_tol=0.1, omp_sparsity_limit=1L),
               inputs=c("input"))

  h <- DataHandle$new(initial_stash=list(input=B), plan=plan, h5=h5,
                      mask_info=list(mask=vol, active_voxels=1))
  h2 <- neuroarchive:::forward_step.basis.empirical_hrbf_compressed("basis.empirical_hrbf_compressed", desc, h)

  vt_path <- "/basis/00_basis.empirical_hrbf_compressed/vt_matrix"
  neuroarchive:::h5_write_dataset(h5[["/"]], vt_path, h2$stash$hrbf_vt)

  inv_desc <- h2$plan$descriptors[[1]]
  handle_inv <- DataHandle$new(initial_stash=list(hrbf_codes=h2$stash$hrbf_codes),
                               h5=h5, mask_info=list(mask=vol, active_voxels=1))
  h3 <- neuroarchive:::invert_step.basis.empirical_hrbf_compressed("basis.empirical_hrbf_compressed", inv_desc, handle_inv)

  expect_true(h3$has_key("input"))
  result <- h3$stash$input
  
  # Convert Matrix objects to regular matrix
  if (inherits(result, "Matrix")) {
    result <- as.matrix(result)
  }
  
  expect_equal(dim(result), dim(B))
  expect_true(is.numeric(result))
  expect_true(all(is.finite(result)))
})
