library(neuroarchive)
library(withr)

FakeSpace <- function(dim, spacing_v, origin_v=c(0,0,0)) {
  structure(list(dim=dim, spacing=spacing_v, origin=origin_v), class="FakeSpace")
}
space.FakeLogicalNeuroVol <- function(x, ...) attr(x,"space")
spacing.FakeSpace <- function(x, ...) x$spacing
origin.FakeSpace <- function(x, ...) x$origin
as.array.FakeLogicalNeuroVol <- function(x, ...) x$arr


test_that("basis.empirical_hrbf_compressed roundtrip", {
  mask <- array(TRUE, dim=c(1,1,2))
  vol <- structure(list(arr=mask), class="LogicalNeuroVol")
  attr(vol,"space") <- FakeSpace(c(1,1,2), c(1,1,1))

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.FakeLogicalNeuroVol", space.FakeLogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("origin.FakeSpace", origin.FakeSpace, envir=.GlobalEnv)
  assign("as.array.FakeLogicalNeuroVol", as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  defer({
    rm(FakeSpace, space.FakeLogicalNeuroVol, spacing.FakeSpace,
       origin.FakeSpace, as.array.FakeLogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  B <- diag(2)
  plan <- Plan$new()
  tmp <- local_tempfile(fileext=".h5")
  h5 <- H5File$new(tmp, mode="w")
  h5$create_group("/transforms")

  dict_desc <- list(type="spat.hrbf", params=list(sigma0=1, levels=0, radius_factor=2.5,
                                                  kernel_type="gaussian", seed=1))
  neuroarchive:::write_json_descriptor(h5[["/" ]], "dict.json", dict_desc)

  desc <- list(type="basis.empirical_hrbf_compressed",
               params=list(hrbf_dictionary_descriptor_path="/dict.json",
                           svd_rank=2L, omp_quant_bits=8L),
               inputs=c("input"))

  h <- DataHandle$new(initial_stash=list(input=B), plan=plan, h5=h5,
                      mask_info=list(mask=vol, active_voxels=2))
  h2 <- neuroarchive:::forward_step.basis.empirical_hrbf_compressed("basis.empirical_hrbf_compressed", desc, h)

  vt_path <- "/basis/00_basis.empirical_hrbf_compressed/vt_matrix"
  neuroarchive:::h5_write_dataset(h5[["/"]], vt_path, h2$stash$hrbf_vt)

  inv_desc <- h2$plan$descriptors[[1]]
  handle_inv <- DataHandle$new(initial_stash=list(hrbf_codes=h2$stash$hrbf_codes),
                               h5=h5, mask_info=list(mask=vol, active_voxels=2))
  h3 <- neuroarchive:::invert_step.basis.empirical_hrbf_compressed("basis.empirical_hrbf_compressed", inv_desc, handle_inv)

  expect_equal(dim(h3$stash$basis_matrix), dim(B))
  expect_lt(mean(abs(h3$stash$basis_matrix - B)), 1)
})
