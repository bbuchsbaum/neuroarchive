library(neuroarchive)
library(withr)

FakeSpace <- function(dim, spacing_v, origin_v=c(0,0,0)) {
  structure(list(dim=dim, spacing=spacing_v, origin=origin_v), class="FakeSpace")
}
space.LogicalNeuroVol <- function(x, ...) attr(x, "space")
spacing.FakeSpace <- function(x, ...) x$spacing
origin.FakeSpace <- function(x, ...) x$origin
as.array.LogicalNeuroVol <- function(x, ...) x$arr


test_that("invert_step.basis.empirical_hrbf_compressed returns basis", {
  mask <- array(TRUE, dim=c(1,1,1))
  vol <- structure(list(arr=mask), class="LogicalNeuroVol")
  attr(vol, "space") <- FakeSpace(c(1,1,1), c(1,1,1))

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.LogicalNeuroVol", space.LogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("origin.FakeSpace", origin.FakeSpace, envir=.GlobalEnv)
  assign("as.array.LogicalNeuroVol", as.array.LogicalNeuroVol, envir=.GlobalEnv)
  defer({
    rm(FakeSpace, space.LogicalNeuroVol, spacing.FakeSpace,
       origin.FakeSpace, as.array.LogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  tmp <- tempfile(fileext=".h5")
  on.exit(unlink(tmp), add = TRUE)
  h5 <- H5File$new(tmp, mode="w")
  vt_mat <- matrix(1, nrow=1, ncol=1)
  neuroarchive:::h5_write_dataset(h5[["/"]], "/basis/vt", vt_mat)
  dict_desc <- list(type="spat.hrbf", params=list(sigma0=1, levels=0, radius_factor=2.5, kernel_type="gaussian", seed=1))
  neuroarchive:::write_json_descriptor(h5[["/"]], "dict.json", dict_desc)

  desc <- list(
    type = "basis.empirical_hrbf_compressed",
    params = list(hrbf_dictionary_descriptor_path="/dict.json"),
    datasets = list(list(path="/basis/vt", role="svd_vt")),
    inputs = c("basis_matrix"),
    outputs = c("hrbf_codes")
  )

  handle <- DataHandle$new(initial_stash=list(hrbf_codes=matrix(1L, nrow=1, ncol=1)),
                           h5=h5, mask_info=list(mask=vol, active_voxels=1))
  h2 <- neuroarchive:::invert_step.basis.empirical_hrbf_compressed("basis.empirical_hrbf_compressed", desc, handle)
  expect_true(h2$has_key("basis_matrix"))
  expect_false(h2$has_key("hrbf_codes"))
  expect_equal(dim(h2$stash$basis_matrix), c(1, length(as.array(vol))))
  h5$close_all()
})
