library(neuroarchive)
library(withr)

FakeSpace <- function(dim, spacing_v, origin_v=c(0,0,0)) {
  structure(list(dim=dim, spacing=spacing_v, origin=origin_v), class="FakeSpace")
}
space.FakeLogicalNeuroVol <- function(x, ...) attr(x,"space")
spacing.FakeSpace <- function(x, ...) x$spacing
origin.FakeSpace <- function(x, ...) x$origin
as.array.LogicalNeuroVol <- function(x, ...) {
  arr <- x$arr
  if (is.null(dim(arr))) {
    # If arr lost its dimensions, restore them from the space
    spc <- attr(x, "space")
    if (!is.null(spc) && !is.null(spc$dim)) {
      dim(arr) <- spc$dim
    }
  }
  arr
}


test_that("embed.transfer_hrbf_basis computes coefficients and inverse", {
  mask <- array(TRUE, dim=c(2,2,2))
  vol <- structure(list(arr=mask), class="LogicalNeuroVol")
  attr(vol,"space") <- FakeSpace(c(2,2,2), c(1,1,1))

  assign("FakeSpace", FakeSpace, envir=.GlobalEnv)
  assign("space.FakeLogicalNeuroVol", space.FakeLogicalNeuroVol, envir=.GlobalEnv)
  assign("spacing.FakeSpace", spacing.FakeSpace, envir=.GlobalEnv)
  assign("origin.FakeSpace", origin.FakeSpace, envir=.GlobalEnv)
  assign("as.array.LogicalNeuroVol", as.array.LogicalNeuroVol, envir=.GlobalEnv)
  defer({
    rm(FakeSpace, space.FakeLogicalNeuroVol, spacing.FakeSpace,
       origin.FakeSpace, as.array.LogicalNeuroVol, envir=.GlobalEnv)
  }, envir = parent.frame())

  B <- diag(8)
  plan <- Plan$new()
  tmp <- local_tempfile(fileext=".h5")
  h5 <- H5File$new(tmp, mode="w")
  h5$create_group("/transforms")

  dict_desc <- list(type="spat.hrbf", params=list(sigma0=0.2, levels=0, radius_factor=2.5,
                                                  kernel_type="gaussian", seed=1))
  neuroarchive:::write_json_descriptor(h5[["/" ]], "dict.json", dict_desc)

  # create compressed basis
  desc_basis <- list(type="basis.empirical_hrbf_compressed",
                     params=list(hrbf_dictionary_descriptor_path="/dict.json",
                                 svd_rank=8L, omp_quant_bits=8L),
                     inputs=c("input"))
  h_write <- DataHandle$new(initial_stash=list(input=B), plan=plan, h5=h5,
                            mask_info=list(mask=vol, active_voxels=8))
  h_write2 <- neuroarchive:::forward_step.basis.empirical_hrbf_compressed("basis.empirical_hrbf_compressed", desc_basis, h_write)
  vt_path <- "/basis/00_basis.empirical_hrbf_compressed/vt_matrix"
  codes_path <- "/basis/00_basis.empirical_hrbf_compressed/hrbf_codes"
  neuroarchive:::h5_write_dataset(h5[["/"]], vt_path, h_write2$stash$hrbf_vt)
  # convert codes to matrix for storage
  B_dict <- neuroarchive:::hrbf_basis_from_params(dict_desc$params, vol)
  k_dict <- nrow(B_dict)
  codes_mat <- matrix(0, nrow=length(h_write2$stash$hrbf_codes), ncol=k_dict)
  for (j in seq_along(h_write2$stash$hrbf_codes)) {
    cinfo <- h_write2$stash$hrbf_codes[[j]]
    if (length(cinfo$indices) > 0) {
      codes_mat[j, cinfo$indices] <- cinfo$q * cinfo$scale
    }
  }
  neuroarchive:::h5_write_dataset(h5[["/"]], codes_path, codes_mat)

  basis_desc_name <- names(h_write2$plan$descriptors)[1]
  neuroarchive:::write_json_descriptor(h5[["/transforms"]], basis_desc_name, h_write2$plan$descriptors[[1]])
  h5$close_all()

  X <- matrix(1:16, nrow=2)
  plan2 <- Plan$new()
  handle <- DataHandle$new(initial_stash=list(input_dense_mat=X), plan=plan2,
                           mask_info=list(mask=vol, active_voxels=8))
  desc_transfer <- list(type="embed.transfer_hrbf_basis",
                        params=list(source_lna_file_path=tmp,
                                    source_transform_descriptor_name=basis_desc_name))
  h_tr <- neuroarchive:::forward_step.embed.transfer_hrbf_basis("embed.transfer_hrbf_basis", desc_transfer, handle)

  B_emp <- neuroarchive:::.load_empirical_hrbf_basis(tmp, basis_desc_name, vol)
  expected_coef <- tcrossprod(X, B_emp)
  expect_true(h_tr$has_key("coefficients"))
  expect_equal(dim(h_tr$stash$coefficients), dim(expected_coef))

  desc_inv <- h_tr$plan$descriptors[[1]]
  handle_inv <- DataHandle$new(initial_stash=list(coefficients=h_tr$stash$coefficients),
                               mask_info=list(mask=vol, active_voxels=8))
  out <- neuroarchive:::invert_step.embed.transfer_hrbf_basis("embed.transfer_hrbf_basis", desc_inv, handle_inv)
  expect_equal(out$stash$input_dense_mat, X, tolerance = 1e-3)
})
