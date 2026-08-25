library(testthat)

qpp_public_exports <- c(
  "LNANeuroVec",
  "analysis_contract",
  "as_lna_neurovec",
  "as_neurovec",
  "assert_split_compatible",
  "lna_affine_identity",
  "lna_close",
  "lna_grid_identity",
  "lna_loaded",
  "lna_profile_analysis_safe",
  "lna_profile_multiresolution_run",
  "lna_support_identity"
)

namespace_exports <- function() {
  getNamespaceExports("neuroarchive")
}

test_that("QPP and analysis-safe public entry points are exported", {
  expect_setequal(
    intersect(namespace_exports(), qpp_public_exports),
    qpp_public_exports
  )
})

test_that("QPP codec and provenance internals stay private", {
  internal_only <- c(
    "qpp_encode_cpp",
    "qpp_decode_cpp",
    "forward_step.qpp",
    "invert_step.qpp",
    "forward_step.parcel.mean",
    "invert_step.parcel.mean",
    ".lna_grid_id",
    ".lna_support_id"
  )
  expect_length(intersect(namespace_exports(), internal_only), 0L)
})
