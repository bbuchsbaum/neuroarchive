library(testthat)

expected_exports <- c(
  "LNAFacade",
  "LNANeuroVec",
  "TransformBuilder",
  "analysis_contract",
  "advance_lna_execution",
  "append_lna_journal_event",
  "apply_template",
  "as_dense_mat",
  "as_latent_neurovec",
  "as_lna_neurovec",
  "as_lna_source_record",
  "as_neurovec",
  "as_pipeline",
  "as_vol4d",
  "assert_split_compatible",
  "basis_space_centers",
  "basis_space_embed",
  "basis_space_reconstruct",
  "basis_time",
  "benchmark_sobel3d_sizes",
  "check_edge_adaptive_performance",
  "check_transform_implementation",
  "cleanup_transform",
  "close_lna_journal",
  "compress_fmri",
  "delta",
  "delta_transform",
  "dequantize_data",
  "embed",
  "fmri_workflow",
  "fit_lna_shared_bases",
  "forward_step",
  "get_openmp_threads",
  "get_transform_capabilities",
  "get_transform_info",
  "hrbf",
  "hrbf_generate_basis",
  "hrbf_project_matrix",
  "hrbf_reconstruct_matrix",
  "heartbeat_lna_journal",
  "invert_step",
  "is_nbit_dataset",
  "is_transform_registered",
  "list_h5_filters",
  "list_shared_bases",
  "list_transforms",
  "lna_close",
  "lna_affine_identity",
  "lna_dataset_declaration",
  "lna_dataset_view",
  "lna_default.aggregate_runs",
  "lna_default.basis",
  "lna_default.basis.empirical_hrbf_compressed",
  "lna_default.delta",
  "lna_default.embed",
  "lna_default.embed.shared_basis",
  "lna_default.quant",
  "lna_default.sparsepca",
  "lna_default.spat.haar_octwave",
  "lna_default.spat.hrbf",
  "lna_default.spat.hrbf_project",
  "lna_default.temporal",
  "lna_get_quant_report",
  "lna_get_transform_report",
  "lna_grid_identity",
  "lna_loaded",
  "lna_execution_state",
  "lna_error_contract",
  "lna_file_checksum",
  "lna_manifest",
  "lna_manifest_from_journal",
  "lna_metadata_policy",
  "lna_options",
  "lna_profile_analysis_safe",
  "lna_profile_multiresolution_run",
  "lna_reference_assertion",
  "lna_source_record",
  "lna_source_records_from_bidser",
  "lna_support_identity",
  "lna_shared_basis_groups",
  "lna_spatial_compatibility",
  "lna_template_registry_env",
  "lna_transfer_inventory",
  "lna_verb_registry_env",
  "lna_write",
  "neuroim2_space_to_lna_header",
  "open_lna",
  "open_lna_dataset",
  "open_lna_journal",
  "override_lna_reference",
  "pca",
  "plan_lna_dataset",
  "plan_lna_dataset_update",
  "prepare_transform",
  "quant",
  "quantize_data",
  "query_lna",
  "rcpp_control",
  "read_lna",
  "read_lna_dataset_declaration",
  "read_lna_dataset_plan",
  "read_lna_dataset_update",
  "read_lna_journal",
  "read_lna_manifest",
  "read_lna_restoration",
  "read_shared_basis",
  "resolve_basis",
  "restore_bids",
  "ready_lna_execution_nodes",
  "resume_lna_journal",
  "scaffold_transform",
  "simple_basis_transform",
  "simple_function_transform",
  "simple_matrix_transform",
  "sobel3d_magnitude_rcpp",
  "sparse_pca",
  "suggest_dpss_fmri",
  "temporal",
  "temporal_basis",
  "temporal_project_denoise",
  "test_sobel3d_performance",
  "transfer_lna_dataset",
  "transform_min_dims",
  "transform_summary",
  "validate_lna",
  "validate_lna_dataset_declaration",
  "validate_lna_dataset_plan",
  "validate_lna_execution",
  "validate_lna_error_contract",
  "validate_lna_manifest",
  "validate_lna_restoration",
  "validate_lna_metadata_closure",
  "validate_lna_metadata_policy",
  "validate_lna_source_records",
  "validate_lna_transfer_paths",
  "validate_params",
  "verify_lna_transfer",
  "write_lna",
  "write_lna_dataset_declaration",
  "write_lna_dataset_plan",
  "write_lna_dataset",
  "write_lna_manifest",
  "write_shared_basis",
  "update_lna_dataset"
)

package_root <- function() {
  path <- getwd()
  while (!file.exists(file.path(path, "NAMESPACE"))) {
    parent <- dirname(path)
    if (identical(parent, path)) {
      path <- system.file(package = "neuroarchive")
      break
    }
    path <- parent
  }
  if (!file.exists(file.path(path, "NAMESPACE"))) {
    stop("Could not find package NAMESPACE")
  }
  path
}

namespace_exports <- function(path = package_root()) {
  ns_lines <- readLines(file.path(path, "NAMESPACE"), warn = FALSE)
  exports <- grep("^export\\(", ns_lines, value = TRUE)
  sub("^export\\((.*)\\)$", "\\1", exports)
}

test_that("public API exports change only intentionally", {
  exports <- namespace_exports()
  expect_equal(sort(exports), sort(expected_exports))
})

test_that("unsafe unreleased dataset API is absent", {
  retired <- c(
    "lna_dataset", "init_lna_dataset", "find_lna_files", "read_subject",
    "write_subject", "subjects", "bases", "participants",
    "dataset_description"
  )
  exports <- namespace_exports()
  expect_length(intersect(exports, retired), 0L)

  namespace <- readLines(file.path(package_root(), "NAMESPACE"), warn = FALSE)
  expect_false(any(grepl(
    "^S3method\\([^,]+,lna_dataset\\)$", namespace
  )))
})

test_that("neuroarchive does not collide with bidser dataset metadata names", {
  fixture <- testthat::test_path("fixtures", "bidser-exports.txt")
  bidser_exports <- trimws(readLines(fixture, warn = FALSE))
  bidser_exports <- bidser_exports[
    nzchar(bidser_exports) & !startsWith(bidser_exports, "#")
  ]

  expect_setequal(bidser_exports, c("dataset_description", "participants"))
  expect_length(intersect(namespace_exports(), bidser_exports), 0L)
})

test_that("retired dataset documentation and Stable claims are absent", {
  root <- package_root()
  retired_rd <- file.path(
    root,
    "man",
    paste0(
      c(
        "lna_dataset", "init_lna_dataset", "find_lna_files", "read_subject",
        "write_subject", "subjects", "bases", "participants",
        "dataset_description"
      ),
      ".Rd"
    )
  )
  expect_false(any(file.exists(retired_rd)))

  rd_files <- list.files(file.path(root, "man"), pattern = "[.]Rd$", full.names = TRUE)
  rd_text <- unlist(lapply(rd_files, readLines, warn = FALSE), use.names = FALSE)
  expect_false(any(grepl(
    "Stable. This is .*BIDS-style LNA dataset",
    rd_text
  )))
})
