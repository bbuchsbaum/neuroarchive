library(testthat)

test_that("transform files do not call low-level math internals directly", {
  pkg_root <- normalizePath(
    file.path(testthat::test_path(), "..", ".."),
    winslash = "/",
    mustWork = TRUE
  )
  transform_files <- sort(Sys.glob(file.path(pkg_root, "R", "transform_*.R")))
  if (length(transform_files) == 0) {
    skip("source transform files are not available in the installed test environment")
  }
  expect_gt(length(transform_files), 0)

  # These functions represent low-level math internals. Transform code should
  # route through backend-aware wrappers instead.
  forbidden_symbols <- c(
    "hrbf_basis_from_params_internal",
    "poisson_disk_sample_neuroim2",
    "generate_hrbf_atom",
    "hrbf_atoms_rcpp",
    "forward_lift_R",
    "forward_lift_rcpp",
    "forward_lift_matrix_rcpp",
    "inverse_lift_R",
    "inverse_lift_rcpp",
    "inverse_lift_matrix_rcpp"
  )

  legacy_exemptions <- list()

  findings <- character()
  for (file in transform_files) {
    lines <- readLines(file, warn = FALSE)
    file_exemptions <- legacy_exemptions[[file]]
    if (is.null(file_exemptions)) {
      file_exemptions <- character()
    }
    for (sym in forbidden_symbols) {
      if (sym %in% file_exemptions) {
        next
      }
      pattern <- paste0("\\b", sym, "\\s*\\(")
      hit_lines <- which(grepl(pattern, lines, perl = TRUE))
      if (length(hit_lines) > 0) {
        rel_file <- sub(paste0("^", pkg_root, "/"), "", file)
        findings <- c(
          findings,
          sprintf("%s:%s -> %s", rel_file, paste(hit_lines, collapse = ","), sym)
        )
      }
    }
  }

  if (length(findings) > 0) {
    fail(
      paste(
        c(
          "Direct low-level math calls found in transform files:",
          findings
        ),
        collapse = "\n"
      )
    )
  } else {
    succeed()
  }
})
