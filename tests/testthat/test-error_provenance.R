library(testthat)

# forward_step error provenance

forward_step.fail <- function(type, desc, handle) {
  stop("boom")
}
assign("forward_step.fail", forward_step.fail, envir = .GlobalEnv)
withr::defer(rm("forward_step.fail", envir = .GlobalEnv))

test_that("core_write reports step provenance", {
  # Mock default_params to prevent schema not found warning for "fail" type
  original_default_params <- if (exists("default_params", envir = asNamespace("neuroarchive"))) {
    get("default_params", envir = asNamespace("neuroarchive"))
  } else { NA }
  
  mocked_dp <- function(type) {
    if (type == "fail") return(list())
    if (is.function(original_default_params)) return(original_default_params(type))
    return(list()) # Fallback if original couldn't be found/isn't function
  }
  
  # Use local_mocked_bindings if default_params is an exported or known binding.
  # If default_params is not exported and local_mocked_bindings fails, 
  # we might need a more direct unlockBinding/assignInNamespace approach for the mock, 
  # though local_mocked_bindings is preferred.
  # For now, assuming default_params can be shimmed by local_mocked_bindings.
  local_mocked_bindings(
    default_params = mocked_dp,
    .env = asNamespace("neuroarchive")
  )

  err_cw_prov <- expect_error(core_write(x = array(1, dim = c(1, 1, 1)), transforms = "fail"))
  # cat("\n--- Diagnostic for core_write reports step provenance ---\n")
  # cat("Error Class: ", paste(class(err_cw_prov), collapse=", "), "\n")
  # cat("Error Message: ", conditionMessage(err_cw_prov), "\n")
  # cat("Error Location: ", err_cw_prov$location, "\n")
  # cat("--- End Diagnostic ---\n")
  expect_true(grepl("forward_step.fail[0]", err_cw_prov$location, fixed = TRUE))
})

# invert_step error provenance
invert_step.fail <- function(type, desc, handle) {
  stop("oops")
}
assign("invert_step.fail", invert_step.fail, envir = .GlobalEnv)
withr::defer(rm("invert_step.fail", envir = .GlobalEnv))

create_fail_file <- function(path) {
  h5 <- neuroarchive:::open_h5(path, mode = "w")
  tf <- h5$create_group("transforms")
  # Ensure a run exists for core_read to attempt processing transforms
  scans_group <- h5$create_group("scans")
  scans_group$create_group("run-01") 
  write_json_descriptor(tf, "00_fail.json", list(type = "fail"))
  neuroarchive:::close_h5_safely(h5)
}

test_that("core_read reports step provenance", {
  tmp <- local_tempfile(fileext = ".h5")
  create_fail_file(tmp)
  err_cr_prov <- expect_error(core_read(tmp))
  # cat("\n--- Diagnostic for core_read reports step provenance ---\n")
  # cat("Error Class: ", paste(class(err_cr_prov), collapse=", "), "\n")
  # cat("Error Message: ", conditionMessage(err_cr_prov), "\n")
  # cat("Error Location: ", err_cr_prov$location, "\n")
  # cat("--- End Diagnostic ---\n")
  expect_true(grepl("invert_step.fail[0]", err_cr_prov$location, fixed = TRUE))
})

create_double_fail_file <- function(path) {
  h5 <- neuroarchive:::open_h5(path, mode = "w")
  tf <- h5$create_group("transforms")
  # Ensure a run exists for core_read to attempt processing transforms
  scans_group <- h5$create_group("scans")
  scans_group$create_group("run-01") 
  write_json_descriptor(tf, "00_dummy.json", list(type = "dummy"))
  write_json_descriptor(tf, "01_fail.json", list(type = "fail"))
  neuroarchive:::close_h5_safely(h5)
}

test_that("core_read error location uses transform index", {
  tmp <- local_tempfile(fileext = ".h5")
  create_double_fail_file(tmp)
  invert_step.dummy <- function(type, desc, handle) handle
  assign("invert_step.dummy", invert_step.dummy, envir = .GlobalEnv)
  withr::defer(rm("invert_step.dummy", envir = .GlobalEnv))

  err_cr_idx <- expect_error(core_read(tmp))
  # cat("\n--- Diagnostic for core_read error location uses transform index ---\n")
  # cat("Error Class: ", paste(class(err_cr_idx), collapse=", "), "\n")
  # cat("Error Message: ", conditionMessage(err_cr_idx), "\n")
  # cat("Error Location: ", err_cr_idx$location, "\n")
  # cat("--- End Diagnostic ---\n")
  expect_true(grepl("invert_step.fail[1]", err_cr_idx$location, fixed = TRUE))
})
