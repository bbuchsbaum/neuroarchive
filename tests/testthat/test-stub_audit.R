library(testthat)

test_that("no exported function aborts with 'not yet implemented'", {
  # Get all exported functions from neuroarchive
  ns <- asNamespace("neuroarchive")
  exports <- getNamespaceExports(ns)

  # Check for functions that contain "not yet implemented" in their body
  stub_fns <- character()
  for (fn_name in exports) {
    fn <- tryCatch(get(fn_name, envir = ns), error = function(e) NULL)
    if (!is.function(fn)) next
    body_text <- tryCatch(paste(deparse(body(fn)), collapse = "\n"), error = function(e) "")
    if (grepl("not yet implemented", body_text, ignore.case = TRUE)) {
      stub_fns <- c(stub_fns, fn_name)
    }
  }

  if (length(stub_fns) > 0) {
    fail(sprintf("Functions with 'not yet implemented': %s",
                 paste(stub_fns, collapse = ", ")))
  } else {
    succeed()
  }
})
