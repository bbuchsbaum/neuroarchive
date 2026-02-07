local_mock_write_lna <- function(mock_fn, envir = parent.frame()) {
  ns <- asNamespace("neuroarchive")
  pkg_env <- as.environment("package:neuroarchive")
  original_ns <- get("write_lna", envir = ns)
  original_pkg <- get("write_lna", envir = pkg_env)

  # Manually replace the binding — avoid local_mocked_bindings which
  # registers its own withr::defer that conflicts with S3 generics.
  if (bindingIsLocked("write_lna", ns)) unlockBinding("write_lna", ns)
  assign("write_lna", mock_fn, envir = ns)
  lockBinding("write_lna", ns)

  if (bindingIsLocked("write_lna", pkg_env)) unlockBinding("write_lna", pkg_env)
  assign("write_lna", mock_fn, envir = pkg_env)
  lockBinding("write_lna", pkg_env)

  withr::defer({
    if (bindingIsLocked("write_lna", ns)) unlockBinding("write_lna", ns)
    assign("write_lna", original_ns, envir = ns)
    lockBinding("write_lna", ns)

    if (bindingIsLocked("write_lna", pkg_env)) unlockBinding("write_lna", pkg_env)
    assign("write_lna", original_pkg, envir = pkg_env)
    lockBinding("write_lna", pkg_env)
  }, envir = envir)
}
