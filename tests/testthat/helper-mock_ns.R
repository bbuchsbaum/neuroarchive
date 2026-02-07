ns_name <- "neuroarchive"
ns_env <- asNamespace(ns_name)
pkg_env <- as.environment(paste0("package:", ns_name))

# Do NOT cache original_bindings at file load time — devtools::load_all()
# may replace them later and restoring stale copies causes "unused argument"
# errors (e.g. core_write missing the h5 parameter).

reset_neuroarchive_bindings <- function(names = NULL, target_env = NULL) {
  # No-op: callers should use local_mock_ns with withr::defer instead.
  invisible(NULL)
}

local_mock_ns <- function(bindings, ns = ns_name, envir = parent.frame()) {

  stopifnot(is.list(bindings))
  ns_env_local <- if (identical(ns, ns_name)) ns_env else asNamespace(ns)
  pkg_env_local <- if (identical(ns, ns_name)) pkg_env else as.environment(paste0("package:", ns))

  # Capture originals at call time, not at file load time
  originals <- lapply(names(bindings), function(name) get(name, envir = ns_env_local))
  names(originals) <- names(bindings)

  # Manually replace bindings instead of using local_mocked_bindings
  # to avoid double-cleanup conflicts with S3 generics.
  for (name in names(bindings)) {
    if (bindingIsLocked(name, ns_env_local)) unlockBinding(name, ns_env_local)
    assign(name, bindings[[name]], envir = ns_env_local)
    lockBinding(name, ns_env_local)
    if (exists(name, envir = pkg_env_local, inherits = FALSE)) {
      assign(name, bindings[[name]], envir = pkg_env_local)
    }
  }

  withr::defer({
    for (name in names(originals)) {
      if (bindingIsLocked(name, ns_env_local)) unlockBinding(name, ns_env_local)
      assign(name, originals[[name]], envir = ns_env_local)
      lockBinding(name, ns_env_local)
      if (exists(name, envir = pkg_env_local, inherits = FALSE)) {
        assign(name, originals[[name]], envir = pkg_env_local)
      }
    }
  }, envir = envir)
}

local_tempfile <- withr::local_tempfile
