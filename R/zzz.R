#' Package Startup Functions
#'
#' @description
#' Functions that run when the package is loaded or attached.
#'
#' @keywords internal

.onLoad <- function(libname, pkgname) {
  # Ensure our S3 shims for hdf5r are registered even when the package is only
  # loaded (not attached). This avoids brittle path handling differences across
  # hdf5r versions in tests and internal helpers.
  if (requireNamespace("hdf5r", quietly = TRUE)) {
    try(registerS3method("[[<-", "H5File", `[[<-.H5File`, envir = asNamespace(pkgname)),
        silent = TRUE)
  }

  # Auto-discover and register transforms
  tryCatch({
    n_transforms <- discover_and_register_transforms(validate = FALSE)
    
    # Log if in verbose mode
    if (isTRUE(getOption("neuroarchive.verbose"))) {
      packageStartupMessage(sprintf(
        "neuroarchive: Registered %d transforms", 
        n_transforms
      ))
    }
  }, error = function(e) {
    # Silently fail during package load to avoid breaking package
    # Users can manually call discover_and_register_transforms() if needed
    if (isTRUE(getOption("neuroarchive.debug"))) {
      warning(sprintf(
        "Failed to auto-register transforms: %s", 
        conditionMessage(e)
      ), call. = FALSE)
    }
  })
  
  invisible()
}

.onAttach <- function(libname, pkgname) {
  # Display package startup message
  packageStartupMessage(sprintf(
    "neuroarchive %s - Latent NeuroArchive (LNA) tools",
    utils::packageVersion(pkgname)
  ))
  
  # Show registry summary if verbose
  if (isTRUE(getOption("neuroarchive.verbose"))) {
    transforms <- list_transforms()
    if (length(transforms) > 0) {
      packageStartupMessage(sprintf(
        "Available transforms: %s",
        paste(transforms, collapse = ", ")
      ))
    }
  }
  
  invisible()
}

.onUnload <- function(libpath) {
  # Clear transform registry on unload
  clear_transform_registry()
  
  invisible()
}
