#' Robust Centering Transform (Plugin)
#'
#' Subtracts a per-column center (median or mean) and stores the centers.
#' Inverse adds centers back to reconstruct.
#'
#' @keywords internal
forward_step.robust_center <- function(type, desc, handle) {
  p <- desc$params %||% list()
  method <- p$method %||% "median"

  input_key <- if (!is.null(desc$inputs) && length(desc$inputs) > 0) {
    desc$inputs[[1]]
  } else {
    inp <- handle$pull_first(c("input", "dense_mat", "aggregated_matrix"))
    inp$key
  }

  X <- handle$get_inputs(input_key)[[1]]
  X <- neuroarchive::as_dense_mat(X)

  centers <- switch(method,
    median = apply(X, 2, median, na.rm = TRUE),
    mean   = colMeans(X, na.rm = TRUE),
    stop("Unknown method: ", method)
  )
  centered <- sweep(X, 2, centers, "-")

  # Descriptor + datasets
  plan <- handle$plan
  fname <- plan$get_next_filename("robust_center")
  base_name <- tools::file_path_sans_ext(fname)

  run_id <- handle$current_run_id %||% "run-01"
  run_id <- neuroarchive:::sanitize_run_id(run_id)
  scans_root <- neuroarchive::lna_options("paths.scans_root")[[1]]
  tr_name <- paste0("robust_center_", base_name)
  coeff_path <- paste0(scans_root, run_id, "/", tr_name, "/coefficients")
  centers_path <- paste0(scans_root, run_id, "/", tr_name, "/centers")

  p$method <- method
  desc$params <- p
  desc$version <- "1.0"
  desc$inputs <- c(input_key)
  desc$outputs <- c("robust_centered")
  desc$datasets <- list(
    list(path = coeff_path, role = "coefficients"),
    list(path = centers_path, role = "centers")
  )

  # Write payloads + defs
  plan$add_descriptor(fname, desc)
  plan$add_payload(coeff_path, centered)
  plan$add_dataset_def(coeff_path, "coefficients", "robust_center",
                       run_id, as.integer(plan$next_index - 1L),
                       jsonlite::toJSON(p, auto_unbox = TRUE),
                       coeff_path, "eager", dtype = "float32")
  plan$add_payload(centers_path, centers)
  plan$add_dataset_def(centers_path, "centers", "robust_center",
                       run_id, as.integer(plan$next_index - 1L),
                       jsonlite::toJSON(p, auto_unbox = TRUE),
                       centers_path, "eager", dtype = "float32")

  handle$plan <- plan
  handle$update_stash(keys = input_key,
                      new_values = setNames(list(centered), "robust_centered"))
}

#' @keywords internal
invert_step.robust_center <- function(type, desc, handle) {
  # Identify dataset paths
  ds <- desc$datasets
  get_path <- function(role) {
    idx <- which(vapply(ds, function(d) d$role, character(1)) == role)
    if (length(idx) > 0) ds[[idx[1]]]$path else NULL
  }
  coeff_path <- get_path("coefficients")
  centers_path <- get_path("centers")

  # Take coefficients from stash if present, else read
  coeff_key <- desc$outputs[[1]] %||% "robust_centered"
  coeff <- if (handle$has_key(coeff_key)) {
    handle$get_inputs(coeff_key)[[coeff_key]]
  } else {
    root <- handle$h5[["/"]]
    neuroarchive:::h5_read(root, coeff_path)
  }
  # Read centers
  root <- handle$h5[["/"]]
  centers <- neuroarchive:::h5_read(root, centers_path)

  dense <- sweep(coeff, 2, centers, "+")
  input_key <- desc$inputs[[1]] %||% "input"
  handle$update_stash(keys = coeff_key,
                      new_values = setNames(list(dense), input_key))
}

#' @export
lna_default.robust_center <- function() {
  list(method = "median")
}
