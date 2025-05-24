#' Convert a neuroim2 NeuroSpace to an LNA header list
#'
#' This helper converts a `neuroim2` `NeuroSpace` object to the
#' named list format expected for the `header` argument of
#' [write_lna()].
#'
#' @param neurospace_obj A `NeuroSpace` object from the `neuroim2` package.
#' @return A named list with `dims`, `spacing`, `origin`, and `transform` fields.
#' @export
neuroim2_space_to_lna_header <- function(neurospace_obj) {
  if (missing(neurospace_obj)) {
    abort_lna("neurospace_obj is required", .subclass = "lna_error_validation",
              location = "neuroim2_space_to_lna_header")
  }

  # Helper function to try global override first, then neuroim2 namespace
  safe_call <- function(fn_name, obj) {
    if (exists(fn_name, envir = .GlobalEnv, mode = "function")) {
      get(fn_name, envir = .GlobalEnv)(obj)
    } else {
      get(fn_name, envir = asNamespace("neuroim2"))(obj)
    }
  }

  dims <- dim(neurospace_obj)
  nd <- length(dims)
  header_list <- list(
    dims = dims[seq_len(min(3L, nd))],
    spacing = safe_call("spacing", neurospace_obj),
    origin = safe_call("origin", neurospace_obj),
    transform = safe_call("trans", neurospace_obj)
  )

  header_list
}

