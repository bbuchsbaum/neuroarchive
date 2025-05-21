#' Package Options for LNA
#'
#' Provides a lightweight mechanism for storing global package defaults.
#' Options are kept in an internal environment and can be retrieved or
#' updated via this helper.  Typical options include
#' `write.compression_level`, `write.chunk_target_mib` and per-transform
#' defaults such as `quant` or `delta` lists.
#'
#' @param ... Named options to set, or character names of options to
#'   retrieve.  If no arguments are provided, the full option list is
#'   returned.
#' @return A list of current options or the requested subset.  When setting
#'   values the updated option list is returned invisibly.
#' @export
lna_options <- function(...) {
  .lna_opts <- get(".lna_opts", envir = lna_options_env)
  args <- list(...)
  if (length(args) == 0) {
    return(as.list(.lna_opts))
  }
  if (is.null(names(args))) {
    return(mget(unlist(args), envir = .lna_opts, ifnotfound = list(NULL)))
  }
  for (nm in names(args)) {
    assign(nm, args[[nm]], envir = .lna_opts)
  }
  invisible(as.list(.lna_opts))
}

lna_options_env <- new.env(parent = emptyenv())
default_opts <- list(
  write.compression_level = 0L,
  write.chunk_target_mib = 1,
  quant = list(),
  delta = list()
)
assign(".lna_opts", list2env(default_opts, parent = emptyenv()),
       envir = lna_options_env)
