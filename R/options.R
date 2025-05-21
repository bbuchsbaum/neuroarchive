#' Package Options for LNA
#'
#' Minimal stub implementation storing options in a package environment.
#'
#' @param ... Named options to set, or names of options to retrieve.
#' @return A list of current options or requested values.
#' @keywords internal
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
assign(".lna_opts", new.env(parent = emptyenv()), envir = lna_options_env)
