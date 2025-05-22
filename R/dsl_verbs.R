#' Initiate an LNA pipeline
#'
#' Creates a new `lna_pipeline` object and sets its input using
#' `lna_pipeline$set_input()`.
#'
#' @param x Data object or list of run data.
#' @param run_ids Optional character vector of run identifiers.
#' @param chunk_mb_suggestion Optional numeric hint for chunk size.
#'
#' @return A configured `lna_pipeline` object.
#' @export
as_pipeline <- function(x, run_ids = NULL, chunk_mb_suggestion = NULL) {
  pipe <- lna_pipeline$new()
  pipe$set_input(x, run_ids = run_ids, chunk_mb_suggestion = chunk_mb_suggestion)
  pipe
}

#' Execute an LNA pipeline
#'
#' Translates an `lna_pipeline` object into a call to `write_lna()` and
#' materialises the resulting archive on disk.
#'
#' @param pipeline_obj An `lna_pipeline` object.
#' @param file Path to the output `.h5` file.
#' @param ... Additional arguments forwarded to `write_lna()` such as
#'   `header`, `mask`, or `plugins`.
#' @param .verbose Logical flag controlling verbosity (currently unused).
#' @param .checksum Checksum mode forwarded to `write_lna()`.
#'
#' @return The result returned by `write_lna()`.
#' @export
lna_write <- function(pipeline_obj, file, ...,
                      .verbose = TRUE, .checksum = "sha256") {
  if (!inherits(pipeline_obj, "lna_pipeline")) {
    abort_lna(
      "pipeline_obj must be an lna_pipeline",
      .subclass = "lna_error_validation",
      location = "lna_write:pipeline_obj"
    )
  }

  transform_types <- vapply(pipeline_obj$steps, function(s) s$type, character(1))
  transform_params_list <- lapply(pipeline_obj$steps, function(s) s$params)
  names(transform_params_list) <- transform_types

  extra_args <- utils::modifyList(pipeline_obj$engine_opts %||% list(), list(...))

  args <- c(
    list(
      x = pipeline_obj$input,
      file = file,
      run_id = pipeline_obj$runs,
      transforms = transform_types,
      transform_params = transform_params_list
    ),
    extra_args
  )

  args$checksum <- .checksum

  result <- tryCatch(
    {
      do.call(write_lna, args)
    },
    lna_error = function(e) {
      step_idx_core <- attr(e, "step_index", exact = TRUE)
      ttype_core <- attr(e, "transform_type", exact = TRUE)

      bullet <- "Pipeline execution failed."
      if (!is.null(step_idx_core) && !is.null(ttype_core)) {
        bullet <- sprintf(
          "Pipeline failure in step %d (type='%s')",
          step_idx_core + 1, ttype_core
        )
      }

      rlang::abort(
        message = c(bullet, "i" = conditionMessage(e)),
        parent = e,
        .subclass = class(e)
      )
    },
    error = function(e) {
      rlang::abort(
        message = c("Pipeline execution failed.", "i" = conditionMessage(e)),
        parent = e
      )
    }
  )
  invisible(result)
}

##' Quantization DSL verb
#'
#' Adds a quantization step to a pipeline. If `data_or_pipe`
#' is not an `lna_pipeline`, a new pipeline is created via
#' `as_pipeline()`.
#'
#' Parameter values are resolved by merging schema defaults,
#' global `lna_options("quant")`, and user-supplied arguments.
#'
#' @param data_or_pipe Data object or `lna_pipeline`.
#' @param bits Optional number of quantization bits.
#' @param ... Additional parameters for the quant transform.
#'
#' @return An `lna_pipeline` object with the quant step appended.
#' @export
quant <- function(data_or_pipe, bits = NULL, ...) {
  pipe <- if (inherits(data_or_pipe, "lna_pipeline")) {
    data_or_pipe
  } else {
    as_pipeline(data_or_pipe)
  }

  user_params <- c(list(bits = bits), list(...))
  user_params <- user_params[!vapply(user_params, is.null, logical(1))]

  pars <- default_params("quant")
  opts <- lna_options("quant")$quant %||% list()
  pars <- utils::modifyList(pars, opts)
  pars <- utils::modifyList(pars, user_params)

  step_spec <- list(type = "quant", params = pars)
  pipe$add_step(step_spec)
  pipe
}

##' PCA DSL verb
#'
#' Adds a PCA basis computation step to a pipeline. If
#' `data_or_pipe` is not an `lna_pipeline`, a new pipeline is
#' created via `as_pipeline()`.
#'
#' Parameter values are resolved by merging schema defaults,
#' global `lna_options("basis.pca")`, and user-supplied arguments.
#'
#' @param data_or_pipe Data object or `lna_pipeline`.
#' @param k Optional number of components to compute.
#' @param ... Additional parameters for the PCA basis transform.
#'
#' @return An `lna_pipeline` object with the PCA step appended.
#' @export
pca <- function(data_or_pipe, k = NULL, ...) {
  pipe <- if (inherits(data_or_pipe, "lna_pipeline")) {
    data_or_pipe
  } else {
    as_pipeline(data_or_pipe)
  }

  user_params <- c(list(k = k), list(...))
  user_params <- user_params[!vapply(user_params, is.null, logical(1))]

  pars <- default_params("basis.pca")
  opts <- lna_options("basis.pca")$`basis.pca` %||% list()
  pars <- utils::modifyList(pars, opts)
  pars <- utils::modifyList(pars, user_params)

  step_spec <- list(type = "basis.pca", params = pars)
  pipe$add_step(step_spec)
  pipe
}
