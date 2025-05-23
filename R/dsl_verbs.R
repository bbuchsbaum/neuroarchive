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
#' @param bits Number of quantization bits (1-16). If `NULL`, the
#'   schema default is used.
#' @param method Method for determining scale/offset (`"range"` or
#'   `"sd"`).
#' @param center Logical indicating whether the data should be
#'   effectively centered before quantisation.
#' @param scale_scope Either `"global"` for one scale/offset or
#'   `"voxel"` for per-voxel parameters.
#' @param allow_clip If `TRUE`, quantisation proceeds even when the
#'   clipping percentage exceeds `lna.quant.clip_abort_pct`.
#' @param ... Additional parameters for the quant transform.
#'
#' @return An `lna_pipeline` object with the quant step appended.
#'
#' @examples
#' # allow over 5% clipping without error
#' pipe <- as_pipeline(matrix(rnorm(10), 5, 2))
#' pipe <- quant(pipe, bits = 4, allow_clip = TRUE)
#'
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


##' Principal Component Analysis DSL verb
#'
#' Adds a PCA basis computation step to a pipeline. If `data_or_pipe`
#' is not an `lna_pipeline`, a new pipeline is created via
#' `as_pipeline()`.
#'
#' Parameter values are resolved by merging schema defaults for the
#' `'basis'` transform, global `lna_options("basis")`, the forced
#' `method = "pca"`, and any user-supplied arguments.
#'
#' @param data_or_pipe Data object or `lna_pipeline`.
#' @param k Optional number of principal components.
#' @param ... Additional parameters for the basis transform.
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

  pars <- default_params("basis")
  opts <- lna_options("basis")$basis %||% list()
  pars <- utils::modifyList(pars, opts)
  pars <- utils::modifyList(pars, list(method = "pca"))
  pars <- utils::modifyList(pars, user_params)

  step_spec <- list(type = "basis", params = pars)

  pipe$add_step(step_spec)
  pipe
}


##' Embed DSL verb
#'
#' Adds an embedding step that projects the input data onto a basis
#' computed by a preceding transform. When called immediately after a
#' `basis` step (e.g., created by `pca()`), the path to the basis matrix
#' is inferred automatically using the conventional HDF5 location
#' `/basis/<NN>_basis/matrix` where `<NN>` is the zero-based index of the
#' previous step.
#'
#' @param data_or_pipe Data object or `lna_pipeline`.
#' @param basis_path Optional explicit HDF5 path to the basis matrix.
#' @param ... Additional parameters for the embed transform.
#'
#' @return An `lna_pipeline` object with the embed step appended.
#' @export
embed <- function(data_or_pipe, basis_path = NULL, ...) {
  pipe <- if (inherits(data_or_pipe, "lna_pipeline")) {
    data_or_pipe
  } else {
    as_pipeline(data_or_pipe)
  }

  prev <- pipe$get_last_step_spec()
  if (is.null(prev)) {
    abort_lna(
      "embed() must follow a basis-producing step",
      .subclass = "lna_error_validation",
      location = "embed:context"
    )
  }

  user_params <- c(list(basis_path = basis_path), list(...))
  user_params <- user_params[!vapply(user_params, is.null, logical(1))]

  embed_type <- "embed"
  default_path <- NULL

  prev_index <- length(pipe$steps)
  if (prev$type == "basis") {
    base_name <- sprintf("%02d_%s", prev_index - 1L, prev$type)
    default_path <- paste0("/basis/", base_name, "/matrix")
  }

  pars <- default_params("embed")
  opts <- lna_options("embed")$embed %||% list()
  pars <- utils::modifyList(pars, opts)

  if (is.null(user_params$basis_path)) {
    if (!is.null(default_path)) {
      user_params$basis_path <- default_path
    } else {
      abort_lna(
        "basis_path must be supplied or inferable from previous step",
        .subclass = "lna_error_validation",
        location = "embed:basis_path"
      )
    }
  }

  pars <- utils::modifyList(pars, user_params)

  step_spec <- list(type = embed_type, params = pars)
  pipe$add_step(step_spec)
  pipe
}

##' Finite Difference (Delta) DSL verb
#'
#' Adds a delta encoding step to a pipeline. If `data_or_pipe`
#' is not an `lna_pipeline`, a new pipeline is created via
#' `as_pipeline()`.
#'
#' Parameter values are resolved by merging schema defaults for
#' the `delta` transform, global `lna_options("delta")`, and any
#' user-supplied arguments.
#'
#' @param data_or_pipe Data object or `lna_pipeline`.
#' @param order Optional difference order.
#' @param ... Additional parameters for the delta transform.
#'
#' @return An `lna_pipeline` object with the delta step appended.
#' @export
delta <- function(data_or_pipe, order = NULL, ...) {

  pipe <- if (inherits(data_or_pipe, "lna_pipeline")) {
    data_or_pipe
  } else {
    as_pipeline(data_or_pipe)
  }

  user_params <- c(list(order = order), list(...))
  user_params <- user_params[!vapply(user_params, is.null, logical(1))]

  pars <- default_params("delta")
  opts <- lna_options("delta")$delta %||% list()
  pars <- utils::modifyList(pars, opts)
  pars <- utils::modifyList(pars, user_params)

  step_spec <- list(type = "delta", params = pars)
  pipe$add_step(step_spec)
  pipe
}

##' Temporal Basis Projection DSL verb
#'
#' Adds a temporal basis transform step to a pipeline. If
#' `data_or_pipe` is not an `lna_pipeline`, a new pipeline is
#' created via `as_pipeline()`.
#'
#' Parameter values are resolved by merging schema defaults for
#' the `temporal` transform, global `lna_options("temporal")`, and
#' any user-supplied arguments.
#'
#' @param data_or_pipe Data object or `lna_pipeline`.
#' @param kind Optional temporal basis type (e.g., "dct").
#' @param ... Additional parameters for the temporal transform.
#'
#' @return An `lna_pipeline` object with the temporal step appended.
#' @export
temporal <- function(data_or_pipe, kind = NULL, ...) {
  pipe <- if (inherits(data_or_pipe, "lna_pipeline")) {
    data_or_pipe
  } else {
    as_pipeline(data_or_pipe)
  }

  user_params <- c(list(kind = kind), list(...))
  user_params <- user_params[!vapply(user_params, is.null, logical(1))]

  pars <- default_params("temporal")
  opts <- lna_options("temporal")[["temporal"]] %||% list()
  pars <- utils::modifyList(pars, opts)
  pars <- utils::modifyList(pars, user_params)

  step_spec <- list(type = "temporal", params = pars)
  pipe$add_step(step_spec)
  pipe
}

##' Hierarchical Radial Basis Function DSL verb
#'
#' Adds a spatial HRBF basis generation step to a pipeline. If
#' `data_or_pipe` is not an `lna_pipeline`, a new pipeline is
#' created via `as_pipeline()`.
#'
#' Parameter values are resolved by merging schema defaults for
#' the `spat.hrbf` transform, global `lna_options("spat.hrbf")`, and
#' any user-supplied arguments.
#'
#' @param data_or_pipe Data object or `lna_pipeline`.
#' @param levels Optional number of HRBF resolution levels.
#' @param ... Additional parameters for the HRBF transform.
#'
#' @return An `lna_pipeline` object with the HRBF step appended.
#' @export
hrbf <- function(data_or_pipe, levels = NULL, ...) {
  pipe <- if (inherits(data_or_pipe, "lna_pipeline")) {
    data_or_pipe
  } else {
    as_pipeline(data_or_pipe)
  }

  user_params <- c(list(levels = levels), list(...))
  user_params <- user_params[!vapply(user_params, is.null, logical(1))]

  pars <- default_params("spat.hrbf")
  opts <- lna_options("spat.hrbf")[["spat.hrbf"]] %||% list()
  pars <- utils::modifyList(pars, opts)
  pars <- utils::modifyList(pars, user_params)

  step_spec <- list(type = "spat.hrbf", params = pars)
  pipe$add_step(step_spec)
  pipe
}
