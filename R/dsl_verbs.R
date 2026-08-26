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

#' Resolve parameters for pipeline steps
#' 
#' Internal helper to merge schema defaults, global options, and user parameters
#' 
#' @param transform_type Character string identifying the transform type
#' @param user_params Named list of user-supplied parameters
#' @return Merged parameter list
#' @keywords internal
resolve_params <- function(transform_type, user_params) {
  pars <- default_params(transform_type)
  opts <- lna_options(transform_type)[[transform_type]] %||% list()
  pars <- utils::modifyList(pars, opts)
  pars <- utils::modifyList(pars, user_params)
  pars
}

#' Append a transform step to a pipeline
#'
#' Internal helper used by DSL verbs to create or coerce a pipeline,
#' resolve parameters and add the step.
#'
#' @param data_or_pipe Data object or `lna_pipeline`.
#' @param type Transform type string.
#' @param param_list Named list of user parameters.
#' @return An updated `lna_pipeline` object.
#' @keywords internal
append_lna_step <- function(data_or_pipe, type, param_list = list()) {
  pipe <- if (inherits(data_or_pipe, "lna_pipeline")) {
    data_or_pipe
  } else {
    as_pipeline(data_or_pipe)
  }

  param_list <- param_list[!vapply(param_list, is.null, logical(1))]
  pars <- resolve_params(type, param_list)
  step_spec <- list(type = type, params = pars)
  pipe$add_step(step_spec)
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

  transform_types <- vapply(pipeline_obj$step_list, function(s) s$type, character(1))
  transform_params_list <- list()
  for (i in seq_along(pipeline_obj$step_list)) {
    step <- pipeline_obj$step_list[[i]]
    step_key <- sprintf("step_%02d", i)
    # Keep positional keys for per-step specificity and duplicate transform types.
    transform_params_list[[step_key]] <- step$params
    # Also expose a type-key alias for backward compatibility with older callers/tests.
    transform_params_list[[step$type]] <- step$params
  }

  transform_params_list <- .normalise_transform_pipeline(
    transform_types,
    transform_params_list,
    warn_semantics = TRUE
  )

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
#' @param ... Additional parameters for the quant transform.
#'
#' @return An `lna_pipeline` object with the quant step appended.
#'
#' @details Named arguments passed through `...` may include:
#'   `method` for scale estimation (`"range"` or `"sd"`),
#'   `center` to center values before quantisation,
#'   `scale_scope` to choose global versus voxel-wise scaling,
#'   and `allow_clip` to permit clipping beyond the configured abort threshold.
#'
#' @examples
#' # allow over 5% clipping without error
#' pipe <- as_pipeline(matrix(rnorm(10), 5, 2))
#' pipe <- quant(pipe, bits = 4, allow_clip = TRUE)
#'
#' @export
quant <- function(data_or_pipe, bits = NULL, ...) {
  append_lna_step(
    data_or_pipe,
    "quant",
    c(list(bits = bits), list(...))
  )
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
  user_params <- c(list(k = k), list(...))
  user_params <- user_params[!vapply(user_params, is.null, logical(1))]
  forced_params <- utils::modifyList(list(method = "pca"), user_params)
  append_lna_step(data_or_pipe, "basis", forced_params)
}


##' Infer Embed Step Type and Default Path
#'
#' Helper used by `embed()` to derive the appropriate transform type and
#' default `basis_path` based on the preceding pipeline step. If no
#' known basis-producing step is detected, returns a generic embed type
#' with a `NULL` path.
#'
#' @param prev_step Step specification from `get_last_step_spec()`.
#' @param prev_index Integer index (1-based) of the previous step.
#' @return List with fields `type` and `params`.
#' @keywords internal
infer_embed_step <- function(prev_step, prev_index) {
  res <- list(type = "embed", params = list(basis_path = NULL))
  if (is.null(prev_step)) return(res)

  zero_idx <- prev_index - 1L

  if (identical(prev_step$type, "basis")) {
    method <- prev_step$params$method %||% "pca"
    res$type <- paste0("embed.", method)
    base_name <- sprintf("%02d_%s", zero_idx, prev_step$type)
    res$params$basis_path <- paste0("/basis/", base_name, "/matrix")
    if (isTRUE(prev_step$params$center %||% FALSE)) {
      res$params$center_data_with <- paste0("/basis/", base_name, "/center")
    }
    if (isTRUE(prev_step$params$scale %||% FALSE)) {
      res$params$scale_data_with <- paste0("/basis/", base_name, "/scale")
    }
    if (!is.null(prev_step$params$storage_order)) {
      res$params$storage_order <- prev_step$params$storage_order
    }
  } else if (identical(prev_step$type, "spat.hrbf")) {
    res$type <- "embed.hrbf_analytic"
    base_name <- sprintf("%02d_%s", zero_idx, prev_step$type)
    res$params$basis_path <- paste0("/basis/", base_name, "/matrix")
  }

  res
}

.basis_step_has_following_embed <- function(transforms, idx) {
  next_basis <- which(transforms == "basis" & seq_along(transforms) > idx)
  stop_idx <- if (length(next_basis) > 0) next_basis[[1]] - 1L else length(transforms)
  if (stop_idx <= idx) {
    return(FALSE)
  }
  window <- transforms[seq.int(idx + 1L, stop_idx)]
  any(startsWith(window, "embed"))
}

.basis_step_window <- function(transforms, idx) {
  next_basis <- which(transforms == "basis" & seq_along(transforms) > idx)
  stop_idx <- if (length(next_basis) > 0) next_basis[[1]] - 1L else length(transforms)
  if (stop_idx <= idx) {
    return(character(0))
  }
  transforms[seq.int(idx + 1L, stop_idx)]
}

.basis_without_embed_mode <- function(mode = NULL) {
  if (!is.null(mode)) {
    return(match.arg(mode, c("error", "warn", "allow")))
  }

  explicit_mode <- lna_options("basis_without_embed")[[1]]
  if (!is.null(explicit_mode)) {
    return(match.arg(explicit_mode, c("error", "warn", "allow")))
  }

  if (isTRUE(lna_options("warn.basis_without_embed")[[1]])) {
    "warn"
  } else {
    "allow"
  }
}

.basis_without_embed_allowed <- function(merged, idx) {
  step_key <- sprintf("step_%02d", idx)
  step_params <- merged[[step_key]] %||% list()
  isTRUE(step_params$allow_dense_passthrough)
}

.strip_semantic_override_params <- function(merged) {
  for (nm in names(merged)) {
    if (is.list(merged[[nm]]) && !is.null(merged[[nm]]$allow_dense_passthrough)) {
      merged[[nm]]$allow_dense_passthrough <- NULL
    }
  }
  merged
}

.validate_basis_without_embed <- function(transforms, merged, mode = NULL) {
  mode <- .basis_without_embed_mode(mode)
  if (identical(mode, "allow")) {
    return(invisible(NULL))
  }

  basis_idx <- which(transforms == "basis")
  if (length(basis_idx) == 0) {
    return(invisible(NULL))
  }

  for (idx in basis_idx) {
    window <- .basis_step_window(transforms, idx)

    if (.basis_step_has_following_embed(transforms, idx) ||
        length(window) == 0L ||
        .basis_without_embed_allowed(merged, idx)) {
      next
    }

    msg <- sprintf(
      "Step %d uses 'basis' without a following 'embed'. This stores basis metadata but keeps dense voxel data; add 'embed' after 'basis' for coefficient-based latent compression. To allow this dense pass-through pipeline explicitly, set basis$params$allow_dense_passthrough = TRUE or lna_options(basis_without_embed = 'allow').",
      idx
    )

    if (identical(mode, "warn")) {
      warn_lna(
        msg,
        .subclass = "lna_warning_validation",
        location = "pipeline_semantics:basis_without_embed"
      )
    } else {
      abort_lna(
        msg,
        .subclass = "lna_error_validation",
        location = "pipeline_semantics:basis_without_embed"
      )
    }
  }

  invisible(NULL)
}

.param_is_missing <- function(x) {
  is.null(x) || (is.character(x) && length(x) == 1L && !nzchar(x))
}

.normalise_transform_pipeline <- function(transforms, transform_params = list(),
                                         warn_semantics = TRUE) {
  checked <- isTRUE(attr(transform_params, "lna_semantics_checked"))
  merged <- resolve_transform_params(transforms, transform_params)

  for (i in seq_along(transforms)) {
    if (i < 2L || !startsWith(transforms[[i]], "embed")) {
      next
    }

    step_key <- sprintf("step_%02d", i)
    prev_key <- sprintf("step_%02d", i - 1L)
    prev_step <- list(
      type = transforms[[i - 1L]],
      params = merged[[prev_key]] %||% list()
    )
    info <- infer_embed_step(prev_step, i - 1L)
    inferred <- info$params %||% list()
    if (length(inferred) == 0) {
      next
    }

    step_params <- merged[[step_key]] %||% list()
    for (nm in names(inferred)) {
      if (.param_is_missing(step_params[[nm]]) && !is.null(inferred[[nm]])) {
        step_params[[nm]] <- inferred[[nm]]
      }
    }
    merged[[step_key]] <- step_params
    merged[[transforms[[i]]]] <- step_params
  }

  if (isTRUE(warn_semantics) && !checked) {
    .validate_basis_without_embed(transforms, merged)
  }

  merged <- .strip_semantic_override_params(merged)
  attr(merged, "lna_semantics_checked") <- TRUE
  merged
}


##' Embed DSL verb
#'
#' Adds an embedding step that projects the input data onto a basis
#' computed by a preceding transform. When called immediately after a
#' `basis` step (e.g., created by `pca()`), the path to the basis matrix
#' is inferred automatically using the conventional HDF5 location
#' `/basis/<NN>_basis/matrix` where `<NN>` is the zero-based index of the
#' previous step. When the preceding step is `basis`, compatible
#' centering, scaling, and storage-order metadata are also inferred.
#'
#' @param data_or_pipe Data object or `lna_pipeline`.
#' @param basis_path Optional explicit HDF5 path to the basis matrix.
#' @param basis_step Optional step index or type to use for basis inference
#'   instead of the immediately preceding step.
#' @param ... Additional parameters for the embed transform.
#'
#' @return An `lna_pipeline` object with the embed step appended.
#' @export
embed <- function(data_or_pipe, basis_path = NULL, basis_step = NULL, ...) {
  pipe <- if (inherits(data_or_pipe, "lna_pipeline")) {
    data_or_pipe
  } else {
    as_pipeline(data_or_pipe)
  }

  # Determine which step to use for basis inference
  if (!is.null(basis_step)) {
    prev <- pipe$get_step(basis_step)
    if (is.null(prev)) {
      abort_lna(
        sprintf("basis_step '%s' not found in pipeline", basis_step),
        .subclass = "lna_error_validation",
        location = "embed:basis_step"
      )
    }
    # Find the index of this step
    step_idx <- if (is.numeric(basis_step)) {
      as.integer(basis_step[1])
    } else {
      # Find last matching step by type
      steps <- pipe$step_list
      matches <- which(vapply(steps, function(s) identical(s$type, basis_step), logical(1)))
      if (length(matches) == 0) NA_integer_ else matches[length(matches)]
    }
  } else {
    prev <- pipe$get_last_step_spec()
    step_idx <- length(pipe$step_list)
  }
  
  if (is.null(prev)) {
    abort_lna(
      "embed() must follow or reference a basis-producing step",
      .subclass = "lna_error_validation",
      location = "embed:context"
    )
  }

  user_params <- c(list(basis_path = basis_path), list(...))
  user_params <- user_params[!vapply(user_params, is.null, logical(1))]

  info <- infer_embed_step(prev, step_idx)
  embed_type <- info$type
  inferred_params <- info$params %||% list()

  for (nm in names(inferred_params)) {
    if (is.null(user_params[[nm]]) && !is.null(inferred_params[[nm]])) {
      user_params[[nm]] <- inferred_params[[nm]]
    }
  }

  if (is.null(user_params$basis_path)) {
    abort_lna(
      "basis_path must be supplied or inferable from previous step",
      .subclass = "lna_error_validation",
      location = "embed:basis_path"
    )
  }

  append_lna_step(pipe, embed_type, user_params)
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
  append_lna_step(
    data_or_pipe,
    "delta",
    c(list(order = order), list(...))
  )
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
  append_lna_step(
    data_or_pipe,
    "temporal",
    c(list(kind = kind), list(...))
  )
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
  append_lna_step(
    data_or_pipe,
    "spat.hrbf",
    c(list(levels = levels), list(...))
  )
}

# S3 methods for seamless verb chaining
# These allow using verbs directly on data without explicit as_pipeline() calls

#' @rdname quant
#' @method quant default
#' @export
quant.default <- quant

#' @rdname pca
#' @method pca default
#' @export
pca.default <- pca

#' @rdname embed
#' @method embed default
#' @export
embed.default <- embed

#' @rdname delta
#' @method delta default
#' @export
delta.default <- delta

#' @rdname temporal
#' @method temporal default
#' @export
temporal.default <- temporal

#' @rdname hrbf
#' @method hrbf default
#' @export
hrbf.default <- hrbf
