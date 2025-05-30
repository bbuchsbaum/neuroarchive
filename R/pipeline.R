#' lna_pipeline Class
#'
#' @description
#' Basic R6 class for constructing LNA pipelines. Stores input data,
#' pipeline steps and optional engine hints. This is an early draft
#' used for experimenting with a tidy DSL facade.
#'
#' @importFrom R6 R6Class
#' @keywords internal
style_subtle <- function(x) {
  if (requireNamespace("pillar", quietly = TRUE)) {
    pillar::style_subtle(x)
  } else {
    x
  }
}
style_bold <- function(x) {
  if (requireNamespace("pillar", quietly = TRUE)) {
    pillar::style_bold(x)
  } else {
    x
  }
}

#' Schema path cache environment
#'
#' Used by `schema_path_find()` to avoid repeatedly scanning installed
#' packages for JSON schemas.
#' @keywords internal
.schema_path_cache <- new.env(parent = emptyenv())

#' Locate a transform schema
#'
#' Searches loaded namespaces for a JSON schema matching `type` and
#' caches the discovered path for future lookups.
#'
#' @param type Character transform type.
#' @return File path string or `NA_character_` when not found.
#' @keywords internal
schema_path_find <- function(type) {
  cache <- .schema_path_cache
  if (exists(type, envir = cache, inherits = FALSE)) {
    return(cache[[type]])
  }

  pkgs <- unique(c("neuroarchive", loadedNamespaces()))
  schema_path <- NA_character_
  for (pkg in pkgs) {
    p <- system.file("schemas", paste0(type, ".schema.json"), package = pkg)
    if (nzchar(p) && file.exists(p)) {
      schema_path <- p
      break
    }
  }

  assign(type, schema_path, envir = cache)
  schema_path
}

lna_pipeline <- R6::R6Class(
  "lna_pipeline",
  public = list(
    #' @field input Data object or list of run data
    input = NULL,
    #' @field input_summary Character summary of input dimensions
    input_summary = "",
    #' @field runs Character vector of run identifiers
    runs = character(),
    #' @field step_list List of transform step specifications
    step_list = list(),
    #' @field engine_opts Optional list of hints for core_write
    engine_opts = list(),

    #' @description
    #' Initialise a new lna_pipeline object
    initialize = function() {
      self$input <- NULL
      self$input_summary <- ""
      self$runs <- character()
      self$step_list <- list()
      self$engine_opts <- list()
    },

    #' @description
    #' Set the pipeline input and related metadata
    #' @param x Data object or list of run data
    #' @param run_ids Optional character vector of run identifiers
    #' @param chunk_mb_suggestion Optional numeric hint for chunk size
    set_input = function(x, run_ids = NULL, chunk_mb_suggestion = NULL) {
      if (is.null(x)) {
        abort_lna(
          "input `x` must not be NULL",
          .subclass = "lna_error_validation",
          location = "lna_pipeline:set_input"
        )
      }

      validate_single <- function(obj) {
        if (!(is.array(obj) || is.matrix(obj) || methods::is(obj, "NeuroVec"))) {
          abort_lna(
            "input must be array, matrix, NeuroVec or list of such objects",
            .subclass = "lna_error_validation",
            location = "lna_pipeline:set_input"
          )
        }
      }

      if (is.list(x) && !methods::is(x, "NeuroVec")) {
        if (length(x) == 0) {
          abort_lna(
            "input list must contain at least one element",
            .subclass = "lna_error_validation",
            location = "lna_pipeline:set_input"
          )
        }
        lapply(x, validate_single)

        ref_dim <- dim(x[[1]])
        if (!all(vapply(x, function(el) identical(dim(el), ref_dim), logical(1)))) {
          abort_lna(
            "all input elements must have identical dimensions",
            .subclass = "lna_error_validation",
            location = "lna_pipeline:set_input"
          )
        }

        run_count <- length(x)
        if (is.null(run_ids)) {
          if (!is.null(names(x)) && all(names(x) != "")) {
            self$runs <- names(x)
          } else {
            self$runs <- sprintf("run-%02d", seq_len(run_count))
          }
        } else {
          run_ids <- as.character(run_ids)
          if (length(run_ids) != run_count) {
            abort_lna(
              "length of run_ids must match number of list elements",
              .subclass = "lna_error_validation",
              location = "lna_pipeline:set_input"
            )
          }
          self$runs <- run_ids
        }
        exemplar <- x[[1]]
      } else {
        validate_single(x)
        run_count <- 1L
        self$runs <- if (is.null(run_ids)) "run-01" else as.character(run_ids[1])
        exemplar <- x
      }

      dims <- dim(exemplar)
      if (is.null(dims)) {
        time_dim <- length(exemplar)
        vox_dim <- 1L
      } else {
        time_dim <- dims[length(dims)]
        vox_dim <- if (length(dims) > 1) prod(dims[-length(dims)]) else dims[1]
      }

      plural <- if (run_count == 1L) "" else "s"
      self$input_summary <- sprintf(
        "%d run%s × (%d TR × %s vox)",
        run_count, plural, as.integer(time_dim), format(as.integer(vox_dim), scientific = FALSE)
      )

      self$input <- x
      if (!is.null(chunk_mb_suggestion)) {
        if (!is.numeric(chunk_mb_suggestion) ||
            length(chunk_mb_suggestion) != 1 ||
            chunk_mb_suggestion <= 0) {
          abort_lna(
            "chunk_mb_suggestion must be a single positive number",
            .subclass = "lna_error_validation",
            location = "lna_pipeline:set_input"
          )
        }
        self$engine_opts$chunk_mb_suggestion <- chunk_mb_suggestion
      } else {
        self$engine_opts$chunk_mb_suggestion <- NULL
      }

      invisible(self)
    },

    #' @description
    #' Append a transform step specification to the pipeline
    #' @param step_spec A list with elements `type` and `params`
    add_step = function(step_spec) {
      step_spec <- validate_step_spec(step_spec, "lna_pipeline:add_step")
      self$step_list[[length(self$step_list) + 1]] <- step_spec
      invisible(self)
    },

    #' @description
    #' Print a human readable summary of the pipeline
    print = function(...) {
      cat("<lna_pipeline>\n")
      if (nzchar(self$input_summary)) {
        cat("  Input:", self$input_summary, "\n")
      } else {
        cat("  Input: (not set)\n")
      }
      step_count <- length(self$step_list)
      cat("  Steps:", step_count, "\n")

      if (step_count > 0) {
        for (i in seq_along(self$step_list)) {
          step <- self$step_list[[i]]
          type <- step$type
          params <- step$params %||% list()

          defaults <- utils::modifyList(
            default_params(type),
            lna_options(type)[[type]] %||% list()
          )

          param_text <- vapply(names(params), function(nm) {
            val <- params[[nm]]
            val_str <- paste0(nm, "=", format(val))
            def <- defaults[[nm]]
            if (!is.null(def) && identical(val, def)) {
              style_subtle(val_str)
            } else {
              style_bold(val_str)
            }
          }, character(1))

          # Clip long parameter lists for readability
          param_line <- paste(param_text, collapse = ", ")
          if (nchar(param_line) > 50) {
            param_line <- paste0(substr(param_line, 1, 47), style_subtle("..."))
          }

          cat(sprintf("  %d: %s [%s]\n", i, type, param_line))
        }
      }

      invisible(self)
    },

    #' @description
    #' Return the internal list of step specifications
    get_steps_list = function() {
      self$step_list
    },

    #' @description
    #' Return the internal list of step specifications
    #' @return List of step specifications
    steps = function() {
      self$step_list
    },

    #' @description
    #' Retrieve a step specification by index or by type name. If a type
    #' string is provided and occurs multiple times, the last matching step
    #' is returned. Returns `NULL` if no matching step exists.
    #' @param index_or_type Integer index or character type string
    get_step = function(index_or_type) {
      if (is.numeric(index_or_type)) {
        idx <- as.integer(index_or_type[1])
        if (idx < 1 || idx > length(self$step_list)) {
          return(NULL)
        }
        return(self$step_list[[idx]])
      } else if (is.character(index_or_type)) {
        typ <- as.character(index_or_type[1])
        matches <- which(vapply(self$step_list, function(s) identical(s$type, typ), logical(1)))
        if (length(matches) == 0) {
          return(NULL)
        }
        return(self$step_list[[matches[length(matches)]]])
      } else {
        abort_lna(
          "index_or_type must be numeric or character",
          .subclass = "lna_error_validation",
          location = "lna_pipeline:get_step"
        )
      }
    },

    #' @description
    #' Return the specification of the most recently added step, or `NULL`
    #' if no steps have been added yet.
    get_last_step_spec = function() {
      if (length(self$step_list) > 0) {
        self$step_list[[length(self$step_list)]]
      } else {
        NULL
      }
    },

    #' @description
    #' Modify parameters of an existing step.
    #' @param index_or_type Integer index or type string identifying the step.
    #' @param new_params_list Named list of parameter updates. `NULL` values
    #'   remove parameters and revert them to defaults/options.
    modify_step = function(index_or_type, new_params_list) {
      if (!is.list(new_params_list)) {
        abort_lna(
          "new_params_list must be a list",
          .subclass = "lna_error_validation",
          location = "lna_pipeline:modify_step"
        )
      }

      idx <- find_step_index(self$step_list, index_or_type)
      if (is.na(idx)) {
        abort_lna(
          "Specified step not found",
          .subclass = "lna_error_validation",
          location = "lna_pipeline:modify_step"
        )
      }

      step <- self$step_list[[idx]]
      merged <- utils::modifyList(step$params, new_params_list)
      merged <- merged[!vapply(merged, is.null, logical(1))]

      base <- utils::modifyList(
        default_params(step$type),
        lna_options(step$type)[[step$type]] %||% list()
      )
      step$params <- utils::modifyList(base, merged)

      self$step_list[[idx]] <- step
      invisible(self)
    },

    #' @description
    #' Remove a step from the pipeline.
    #' @param index_or_type Integer index or type string identifying the step.
    remove_step = function(index_or_type) {
      idx <- find_step_index(self$step_list, index_or_type)
      if (is.na(idx)) {
        abort_lna(
          "Specified step not found",
          .subclass = "lna_error_validation",
          location = "lna_pipeline:remove_step"
        )
      }

      self$step_list[[idx]] <- NULL
      invisible(self)
    },

    #' @description
    #' Insert a new step at a specific position.
    #' @param step_spec Step specification list with `type` and `params`.
    #' @param after_index_or_type Insert after this step. Mutually exclusive with
    #'   `before_index_or_type`.
    #' @param before_index_or_type Insert before this step.
    insert_step = function(step_spec,
                           after_index_or_type = NULL,
                           before_index_or_type = NULL) {
      if (!is.null(after_index_or_type) && !is.null(before_index_or_type)) {
        abort_lna(
          "Specify only one of after_index_or_type or before_index_or_type",
          .subclass = "lna_error_validation",
          location = "lna_pipeline:insert_step"
        )
      }

      step_spec <- validate_step_spec(step_spec, "lna_pipeline:insert_step")

      if (!is.null(after_index_or_type)) {
        idx <- find_step_index(self$step_list, after_index_or_type)
        if (is.na(idx)) {
          abort_lna(
            "Specified step not found",
            .subclass = "lna_error_validation",
            location = "lna_pipeline:insert_step"
          )
        }
        self$step_list <- append(self$step_list, list(step_spec), after = idx)
      } else if (!is.null(before_index_or_type)) {
        idx <- find_step_index(self$step_list, before_index_or_type)
        if (is.na(idx)) {
          abort_lna(
            "Specified step not found",
            .subclass = "lna_error_validation",
            location = "lna_pipeline:insert_step"
          )
        }
        self$step_list <- append(self$step_list, list(step_spec), after = idx - 1L)
      } else {
        self$step_list <- append(self$step_list, list(step_spec), after = length(self$step_list))
      }
      invisible(self)
    },

    #' @description
    #' Validate all step parameters against their JSON schemas.
    #' @param strict Logical flag. If `TRUE`, abort on the first validation
    #'   failure. If `FALSE` (default), collect all issues and return them.
    validate_params = function(strict = FALSE) {
      stopifnot(is.logical(strict), length(strict) == 1)

      issues <- character()

      fail <- function(msg, type) {
        loc <- sprintf("lna_pipeline:validate_params:%s", type)
        if (strict) {
          abort_lna(msg, .subclass = "lna_error_validation", location = loc)
        } else {
          warning(msg, call. = FALSE)
          issues <<- c(issues, msg)
        }
      }

      for (i in seq_along(self$step_list)) {
        step <- self$step_list[[i]]
        type <- step$type
        params <- step$params %||% list()

        schema_path <- schema_path_find(type)
        if (is.na(schema_path) || !nzchar(schema_path)) {
          fail(sprintf("Schema for transform '%s' not found", type), type)
          next
        }

        json <- jsonlite::toJSON(params, auto_unbox = TRUE)
        valid <- tryCatch(
          jsonvalidate::json_validate(json, schema_path, verbose = TRUE),
          error = function(e) e
        )

        if (!isTRUE(valid)) {
          msg <- sprintf("Step %d (type='%s') parameters failed schema validation", i, type)
          fail(msg, type)
        }
      }

      if (length(issues) == 0) TRUE else issues
    },

    #' @description
    #' Produce a diagram of the pipeline.
    #' @param engine Output engine: one of "grViz", "ascii", or "dot".
    #' @return DOT string, a `DiagrammeR` graph object, or ASCII text.
    diagram = function(engine = c("grViz", "ascii", "dot")) {
      engine <- match.arg(engine)

      clip_text <- function(x, limit = 30) {
        if (nchar(x) > limit) paste0(substr(x, 1, limit), "...") else x
      }

      param_summary <- function(params) {
        if (length(params) == 0) return("")
        kv <- vapply(names(params), function(nm) {
          val <- params[[nm]]
          if (is.atomic(val) && length(val) == 1) {
            paste0(nm, "=", clip_text(as.character(val)))
          } else {
            paste0(nm, "=[...]")
          }
        }, character(1))
        paste(kv, collapse = "\n")
      }

      nodes <- list(sprintf('n0 [label="Input\n%s"];', self$input_summary))
      for (i in seq_along(self$step_list)) {
        step <- self$step_list[[i]]
        lbl <- sprintf("%d: %s", i, step$type)
        psum <- param_summary(step$params %||% list())
        if (nzchar(psum)) lbl <- paste(lbl, psum, sep = "\n")
        nodes[[length(nodes) + 1]] <- sprintf("n%d [label=\"%s\"];", i, lbl)
      }
      out_idx <- length(self$step_list) + 1L
      nodes[[length(nodes) + 1]] <- sprintf("n%d [label=\"Output\"];", out_idx)

      edges <- vapply(seq_len(out_idx), function(j) {
        sprintf("n%d -> n%d;", j - 1L, j)
      }, character(1))

      dot <- paste(c(
        "digraph pipeline {",
        "  rankdir=LR;",
        paste0("  ", unlist(nodes)),
        paste0("  ", edges),
        "}"), collapse = "\n")

      if (engine == "dot") {
        return(dot)
      }

      if (engine == "grViz") {
        if (requireNamespace("DiagrammeR", quietly = TRUE)) {
          return(DiagrammeR::grViz(dot))
        } else {
          warning("DiagrammeR not installed; returning DOT string.", call. = FALSE)
          return(dot)
        }
      }

      if (engine == "ascii") {
        if (requireNamespace("DiagrammeR", quietly = TRUE) &&
            requireNamespace("DiagrammeRsvg", quietly = TRUE) &&
            requireNamespace("asciiSVG", quietly = TRUE)) {
          gr <- DiagrammeR::grViz(dot)
          svg <- DiagrammeRsvg::export_svg(gr)
          asc <- asciiSVG::ascii_svg(svg)
          return(paste(asc, collapse = "\n"))
        } else {
          warning("ASCII engine unavailable; returning DOT string.", call. = FALSE)
          return(dot)
        }
      }
    }
  )
)

#' Find the index of a pipeline step
#'
#' Internal helper used by lna_pipeline methods to resolve a step
#' by numeric position or by its `type` name.
#'
#' @param steps List of step specifications.
#' @param key Numeric index or type string identifying the step.
#' @return Integer index or `NA_integer_` if no match is found.
#' @keywords internal
find_step_index <- function(steps, key) {
  if (is.numeric(key)) {
    idx <- as.integer(key[1])
    if (idx < 1 || idx > length(steps)) return(NA_integer_)
    idx
  } else if (is.character(key)) {
    typ <- as.character(key[1])
    matches <- which(vapply(steps, function(s) identical(s$type, typ), logical(1)))
    if (length(matches) == 0) return(NA_integer_)
    matches[length(matches)]
  } else {
    abort_lna(
      "index_or_type must be numeric or character",
      .subclass = "lna_error_validation"
    )
  }
}

#' Validate a step specification
#'
#' Internal helper used by lna_pipeline methods to verify that a
#' step specification contains the required `type` field and an
#' optional `params` list.
#'
#' @param step_spec List describing the step.
#' @param location Character string identifying the caller for errors.
#' @return The validated (and normalised) step specification.
#' @keywords internal
validate_step_spec <- function(step_spec, location) {
  if (!is.list(step_spec) || is.null(step_spec$type)) {
    abort_lna(
      "step_spec must be a list with element `type`",
      .subclass = "lna_error_validation",
      location = location
    )
  }

  if (!is.character(step_spec$type) || length(step_spec$type) != 1) {
    abort_lna(
      "step_spec$type must be a single character string",
      .subclass = "lna_error_validation",
      location = location
    )
  }

  if (!is.null(step_spec$params) && !is.list(step_spec$params)) {
    abort_lna(
      "step_spec$params must be a list or NULL",
      .subclass = "lna_error_validation",
      location = location
    )
  }

  if (is.null(step_spec$params)) step_spec$params <- list()
  step_spec
}

