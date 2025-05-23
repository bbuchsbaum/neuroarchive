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

lna_pipeline <- R6::R6Class(
  "lna_pipeline",
  public = list(
    #' @field input Data object or list of run data
    input = NULL,
    #' @field input_summary Character summary of input dimensions
    input_summary = "",
    #' @field runs Character vector of run identifiers
    runs = character(),
    #' @field steps List of transform step specifications
    steps = list(),
    #' @field engine_opts Optional list of hints for core_write
    engine_opts = list(),

    #' @description
    #' Initialise a new lna_pipeline object
    initialize = function() {
      self$input <- NULL
      self$input_summary <- ""
      self$runs <- character()
      self$steps <- list()
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
      if (!is.list(step_spec) || is.null(step_spec$type)) {
        abort_lna(
          "step_spec must be a list with element `type`",
          .subclass = "lna_error_validation",
          location = "lna_pipeline:add_step"
        )
      }

      if (!is.character(step_spec$type) || length(step_spec$type) != 1) {
        abort_lna(
          "step_spec$type must be a single character string",
          .subclass = "lna_error_validation",
          location = "lna_pipeline:add_step"
        )
      }

      if (!is.null(step_spec$params) && !is.list(step_spec$params)) {
        abort_lna(
          "step_spec$params must be a list or NULL",
          .subclass = "lna_error_validation",
          location = "lna_pipeline:add_step"
        )
      }

      if (is.null(step_spec$params)) step_spec$params <- list()

      self$steps[[length(self$steps) + 1]] <- step_spec
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
      step_count <- length(self$steps)
      cat("  Steps:", step_count, "\n")

      if (step_count > 0) {
        for (i in seq_along(self$steps)) {
          step <- self$steps[[i]]
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

          cat(sprintf("  %d: %s [%s]\n", i, type,
                      paste(param_text, collapse = ", ")))
        }
      }

      invisible(self)
    },

    #' @description
    #' Return the internal list of step specifications
    get_steps_list = function() {
      self$steps
    },

    #' @description
    #' Retrieve a step specification by index or by type name. If a type
    #' string is provided and occurs multiple times, the last matching step
    #' is returned. Returns `NULL` if no matching step exists.
    #' @param index_or_type Integer index or character type string
    get_step = function(index_or_type) {
      if (is.numeric(index_or_type)) {
        idx <- as.integer(index_or_type[1])
        if (idx < 1 || idx > length(self$steps)) {
          return(NULL)
        }
        return(self$steps[[idx]])
      } else if (is.character(index_or_type)) {
        typ <- as.character(index_or_type[1])
        matches <- which(vapply(self$steps, function(s) identical(s$type, typ), logical(1)))
        if (length(matches) == 0) {
          return(NULL)
        }
        return(self$steps[[matches[length(matches)]]])
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
      if (length(self$steps) > 0) {
        self$steps[[length(self$steps)]]
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

      find_idx <- function(key) {
        if (is.numeric(key)) {
          idx <- as.integer(key[1])
          if (idx < 1 || idx > length(self$steps)) return(NA_integer_)
          idx
        } else if (is.character(key)) {
          typ <- as.character(key[1])
          matches <- which(vapply(self$steps, function(s) identical(s$type, typ), logical(1)))
          if (length(matches) == 0) return(NA_integer_)
          matches[length(matches)]
        } else {
          abort_lna(
            "index_or_type must be numeric or character",
            .subclass = "lna_error_validation",
            location = "lna_pipeline:modify_step"
          )
        }
      }

      idx <- find_idx(index_or_type)
      if (is.na(idx)) {
        abort_lna(
          "Specified step not found",
          .subclass = "lna_error_validation",
          location = "lna_pipeline:modify_step"
        )
      }

      step <- self$steps[[idx]]
      merged <- utils::modifyList(step$params, new_params_list)
      merged <- merged[!vapply(merged, is.null, logical(1))]

      base <- utils::modifyList(
        default_params(step$type),
        lna_options(step$type)[[step$type]] %||% list()
      )
      step$params <- utils::modifyList(base, merged)

      self$steps[[idx]] <- step
      invisible(self)
    },

    #' @description
    #' Remove a step from the pipeline.
    #' @param index_or_type Integer index or type string identifying the step.
    remove_step = function(index_or_type) {
      find_idx <- function(key) {
        if (is.numeric(key)) {
          idx <- as.integer(key[1])
          if (idx < 1 || idx > length(self$steps)) return(NA_integer_)
          idx
        } else if (is.character(key)) {
          typ <- as.character(key[1])
          matches <- which(vapply(self$steps, function(s) identical(s$type, typ), logical(1)))
          if (length(matches) == 0) return(NA_integer_)
          matches[length(matches)]
        } else {
          abort_lna(
            "index_or_type must be numeric or character",
            .subclass = "lna_error_validation",
            location = "lna_pipeline:remove_step"
          )
        }
      }

      idx <- find_idx(index_or_type)
      if (is.na(idx)) {
        abort_lna(
          "Specified step not found",
          .subclass = "lna_error_validation",
          location = "lna_pipeline:remove_step"
        )
      }

      self$steps[[idx]] <- NULL
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
      if (!is.list(step_spec) || is.null(step_spec$type)) {
        abort_lna(
          "step_spec must be a list with element `type`",
          .subclass = "lna_error_validation",
          location = "lna_pipeline:insert_step"
        )
      }

      find_idx <- function(key) {
        if (is.numeric(key)) {
          idx <- as.integer(key[1])
          if (idx < 1 || idx > length(self$steps)) return(NA_integer_)
          idx
        } else if (is.character(key)) {
          typ <- as.character(key[1])
          matches <- which(vapply(self$steps, function(s) identical(s$type, typ), logical(1)))
          if (length(matches) == 0) return(NA_integer_)
          matches[length(matches)]
        } else {
          abort_lna(
            "index_or_type must be numeric or character",
            .subclass = "lna_error_validation",
            location = "lna_pipeline:insert_step"
          )
        }
      }

      if (!is.null(after_index_or_type)) {
        idx <- find_idx(after_index_or_type)
        if (is.na(idx)) {
          abort_lna(
            "Specified step not found",
            .subclass = "lna_error_validation",
            location = "lna_pipeline:insert_step"
          )
        }
        self$steps <- append(self$steps, list(step_spec), after = idx)
      } else if (!is.null(before_index_or_type)) {
        idx <- find_idx(before_index_or_type)
        if (is.na(idx)) {
          abort_lna(
            "Specified step not found",
            .subclass = "lna_error_validation",
            location = "lna_pipeline:insert_step"
          )
        }
        self$steps <- append(self$steps, list(step_spec), after = idx - 1L)
      } else {
        self$steps <- append(self$steps, list(step_spec), after = length(self$steps))
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

      pkgs <- unique(c("neuroarchive", loadedNamespaces()))

      for (i in seq_along(self$steps)) {
        step <- self$steps[[i]]
        type <- step$type
        params <- step$params %||% list()

        schema_path <- ""
        for (pkg in pkgs) {
          p <- system.file("schemas", paste0(type, ".schema.json"), package = pkg)
          if (nzchar(p) && file.exists(p)) {
            schema_path <- p
            break
          }
        }

        if (!nzchar(schema_path)) {
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
      for (i in seq_along(self$steps)) {
        step <- self$steps[[i]]
        lbl <- sprintf("%d: %s", i, step$type)
        psum <- param_summary(step$params %||% list())
        if (nzchar(psum)) lbl <- paste(lbl, psum, sep = "\n")
        nodes[[length(nodes) + 1]] <- sprintf("n%d [label=\"%s\"];", i, lbl)
      }
      out_idx <- length(self$steps) + 1L
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

