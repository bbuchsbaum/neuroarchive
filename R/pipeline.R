#' lna_pipeline Class
#'
#' @description
#' Basic R6 class for constructing LNA pipelines. Stores input data,
#' pipeline steps and optional engine hints. This is an early draft
#' used for experimenting with a tidy DSL facade.
#'
#' @importFrom R6 R6Class
#' @keywords internal
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
        if (!(is.array(obj) || is.matrix(obj) || inherits(obj, "NeuroVec"))) {
          abort_lna(
            "input must be array, matrix, NeuroVec or list of such objects",
            .subclass = "lna_error_validation",
            location = "lna_pipeline:set_input"
          )
        }
      }

      if (is.list(x) && !inherits(x, "NeuroVec")) {
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
    }
  )
)

