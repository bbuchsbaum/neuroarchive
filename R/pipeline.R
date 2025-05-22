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
    }
  )
)

