#' LNAFacade class
#'
#' A lightweight wrapper around \code{write_lna()} and \code{read_lna()}.
#' Default transform parameters can be supplied when creating the object.
#' The path of the last written file is stored in \code{last_output}.
#'
#' @examples
#' fac <- LNAFacade$new()
#' tmp <- tempfile(fileext = ".h5")
#' fac$write(array(1, dim = c(1, 1, 1)), tmp, transforms = character())
#' fac$read(tmp)
#' @export
LNAFacade <- R6::R6Class(
  "LNAFacade",
  public = list(
    #' @field default_transform_params Default transform parameters
    default_transform_params = NULL,
    #' @field last_output Path of the last written file
    last_output = NULL,

    #' @description
    #' Create a new LNAFacade
    #' @param transform_params Named list of default transform parameters
    initialize = function(transform_params = list()) {
      stopifnot(is.list(transform_params))
      self$default_transform_params <- transform_params
    },

    #' @description
    #' Replace the default transform parameters
    #' @param params Named list of default transform parameters
    #' @return The \code{LNAFacade} object, invisibly
    set_defaults = function(params) {
      stopifnot(is.list(params))
      self$default_transform_params <- params
      invisible(self)
    },

    #' @description
    #' Write data to an LNA file
    #' @param x Array or list of arrays
    #' @param file Output path
    #' @param transforms Character vector of transforms
    #' @param transform_params Optional named list overriding defaults
    #' @param ... Additional arguments forwarded to \code{write_lna()}
    #' @return Invisibly returns the \code{lna_write_result} from
    #'   \code{write_lna()}, including elements \code{file}, \code{plan},
    #'   and \code{header}.
    write = function(x, file, transforms, transform_params = NULL, ...) {
      params <- utils::modifyList(
        self$default_transform_params,
        transform_params %||% list(),
        keep.null = TRUE
      )
      res <- write_lna(x = x, file = file, transforms = transforms,
                       transform_params = params, ...)
      self$last_output <- res$file
      invisible(res)
    },

    #' @description
    #' Read data from an LNA file
    #' @param file Path to an LNA file. If missing, uses \code{last_output}
    #'   from the most recent call to \code{$write()}.
    #' @param ... Arguments forwarded to \code{read_lna()}
    read = function(file = self$last_output, ...) {
      read_lna(file = file, ...)
    },

    #' @description
    #' Execute a pipeline and write an LNA file
    #' @param pipe An \code{lna_pipeline} object
    #' @param file Output path
    #' @param ... Additional arguments forwarded to \code{lna_write()}
    write_pipeline = function(pipe, file, ...) {
      res <- lna_write(pipe, file = file, ...)
      self$last_output <- res$file
      invisible(res$file)
    }
  )
)
