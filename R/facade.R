#' LNAFacade helper class
#'
#' @description
#' A lightweight wrapper exposing `compress()` and `open()` methods which
#' simply call [write_lna()] and [read_lna()] respectively.
#'
#' @examples
#' f <- LNAFacade$new()
#' f$compress(array(rnorm(8), dim = c(2,2,2)), "demo.h5", transforms = "quant")
#' f$open("demo.h5")
#' @importFrom R6 R6Class
#'
#' @export
LNAFacade <- R6::R6Class(
  "LNAFacade",
  public = list(
    #' @description Compress data to an LNA file
    #' @param x Input array or list of arrays
    #' @param file Output file path
    #' @param ... Additional arguments passed to [write_lna()]
    compress = function(x, file, ...) {
      write_lna(x, file = file, ...)
    },
    #' @description Open an LNA file
    #' @param file Path to the LNA file
    #' @param ... Additional arguments passed to [read_lna()]
    open = function(file, ...) {
      read_lna(file, ...)
    }
  )
)

