#' S3 Dispatch for Transform Steps
#'
#' @description Defines S3 generic functions for forward and inverse transform steps.
#'   Methods should be implemented for specific transform types.
#'
#' @keywords internal

#' Apply a forward transform step.
#'
#' @param type (character) The type identifier of the transform (e.g., "mask", "pca").
#' @param desc (list) The parsed JSON descriptor for this transform step.
#' @param handle (DataHandle) The current data handle.
#'
#' @return An updated `DataHandle` object after applying the forward step.
#' @export
forward_step <- function(type, desc, handle) {
  UseMethod("forward_step", type)
}

#' Default method for forward_step.
#'
#' @param type (character) The transform type.
#' @param desc (list) The transform descriptor.
#' @param handle (DataHandle) The data handle.
#' @return Throws an error because no specific method is defined.
#' @export
#' @keywords internal
forward_step.default <- function(type, desc, handle) {
  abort_lna(
    sprintf(
      "No forward_step method implemented for transform type: %s",
      type
    ),
    .subclass = "lna_error_no_method"
  )
}

#' Apply an inverse transform step.
#'
#' @param type (character) The type identifier of the transform (e.g., "mask", "pca").
#' @param desc (list) The parsed JSON descriptor for this transform step.
#' @param handle (DataHandle) The current data handle.
#'
#' @return An updated `DataHandle` object after applying the inverse step.
#' @export
invert_step <- function(type, desc, handle) {
  UseMethod("invert_step", type)
}

#' Default method for invert_step.
#'
#' @param type (character) The transform type.
#' @param desc (list) The transform descriptor.
#' @param handle (DataHandle) The data handle.
#' @return Throws an error because no specific method is defined.
#' @export
#' @keywords internal
invert_step.default <- function(type, desc, handle) {
  abort_lna(
    sprintf(
      "No invert_step method implemented for transform type: %s",
      type
    ),
    .subclass = "lna_error_no_method"
  )
}
