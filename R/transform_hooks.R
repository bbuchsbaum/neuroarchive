#' Transform lifecycle hooks
#'
#' Optional S3 generics that transforms can implement to validate parameters
#' or perform setup/teardown work. These are primarily used by the transform
#' registry and advanced read/write workflows.
#'
#' @param type Character scalar transform type name.
#' @param desc Transform descriptor list.
#' @param handle A `DataHandle` object.
#' @param ... Reserved for future extensions.
#'
#' @return
#' `validate_params()` returns `TRUE` when parameters are valid (or a vector of
#' issues for non-strict workflows). `prepare_transform()` returns an updated
#' `DataHandle`. `cleanup_transform()` returns `TRUE` invisibly.
#'
#' @keywords internal
NULL

#' Validate transform parameters (hook)
#'
#' @export
validate_params <- function(type, desc, handle, ...) {
  type_obj <- structure(type, class = as.character(type))
  UseMethod("validate_params", type_obj)
}

#' @export
#' @keywords internal
validate_params.default <- function(type, desc, handle, ...) {
  TRUE
}

#' Prepare a transform before execution (hook)
#'
#' @export
prepare_transform <- function(type, desc, handle, ...) {
  type_obj <- structure(type, class = as.character(type))
  UseMethod("prepare_transform", type_obj)
}

#' @export
#' @keywords internal
prepare_transform.default <- function(type, desc, handle, ...) {
  handle
}

#' Cleanup after a transform (hook)
#'
#' @export
cleanup_transform <- function(type, desc, handle, ...) {
  type_obj <- structure(type, class = as.character(type))
  UseMethod("cleanup_transform", type_obj)
}

#' @export
#' @keywords internal
cleanup_transform.default <- function(type, desc, handle, ...) {
  invisible(TRUE)
}

