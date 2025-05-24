#' Check transform implementation for namespace collisions
#'
#' Warns if the provided transform type name collides with
#' core LNA transforms or with names of base R packages.
#'
#' @param type Character scalar transform name.
#' @return Logical `TRUE` invisibly. Called for side effects (warnings).
#' @export
check_transform_implementation <- function(type) {
  stopifnot(is.character(type), length(type) == 1)

  core <- c("quant", "basis", "embed", "temporal", "delta")
  base_pkgs <- rownames(installed.packages(priority = "base"))

  msgs <- character()
  if (type %in% core) {
    msgs <- c(msgs, "core LNA transform")
  }
  if (type %in% base_pkgs) {
    msgs <- c(msgs, "base R package")
  }
  if (length(msgs) > 0) {
    warning(sprintf(
      "Transform type '%s' collides with %s namespace",
      type,
      paste(msgs, collapse = " and ")
    ), call. = FALSE)
  }

  invisible(TRUE)
}

#' Handle missing transform implementations
#'
#' Internal helper used by `core_read` and `lna_reader` to process cases
#' where a transform's S3 methods are unavailable. Behaviour depends on
#' the `allow_plugins` mode.
#'
#' @param missing_types Character vector of transform types lacking
#'   implementations.
#' @param allow_plugins One of "installed", "none", or "prompt".
#' @param location Optional string used in error conditions.
#' @keywords internal
handle_missing_methods <- function(missing_types, allow_plugins, location = NULL) {
  if (length(missing_types) == 0) return(character())

  msg <- paste0(
    "Missing invert_step implementation for transform(s): ",
    paste(unique(missing_types), collapse = ", ")
  )

  if (identical(allow_plugins, "none")) {
    abort_lna(msg, .subclass = "lna_error_no_method", location = location)
  } else if (identical(allow_plugins, "prompt") && rlang::is_interactive()) {
    response <- tolower(trimws(readline(paste0(msg, " Continue anyway? [y/N]: "))))
    if (!response %in% c("y", "yes")) {
      abort_lna(msg, .subclass = "lna_error_no_method", location = location)
    }
    warning(msg, call. = FALSE)
  } else {
    warning(msg, call. = FALSE)
  }

  invisible(missing_types)
}

#' Run a single transform step (forward or inverse)
#'
#' @param direction Character, either "forward" or "invert".
#' @param type Character, the transform type name.
#' @param desc List, the transform descriptor.
#' @param handle DataHandle, the current data state.
#' @param step_idx Integer, the index of this step in the sequence.
#' @return Updated DataHandle.
#' @keywords internal
run_transform_step <- function(direction, type, desc, handle, step_idx) {
  
  stopifnot(is.character(direction), length(direction) == 1)
  stopifnot(is.character(type), length(type) == 1)
  stopifnot(is.list(desc))
  stopifnot(inherits(handle, "DataHandle"))
  stopifnot(is.numeric(step_idx), length(step_idx) == 1)

  fun_name <- if (identical(direction, "forward")) "forward_step" else "invert_step"
  
  method_specific_fun <- getS3method(fun_name, type, optional = TRUE)
  
  if (is.null(method_specific_fun)) {
    msg <- sprintf(
      "No S3 method '%s.%s' found for transform type '%s' during %s step.",
      fun_name, type, type, direction
    )
    opposite_fun_name <- if (identical(direction, "forward")) "invert_step" else "forward_step"
    opposite_method <- getS3method(opposite_fun_name, type, optional = TRUE)
    if (!is.null(opposite_method)) {
      msg <- paste0(msg, sprintf(" However, an '%s' method does exist.", opposite_fun_name))
    }
    abort_lna(msg, .subclass = "lna_error_no_method",
              location = sprintf("%s:%s", fun_name, type))
  } else if (!is.function(method_specific_fun)) {
    abort_lna(
      sprintf("S3 method %s.%s found but is not a function. Object class: %s",
              fun_name, type, class(method_specific_fun)[1]),
      .subclass = "lna_error_internal",
      location = sprintf("run_transform_step:%s:%s", direction, type)
    )
  }

  res_handle <- tryCatch({
    method_specific_fun(type = type, desc = desc, handle = handle)
  }, error = function(e) {
    detailed_msg <- sprintf("Error in %s for transform '%s' (step %d): %s", 
                            fun_name, type, step_idx, conditionMessage(e))
    abort_lna(detailed_msg, 
              .subclass = "lna_error_transform_step", 
              parent_error = e, 
              location = sprintf("%s:%s:step%d", fun_name, type, step_idx))
  })
  
  if (!inherits(res_handle, "DataHandle")) {
    abort_lna(
      sprintf("Transform step %s for type '%s' did not return a DataHandle object. Returned class: %s", 
              fun_name, type, class(res_handle)[1]),
      .subclass = "lna_error_transform_step_return",
      location = sprintf("%s:%s:return_check", fun_name, type)
    )
  }
  
  res_handle
}
