#' Coerce object to a dense matrix
#'
#' Provides a simple S3 generic used by write adapters to obtain
#' a time-by-voxel matrix representation. Packages can implement
#' methods for their own classes.
#'
#' @param obj Object to coerce.
#' @return A matrix.
#' @keywords internal
#' @export
as_dense_mat <- function(obj) {
  UseMethod("as_dense_mat")
}

#' @export
#' @keywords internal
as_dense_mat.default <- function(obj) {
  if (is.matrix(obj)) {
    obj
  } else {
    as.matrix(obj)
  }
}

#' @export
#' @keywords internal
as_dense_mat.array <- function(obj) {
  d <- dim(obj)
  if (length(d) <= 2) {
    return(as.matrix(obj))
  }
  time_dim <- d[length(d)]
  vox_dim <- prod(d[-length(d)])
  mat <- matrix(as.numeric(aperm(obj, c(length(d), seq_len(length(d) - 1)))),
                nrow = time_dim, ncol = vox_dim)
  mat
}

#' Coerce object to a 4-D array
#'
#' Provides a simple S3 generic used by write adapters when
#' volumetric data is required.
#'
#' @param obj Object to coerce.
#' @return A 4-D array.
#' @keywords internal
#' @export
as_vol4d <- function(obj) {
  UseMethod("as_vol4d")
}

#' @export
#' @keywords internal
as_vol4d.default <- function(obj) {
  if (is.array(obj) && length(dim(obj)) == 4) {
    obj
  } else {
    abort_lna(
      "cannot coerce object to 4-D array",
      .subclass = "lna_error_validation",
      location = "as_vol4d.default"
    )
  }
}
