#' Run-Length Encode a Vector
#'
#' Converts a vector into a two-column matrix with columns
#' "lengths" and "values" using base `rle`.
#'
#' @param vec Vector to encode.
#' @return A matrix with columns "lengths" and "values".
#' @keywords internal
.encode_rle <- function(vec) {
  r_obj <- rle(vec)
  matrix(c(r_obj$lengths, r_obj$values), ncol = 2,
         dimnames = list(NULL, c("lengths", "values")))
}

#' Run-Length Decode a Matrix
#'
#' Decodes a two-column matrix produced by `.encode_rle` back to a vector.
#'
#' @param mat Matrix or vector representing run-length encoded data.
#' @param expected_length Optional integer expected length of the decoded vector.
#' @param location Optional string used for error messages.
#' @return Decoded vector.
#' @keywords internal
.decode_rle <- function(mat, expected_length = NULL, location = NULL) {
  if (!is.matrix(mat)) {
    if (length(mat) == 0) {
      mat <- matrix(numeric(0), ncol = 2,
                    dimnames = list(NULL, c("lengths", "values")))
    } else if (length(mat) %% 2 == 0) {
      mat <- matrix(mat, ncol = 2,
                    dimnames = list(NULL, c("lengths", "values")))
    } else {
      abort_lna(
        "RLE matrix has incorrect number of elements to form a 2-column matrix",
        .subclass = "lna_error_runtime",
        location = paste0(location, ":rle_matrix_reshape")
      )
    }
  }
  r_obj <- structure(list(lengths = mat[, 1], values = mat[, 2]), class = "rle")
  vec <- inverse.rle(r_obj)
  if (!is.null(expected_length) && length(vec) != expected_length) {
    abort_lna(
      sprintf(
        "RLE decoded data length (%d) mismatch. Expected %d element(s).",
        length(vec), expected_length
      ),
      .subclass = "lna_error_runtime",
      location = paste0(location, ":rle_decode")
    )
  }
  vec
}
