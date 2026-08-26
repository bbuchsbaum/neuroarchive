# Checked production boundary for the QPP1 integer tile codec.

.qpp_default_max_values <- 16 * 1024^2
.qpp_default_max_bytes <- 512 * 1024^2

qpp_encode_tile <- function(q,
                            max_values = .qpp_default_max_values,
                            max_bytes = .qpp_default_max_bytes) {
  q <- as.matrix(q)
  if (!is.numeric(q)) {
    abort_lna("q must be a numeric matrix",
              .subclass = "lna_error_validation",
              location = "qpp_encode_tile:q")
  }
  qpp_encode_tile_rcpp(q, as.numeric(max_values), as.numeric(max_bytes))
}

qpp_decode_tile <- function(bytes,
                            max_time = 65535L,
                            max_parcels = 65535L,
                            max_values = .qpp_default_max_values,
                            max_bytes = .qpp_default_max_bytes) {
  if (!is.raw(bytes)) {
    abort_lna("bytes must be a raw vector",
              .subclass = "lna_error_validation",
              location = "qpp_decode_tile:bytes")
  }
  qpp_decode_tile_rcpp(
    bytes,
    as.integer(max_time),
    as.integer(max_parcels),
    as.numeric(max_values),
    as.numeric(max_bytes)
  )
}
