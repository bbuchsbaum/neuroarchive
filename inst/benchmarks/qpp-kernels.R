# Kernel-only QPP benchmark.
#
# This measures qpp_encode_tile()/qpp_decode_tile() on an in-memory
# time-by-parcel integer tile. It deliberately excludes quantization, HDF5 I/O,
# outer compression, file indexing, and voxel gather. Report those boundaries
# separately; do not generalize this result to end-to-end LNA throughput.

benchmark_qpp_kernels <- function(n_time = 128L, n_parcel = 128L,
                                  iterations = 100L, seed = 20260825L) {
  stopifnot(n_time > 0L, n_parcel > 0L, iterations > 0L)
  set.seed(seed)
  innovations <- matrix(sample(-3:3, n_time * n_parcel, replace = TRUE),
                        nrow = n_time)
  q <- apply(innovations, 2L, cumsum)
  if (is.null(dim(q))) q <- matrix(q, ncol = 1L)

  encode_once <- function() neuroarchive:::qpp_encode_tile(q)
  bytes <- encode_once()
  decode_once <- function() neuroarchive:::qpp_decode_tile(bytes)
  stopifnot(identical(unname(decode_once()), unname(q * 1)))

  elapsed_encode <- system.time(for (i in seq_len(iterations)) encode_once())[["elapsed"]]
  elapsed_decode <- system.time(for (i in seq_len(iterations)) decode_once())[["elapsed"]]

  allocation_bytes <- function(fun) {
    path <- tempfile("qpp-rprofmem-")
    on.exit(unlink(path), add = TRUE)
    Rprofmem(path)
    on.exit(Rprofmem(NULL), add = TRUE)
    fun()
    Rprofmem(NULL)
    lines <- readLines(path, warn = FALSE)
    sizes <- suppressWarnings(as.numeric(sub(" .*", "", lines)))
    sum(sizes, na.rm = TRUE)
  }

  values <- length(q) * iterations
  elapsed <- pmax(c(elapsed_encode, elapsed_decode), .Machine$double.eps)
  data.frame(
    boundary = c("QPP1 encode in memory", "QPP1 decode in memory"),
    values = values,
    seconds = c(elapsed_encode, elapsed_decode),
    million_values_per_second = values / elapsed / 1e6,
    R_allocation_bytes_per_call = c(allocation_bytes(encode_once),
                                    allocation_bytes(decode_once)),
    encoded_bytes = length(bytes),
    stringsAsFactors = FALSE
  )
}
