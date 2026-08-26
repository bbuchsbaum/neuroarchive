#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(hdf5r)
  library(neuroarchive)
})

bench_case <- function(bits = 3L, dims = c(256L, 256L), pattern = c("uniform", "structured"),
                       shuffle = FALSE, deflate_level = 0L) {
  pattern <- match.arg(pattern)
  n <- prod(dims)
  maxv <- bitwShiftL(1L, bits) - 1L

  set.seed(42)
  codes <- switch(pattern,
                  uniform = as.integer(sample.int(maxv + 1L, n, replace = TRUE) - 1L),
                  structured = as.integer(rep(seq_len(maxv + 1L) - 1L, length.out = n)))

  # Apply user-specified options for N-bit pipeline
  lna_options(write.quant.nbit.shuffle = isTRUE(shuffle),
              write.quant.nbit.deflate_level = as.integer(deflate_level))

  # Create a temp file and write dataset
  tf <- tempfile(fileext = ".h5")
  on.exit(unlink(tf), add = TRUE)
  f <- H5File$new(tf, mode = "w")
  root <- f[["/"]]

  # Write and time
  t_write <- system.time({
    dset <- neuroarchive:::h5_write_qbits_nbit(root, "/codes", matrix(codes, nrow = dims[1], ncol = dims[2]), bits = bits)
  })["elapsed"]

  # Inspect filters and dataset storage size
  filters <- tryCatch(neuroarchive::list_h5_filters(dset), error = function(e) character(0))
  dset_size <- tryCatch(dset$get_storage_size(), error = function(e) NA_real_)

  # Read and time
  t_read <- system.time({
    back <- dset$read()
  })["elapsed"]

  # Basic validation
  stopifnot(identical(as.integer(back), as.integer(matrix(codes, nrow = dims[1], ncol = dims[2]))))

  dset$close(); f$close_all()

  data.frame(
    bits = bits,
    pattern = pattern,
    dims = paste(dims, collapse = "x"),
    n_elems = n,
    shuffle = shuffle,
    deflate_level = deflate_level,
    write_s = as.numeric(t_write),
    read_s = as.numeric(t_read),
    dset_storage_bytes = as.numeric(dset_size),
    file_size_bytes = as.numeric(file.info(tf)$size),
    filters = paste(filters, collapse = ",")
  )
}

run_suite <- function() {
  cases <- expand.grid(bits = c(3L, 7L),
                       pattern = c("uniform", "structured"),
                       shuffle = c(FALSE, TRUE),
                       deflate_level = c(0L, 6L),
                       stringsAsFactors = FALSE)

  rows <- vector("list", nrow(cases))
  for (i in seq_len(nrow(cases))) {
    cat(sprintf("[%d/%d] bits=%d, pattern=%s, shuffle=%s, deflate=%d\n",
                i, nrow(cases), cases$bits[i], cases$pattern[i], cases$shuffle[i], cases$deflate_level[i]))
    rows[[i]] <- bench_case(bits = cases$bits[i], pattern = cases$pattern[i],
                            shuffle = cases$shuffle[i], deflate_level = cases$deflate_level[i])
  }
  do.call(rbind, rows)
}

if (sys.nframe() == 0L) {
  res <- run_suite()
  print(res)
  # Simple summary by bits/pattern
  agg <- aggregate(cbind(file_size_bytes, dset_storage_bytes, write_s, read_s) ~ bits + pattern + shuffle + deflate_level,
                   data = res, FUN = median)
  cat("\nMedian summary (by bits/pattern/shuffle/deflate):\n")
  print(agg)
}

