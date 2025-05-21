# neuroarchive

Tools for reading and writing Latent NeuroArchive (LNA) files.

This package implements the core pieces of the LNA v2.0 format.  The
`write_lna()` and `read_lna()` helpers allow compressing numerical
arrays using a sequence of transforms.  The first available transform is
`quant`, which performs simple integer quantisation.

```r
library(neuroarchive)

# Example data
x <- array(rnorm(64), dim = c(4, 4, 4, 1))

# Write data with the quant transform
res <- write_lna(x,
                 file = "example.lna.h5",
                 transforms = "quant")

# Validate the file
validate_lna("example.lna.h5")

# Lazy reading using lna_reader
r <- read_lna("example.lna.h5", lazy = TRUE)
str(r$data())
r$close()
```

### Parameter Merging

`transform_params` are resolved in the following order:

1. defaults from each transform's JSON schema via
   `default_params()`;
2. package level options set with `lna_options()`;
3. parameters supplied in the `transform_params` argument.

Later sources override earlier ones using a deep merge.

### Parallel Writing

`write_lna()` opens files using HDF5 mode `"w"`, which truncates the
target path if it already exists.  When writing in parallel each writer
should therefore create a unique temporary file and use
`file.rename()` to atomically move it into place after a successful
write.  This avoids collisions between concurrent writers.

See `vignettes/cookbook.Rmd` for additional examples including ROI/time
slicing with the lazy reader and scaffolding new transforms.
