# neuroarchive

Tools for reading and writing Latent NeuroArchive (LNA) files.

This package implements the core pieces of the LNA v2.0 format.  The
`write_lna()` and `read_lna()` helpers allow compressing numerical
arrays using a sequence of transforms.  The simplest transform is
`quant`, which performs integer quantisation, but Phase 2 adds
`basis`/`embed` for PCA-style dimensionality reduction.

Arrays fed to `write_lna()` should have the time dimension last
(for example `dim = c(nx, ny, nz, nt)` for a 4D volume).  Before a
transform runs, these arrays are reshaped into a matrix with rows
corresponding to time points and columns corresponding to voxels.
For instance a `10×4×1` array becomes a `10×4` matrix prior to
applying the Sparse PCA transform.

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

### PCA Compression Pipeline

```r
# Compress with basis -> embed -> quant
pca_file <- "pca_example.h5"
write_lna(
  x,
  pca_file,
  transforms = c("basis", "embed", "quant"),
  transform_params = list(basis = list(k = 5))
)

# Select coefficients using run_id globbing and read lazily
r2 <- read_lna(pca_file, run_id = "run-*", lazy = TRUE)
r2$subset(roi_mask = array(TRUE, dim = c(4,4,4)), time_idx = 1:2)
str(r2$data())
r2$close()
```

`read_lna()` supports glob patterns in `run_id` for selecting
multiple runs, and both `roi_mask` and `time_idx` can be specified at
read time or via the lazy reader.

### Parameter Merging

`transform_params` are resolved in the following order:

1. defaults from each transform's JSON schema via
   `default_params()`;
2. package level options set with `lna_options()`;
3. parameters supplied in the `transform_params` argument.

Later sources override earlier ones using a deep merge.

`write_lna()` also accepts a `block_table` data frame describing spatial
blocks to store under `/spatial/block_table` and an optional `plugins`
list saved in `/plugins`.

### Parallel Writing

`write_lna()` opens files using HDF5 mode `"w"`, which truncates the
target path if it already exists.  When writing in parallel each writer
should therefore create a unique temporary file and use
`file.rename()` to atomically move it into place after a successful
write.  This avoids collisions between concurrent writers.

See `vignettes/cookbook.Rmd` for additional examples including ROI/time
slicing with the lazy reader and scaffolding new transforms.

### Extending Temporal Basis Types

Packages can provide additional temporal basis kinds by defining functions of
the form `temporal_basis.<kind>()`. These methods should return a matrix with
`n_time` rows and `n_basis` columns. The generic `temporal_basis()` is used by
the core `temporal` transform and will dispatch to these user-defined methods
when `kind` matches `<kind>`.

### Standalone HRBF Utilities

`hrbf_generate_basis()`, `hrbf_project_matrix()` and
`hrbf_reconstruct_matrix()` allow working with analytic HRBF dictionaries
outside of the LNA pipeline.  Given a mask and HRBF parameter list these
helpers generate the basis, project dense data to coefficients and
reconstruct dense matrices.

```r
coeff <- hrbf_project_matrix(X, mask, params)
reconstructed <- hrbf_reconstruct_matrix(coeff, mask, params)
```

### Validation and Fork Safety

`validate_lna()` uses cached compiled JSON schemas. When running validation inside forked workers (e.g., with `future::plan(multicore)`), clear this cache in each worker using `lna:::schema_cache_clear()` to avoid potential fork-safety issues.

### Running Unit Tests

Unit tests require an R installation. Use the `run-tests.sh` helper script to
execute them from the project root. If `R` is unavailable, the script will
gracefully skip the tests.
