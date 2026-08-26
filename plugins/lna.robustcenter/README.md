# lna.robustcenter

An example neuroarchive plugin package that ships a schema and S3 methods for a
simple `robust_center` transform (median/mean centering per column).

- Schema file: `inst/schemas/robust_center.schema.json`
- S3 methods: `forward_step.robust_center`, `invert_step.robust_center`

## Install (from this repository)

```r
# from the repository root
# install.packages("devtools")
devtools::install("plugins/lna.robustcenter")
```

## Use in a pipeline

```r
library(neuroarchive)
library(lna.robustcenter)

x <- array(rnorm(4*4*4*5), dim = c(4,4,4,5))
res <- write_lna(
  x,
  tempfile(fileext = ".h5"),
  transforms = "robust_center",
  transform_params = list(robust_center = list(method = "median"))
)

# With the plugin loaded, strict validation finds the shipped schema
validate_lna(res$file)
```

This pattern makes custom transforms portable and strictly valid across machines.
