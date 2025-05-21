library(testthat)
library(withr)

# Ensure caches cleared between tests
teardown({
  cache_env <- get(".default_param_cache", envir = asNamespace("neuroarchive"))
  rm(list = ls(envir = cache_env), envir = cache_env)
  schema_env <- get(".schema_cache", envir = asNamespace("neuroarchive"))
  rm(list = ls(envir = schema_env), envir = schema_env)
})

skip_if_not_installed("pkgload")

# Create a temporary plugin package with a schema
local_tempdir <- withr::local_tempdir()
pkg_dir <- file.path(local_tempdir, "plugpkg")
dir.create(file.path(pkg_dir, "R"), recursive = TRUE)
dir.create(file.path(pkg_dir, "inst", "schemas"), recursive = TRUE)

writeLines("Package: plugpkg\nVersion: 0.0.1\n", file.path(pkg_dir, "DESCRIPTION"))
writeLines("S3method(forward_step,plug)\nS3method(invert_step,plug)\nexport(forward_step.plug)\nexport(invert_step.plug)", file.path(pkg_dir, "NAMESPACE"))

writeLines("forward_step.plug <- function(type, desc, handle) handle\ninvert_step.plug <- function(type, desc, handle) handle", file.path(pkg_dir, "R", "plug.R"))
writeLines('{"type":"object","properties":{"foo":{"type":"integer","default":5}}}', file.path(pkg_dir, "inst", "schemas", "plug.schema.json"))

pkgload::load_all(pkg_dir, quiet = TRUE)
on.exit(unloadNamespace("plugpkg"), add = TRUE)

defaults <- neuroarchive:::default_params("plug")

test_that("default_params finds schema in loaded plugin", {
  expect_equal(defaults, list(foo = 5L))
})
