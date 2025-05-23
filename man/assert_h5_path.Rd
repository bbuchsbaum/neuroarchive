% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_hdf5.R
\name{assert_h5_path}
\alias{assert_h5_path}
\title{Assert that an HDF5 path exists}
\usage{
assert_h5_path(h5, path)
}
\arguments{
\item{h5}{An `H5File` or `H5Group` object.}

\item{path}{Character path to check.}
}
\value{
Invisibly returns `NULL` when the path exists.
}
\description{
Convenience helper to verify that a dataset or group is present
at the given path relative to `h5`. Throws an `lna_error_missing_path`
error if the path does not exist.
}
\keyword{internal}
