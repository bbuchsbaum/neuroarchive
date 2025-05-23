% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_hdf5.R
\name{reduce_chunk_dims}
\alias{reduce_chunk_dims}
\title{Reduce chunk dimensions toward a byte target}
\usage{
reduce_chunk_dims(chunk, dtype_size, target_bytes)
}
\arguments{
\item{chunk}{Integer vector of current chunk dimensions.}

\item{dtype_size}{Size in bytes of the datatype being stored.}

\item{target_bytes}{Target maximum chunk size in bytes.}
}
\value{
Integer vector of reduced chunk dimensions.
}
\description{
Helper used when retrying dataset writes. Starting from an existing
chunk dimension vector, halves the first dimension until the
estimated chunk size is below `target_bytes` or the dimension would
drop below 1. Returns the adjusted chunk vector.
}
\keyword{internal}
