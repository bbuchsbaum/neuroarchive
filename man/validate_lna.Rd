\name{validate_lna}
\alias{validate_lna}
\title{Validate an LNA file}
\usage{
validate_lna(file, strict = TRUE, checksum = TRUE)
}
\arguments{
  \item{file}{Path to the file to validate.}
  \item{strict}{If \code{TRUE} validation errors abort.}
  \item{checksum}{If \code{TRUE} verify the stored SHA256 checksum.}
}
\value{\code{TRUE} if validation succeeds.}
\description{Checks the file version and, when requested, verifies the
stored checksum and transform descriptors.}
\details{Validation caches compiled JSON schemas for speed. When running
in forked processes (for example with \code{future::plan(multicore)}),
these cached objects may not be fork-safe. Clear the cache inside each
worker with \code{lna:::schema_cache_clear()} before calling
\code{validate_lna()}.}
\examples{
validate_lna("ex.h5")
}
