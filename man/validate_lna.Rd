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
\examples{
validate_lna("ex.h5")
}
