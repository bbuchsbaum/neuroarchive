\name{runtime_validate_step}
\alias{runtime_validate_step}
\title{Runtime validation for a transform step}
\usage{
runtime_validate_step(type, desc, h5)
}
\arguments{
  \item{type}{Transform type name.}
  \item{desc}{Descriptor list parsed from JSON.}
  \item{h5}{An open \code{H5File} object.}
}
\value{Invisibly returns \code{TRUE} or throws an error on validation failure.}
\description{
Checks dataset paths referenced in a descriptor and verifies that all required parameters are present before a transform is executed.
}
\seealso{\code{\link{validate_lna}}}
\keyword{internal}
