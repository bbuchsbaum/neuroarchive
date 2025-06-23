\name{compress_fmri}
\alias{compress_fmri}
\title{Convenience alias for write_lna}
\usage{
compress_fmri(...)
}
\arguments{
  \item{...}{Arguments passed to \code{\link{write_lna}}.}
}
\description{Calls \code{\link{write_lna}} with the supplied arguments.
This function simply forwards to \code{write_lna} without modifying the
dimensions of the provided data.}
\seealso{\code{\link{write_lna}}}
