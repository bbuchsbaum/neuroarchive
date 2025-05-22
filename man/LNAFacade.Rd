\name{LNAFacade}
\alias{LNAFacade}
\title{Convenience interface for LNA}
\description{A lightweight R6 class providing \code{compress()} and \code{open()} methods that call \code{\link{write_lna}} and \code{\link{read_lna}}.}
\section{Methods}{\itemize{\item \code{compress(x, file, ...)} compress data using \code{write_lna}.\item \code{open(file, ...)} read an LNA file via \code{read_lna}.}}
\examples{
f <- LNAFacade$new()
f$compress(array(rnorm(8), dim=c(2,2,2)), "demo.h5", transforms = "quant")
f$open("demo.h5")
}
