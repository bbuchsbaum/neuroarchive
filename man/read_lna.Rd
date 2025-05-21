\name{read_lna}
\alias{read_lna}
\title{Read data from an LNA file}
\usage{
read_lna(file, allow_plugins = c("warn", "off", "on"),
         validate = FALSE,
         output_dtype = c("float32", "float64", "float16"),
         lazy = FALSE)
}
\arguments{
  \item{file}{Path to an LNA file.}
  \item{allow_plugins}{Forwarded to \code{core_read}.}
  \item{validate}{Logical flag for validation.}
  \item{output_dtype}{Desired output data type.}
  \item{lazy}{If \code{TRUE}, returns an \code{lna_reader} for lazy access.}
}
\value{A \code{DataHandle} or \code{lna_reader}.}
\description{Loads data from an LNA file. With \code{lazy = TRUE} the file
remains open and data are reconstructed on demand.}
\examples{
read_lna("ex.h5", lazy = TRUE)
}
