\name{read_lna}
\alias{read_lna}
\title{Read data from an LNA file}
\usage{
read_lna(file, allow_plugins = c("installed", "none", "prompt"),
         validate = FALSE,
         output_dtype = c("float32", "float64", "float16"),
         lazy = FALSE)
}
\arguments{
  \item{file}{Path to an LNA file.}
  \item{allow_plugins}{How to handle optional-package transforms.
  One of \code{"installed"} (default), \code{"none"}, or \code{"prompt"}.
  Non-interactive sessions treat \code{"prompt"} as \code{"installed"}.}
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
