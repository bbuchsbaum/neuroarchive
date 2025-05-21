\name{write_lna}
\alias{write_lna}
\title{Write data to an LNA file}
\usage{
write_lna(x, file = NULL, transforms = character(),
          transform_params = list(), mask = NULL, header = NULL,
          plugins = NULL)
}
\arguments{
  \item{x}{Input object passed to \code{core_write}.}
  \item{file}{Path to output \code{.h5} file. If \code{NULL}, writing is
    performed in memory using the HDF5 core driver and no file is created.
    The returned list then contains \code{file = NULL}.}
  \item{transforms}{Character vector of transform types.}
  \item{transform_params}{Named list of transform parameters.}
  \item{mask}{Optional mask passed to \code{core_write}.}
  \item{header}{Optional named list of header attributes.}
  \item{plugins}{Optional named list saved to the \code{/plugins} group.}
}
\value{Invisibly returns a list with the output file and plan.}
\description{Compresses and stores a neuroimaging object using the
specified transform pipeline. Parameter values are merged from schema
defaults, \code{lna_options()} and \code{transform_params}.}
\examples{
x <- array(rnorm(64), dim = c(4,4,4,1))
write_lna(x, "ex.h5", transforms = "quant")
}
