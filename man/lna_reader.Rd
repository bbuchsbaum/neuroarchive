\name{lna_reader}
\alias{lna_reader}
\title{Class for lazy LNA reading}
\description{An R6 class returned by \code{read_lna(lazy = TRUE)}. It keeps the
HDF5 file handle open and reconstructs data on demand.}
\section{Methods}{
\itemize{
  \item \code{subset(...)} store subsetting parameters such as
        \code{roi_mask} or \code{time_idx}
  \item \code{data(...)} load data using current or supplied subset
  \item \code{close()} close the HDF5 handle
}}
\examples{
r <- read_lna("ex.h5", lazy = TRUE)
r$subset(roi_mask = array(TRUE, dim = c(4,4,4)), time_idx = 1)
data <- r$data()
r$close()
}
