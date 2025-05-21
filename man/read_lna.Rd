\name{read_lna}
\alias{read_lna}
\title{Read data from an LNA file}
\usage{
read_lna(file, run_id = NULL,
         allow_plugins = c("installed", "none", "prompt"),
         validate = FALSE,
         output_dtype = c("float32", "float64", "float16"),
         roi_mask = NULL, time_idx = NULL,
         lazy = FALSE)
}
\arguments{
  \item{file}{Path to an LNA file.}
  \item{run_id}{Character vector of run identifiers or glob patterns.}
  \item{allow_plugins}{How to handle optional-package transforms.
  One of \code{"installed"} (default), \code{"none"}, or \code{"prompt"}.
  Non-interactive sessions treat \code{"prompt"} as \code{"installed"}.
  Missing implementations are skipped with a warning unless
  \code{allow_plugins = "none"}.}
  \item{validate}{Logical flag for validation.}
  \item{output_dtype}{Desired output data type.}
  \item{roi_mask}{Optional ROI mask for subsetting.}
  \item{time_idx}{Optional integer vector of time indices.}
  \item{lazy}{If \code{TRUE}, returns an \code{lna_reader} for lazy access.}
}
\value{When \code{lazy = TRUE}, an \code{lna_reader}. Otherwise a
\code{DataHandle} for a single run or a list of them when multiple runs are
loaded.}
\description{Loads data from an LNA file. With \code{lazy = TRUE} the file
remains open and data are reconstructed on demand. Run identifiers may
be matched using glob patterns and subsetting via \code{roi_mask} and
\code{time_idx} is supported.}
\seealso{\code{\link{write_lna}}, \code{\link{validate_lna}}}
\examples{
r <- read_lna("ex.h5", run_id = "run-*", lazy = TRUE,
              roi_mask = array(TRUE, dim = c(4,4,4)), time_idx = 1:2)
r$data()
r$close()
}
