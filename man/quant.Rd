% Manual page for the quant DSL verb
% Please edit documentation in R/dsl_verbs.R
\name{quant}
\alias{quant}
\title{Add a quantisation step to a pipeline}
\description{
Add a quantisation step to a pipeline that reduces the precision of the data
by converting it to a lower bit representation. This can significantly reduce
storage requirements while preserving most of the signal information.
}
\usage{
quant(data_or_pipe, bits = NULL, method = "range", center = TRUE,
      scale_scope = "global", allow_clip = FALSE, ...)
}
\arguments{
  \item{data_or_pipe}{Data object or \code{lna_pipeline}.}
  \item{bits}{Number of quantisation bits (1-16). If \code{NULL}, the schema default is used.}
  \item{method}{Method for deriving scale/offset: \code{"range"} or \code{"sd"}.}
  \item{center}{Logical indicating whether the data should be centred before quantisation.}
  \item{scale_scope}{Either \code{"global"} for one scale/offset or \code{"voxel"} for per-voxel parameters.}
  \item{allow_clip}{If \code{TRUE}, quantisation proceeds even when clipping exceeds \code{lna.quant.clip_abort_pct}.}
  \item{...}{Additional parameters forwarded to the transform.}
}
\value{An \code{lna_pipeline} object with the quant step appended.}
\examples{\donttest{
pipe <- as_pipeline(matrix(rnorm(10), 5, 2))
pipe <- quant(pipe, bits = 4, allow_clip = TRUE)
}}
\keyword{manip}
