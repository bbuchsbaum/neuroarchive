\name{quant.default}
\alias{quant.default}
\title{Default method for quant DSL verb}
\usage{
quant.default(data_or_pipe, bits = NULL, method = "range", center = TRUE,
      scale_scope = "global", allow_clip = FALSE, ...)
}
\description{
This is an alias for \code{\link{quant}} that provides the default method 
for the quant DSL verb. It is exported to support S3 method dispatch.
}
\seealso{\code{\link{quant}}}
\keyword{internal}