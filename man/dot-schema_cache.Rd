% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_json.R
\docType{data}
\name{.schema_cache}
\alias{.schema_cache}
\title{Schema Cache Environment}
\format{
An object of class \code{environment} of length 0.
}
\usage{
.schema_cache
}
\description{
Internal environment used to store compiled JSON schema objects for
transform validation.  It is not intended for direct use but can be
emptied via [schema_cache_clear()] when needed (e.g. during unit
testing).
}
\keyword{internal}
