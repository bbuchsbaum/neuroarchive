% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_transform.R
\name{handle_missing_methods}
\alias{handle_missing_methods}
\title{Handle missing transform implementations}
\usage{
handle_missing_methods(missing_types, allow_plugins, location = NULL)
}
\arguments{
\item{missing_types}{Character vector of transform types lacking
implementations.}

\item{allow_plugins}{One of "installed", "none", or "prompt".}

\item{location}{Optional string used in error conditions.}
}
\description{
Internal helper used by `core_read` and `lna_reader` to process cases
where a transform's S3 methods are unavailable. Behaviour depends on
the `allow_plugins` mode.
}
\keyword{internal}
