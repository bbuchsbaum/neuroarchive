% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_defaults.R
\name{resolve_transform_params}
\alias{resolve_transform_params}
\title{Resolve Transform Parameters}
\usage{
resolve_transform_params(transforms, transform_params = list())
}
\arguments{
\item{transforms}{Character vector of transform types.}

\item{transform_params}{Named list of user-supplied parameters.}
}
\value{
Named list of merged parameter lists.
}
\description{
Merges transform parameters from schema defaults, package options, and
user supplied values (in that order). Performs a deep merge using
`utils::modifyList` with left-to-right precedence.
}
\keyword{internal}
