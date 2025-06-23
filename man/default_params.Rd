\name{default_params}
\alias{default_params}
\title{Get default parameters for a transform}
\usage{
default_params(type)
}
\arguments{
  \item{type}{Character scalar specifying the transform type (e.g., "quant", "delta", "temporal").}
}
\value{
A list of default parameter values extracted from the transform's JSON schema.
}
\description{
Retrieves the default parameter values for a specified transform type by
parsing its JSON schema. Results are cached for efficiency.
}
\keyword{internal}