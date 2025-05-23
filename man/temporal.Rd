% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dsl_verbs.R
\name{temporal}
\alias{temporal}
\title{Temporal Basis Projection DSL verb}
\usage{
temporal(data_or_pipe, kind = NULL, ...)
}
\arguments{
\item{data_or_pipe}{Data object or `lna_pipeline`.}

\item{kind}{Optional temporal basis type (e.g., "dct").}

\item{...}{Additional parameters for the temporal transform.}
}
\value{
An `lna_pipeline` object with the temporal step appended.
}
\description{
Adds a temporal basis transform step to a pipeline. If
`data_or_pipe` is not an `lna_pipeline`, a new pipeline is
created via `as_pipeline()`.
}
\details{
Parameter values are resolved by merging schema defaults for
the `temporal` transform, global `lna_options("temporal")`, and
any user-supplied arguments.
}
