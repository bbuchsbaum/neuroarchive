% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dsl_registry.R
\name{register_lna_verb}
\alias{register_lna_verb}
\title{Register LNA Verb}
\usage{
register_lna_verb(
  verb_name = NULL,
  lna_transform_type,
  default_slug = TRUE,
  force = FALSE
)
}
\arguments{
\item{verb_name}{A symbol or single string naming the verb.  If
`NULL` and `default_slug` is `TRUE`, the name is derived from
`lna_transform_type` by replacing non-alphanumeric characters with
underscores.}

\item{lna_transform_type}{Character string identifying the LNA
transform.}

\item{default_slug}{Logical flag controlling whether a missing
`verb_name` is slugged from `lna_transform_type`.}

\item{force}{Overwrite an existing registration with the same name.}
}
\value{
Invisibly returns a list with the registered `name` and
  `type`.
}
\description{
Register a DSL verb that maps to an LNA transform type.
}
\keyword{internal}
