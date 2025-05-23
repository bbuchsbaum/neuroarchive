% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dsl_templates.R
\docType{data}
\name{lna_template_registry_env}
\alias{lna_template_registry_env}
\title{DSL Template Registry}
\format{
An object of class \code{environment} of length 1.
}
\usage{
lna_template_registry_env
}
\arguments{
\item{template_name}{A single string naming the template.}

\item{template_function}{A function taking `(pipeline_obj, ...)` and
returning a modified `lna_pipeline`.}

\item{force}{Overwrite an existing registration.}
}
\value{
Invisibly returns a list with the registered name.
}
\description{
Provides functions to register and apply pipeline templates.
}
\keyword{datasets}
