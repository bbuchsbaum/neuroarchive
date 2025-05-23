% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_write.R
\name{validate_named_list}
\alias{validate_named_list}
\title{Validate optional named lists}
\usage{
validate_named_list(lst, field)
}
\arguments{
\item{lst}{List or `NULL`.}

\item{field}{Field name used in error messages.}
}
\value{
The validated list or an empty list if `NULL` or empty.
}
\description{
Used for the `header` and `plugins` arguments in `core_write`.
}
\keyword{internal}
