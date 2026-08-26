\name{override_lna_reference}
\alias{override_lna_reference}
\title{Override a source record reference assertion with a receipt}
\usage{
override_lna_reference(
  record,
  identifier,
  receipt,
  kind = "template",
  evidence = list()
)
}
\arguments{
\item{record}{An \code{lna_source_record}.}

\item{identifier}{New asserted reference identifier.}

\item{receipt}{Named list with \code{receipt_id}, \code{asserted_by}, \code{asserted_at}, and \code{reason}.}

\item{kind}{New reference kind, normally \code{"template"}.}

\item{evidence}{Additional JSON-safe user evidence.}
}
\value{
A new \code{lna_source_record}; the input is unchanged.
}
\description{
Creates a user-receipt assertion while preserving the prior assertion and all
earlier evidence as immutable history snapshots.
}
