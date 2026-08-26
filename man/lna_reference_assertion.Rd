\name{lna_reference_assertion}
\alias{lna_reference_assertion}
\title{Create a provenance-asserted reference identity}
\usage{
lna_reference_assertion(
  record_id,
  identifier = NULL,
  kind = NULL,
  origin = "unknown",
  evidence = list(),
  status = NULL,
  receipt = NULL,
  history = list()
)
}
\arguments{
\item{record_id}{Record identifier that owns the assertion.}

\item{identifier}{Asserted reference identifier. Missing information becomes \code{unknown-native:<record_id>}.}

\item{kind}{One of \code{"template"}, \code{"named-reference"}, or \code{"unknown-native"}.}

\item{origin}{One of \code{"bids-space-entity"}, \code{"spatial-reference-sidecar"}, \code{"user-receipt"}, or \code{"unknown"}.}

\item{evidence}{Named JSON-safe provenance evidence.}

\item{status}{Assertion status. Usually inferred from \code{origin}.}

\item{receipt}{Explicit user-override receipt.}

\item{history}{Prior assertion snapshots. Prefer \code{override_lna_reference()} for overrides.}
}
\value{
An object of class \code{lna_reference_assertion}.
}
\description{
Reference identity records the claimed spatial meaning of an image separately
from computed grid and support hashes. It never infers a template identity from
dimensions or an affine.
}
\details{
A BIDS-space origin requires matching \code{evidence$space_entity}; a sidecar
origin requires matching \code{evidence$SpatialReference}. User overrides
require a receipt and prior assertion history.
}
