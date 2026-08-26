\name{lna_spatial_compatibility}
\alias{lna_spatial_compatibility}
\title{Compare asserted reference and computed spatial identities}
\usage{
lna_spatial_compatibility(
  x,
  y,
  support_policy = c("same-support", "grid-only")
)
}
\arguments{
\item{x, y}{LNA source records.}

\item{support_policy}{Either \code{"same-support"} or \code{"grid-only"}.}
}
\value{
An \code{lna_spatial_compatibility} receipt with separate reference, grid, and
support results.
}
\description{
Compatibility always requires a matching provenance assertion and a matching
computed grid. The default policy additionally requires two present, equal
support identities. The receipt labels asserted and computed identity evidence
separately.
}
