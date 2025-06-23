\name{schema_cache_clear}
\alias{schema_cache_clear}
\title{Schema cache environment}
\usage{
schema_cache_clear()
}
\description{
\code{.schema_cache} is an internal environment that stores compiled JSON
schemas used during validation.  It can be emptied with
\code{schema_cache_clear()} when needed, typically for testing or to
avoid stale cached objects.
}
\value{Invisibly returns \code{NULL}.}
\keyword{internal}
