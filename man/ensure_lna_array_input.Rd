% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_write.R
\name{ensure_lna_array_input}
\alias{ensure_lna_array_input}
\title{Ensure LNA-Compatible Array Input}
\usage{
ensure_lna_array_input(obj)
}
\arguments{
\item{obj}{Input object.}
}
\value{
A 4D array with possible `lna.was_3d` and `lna.orig_dims` attributes.
}
\description{
Detects common `neuroim2` objects and converts them to 4D
arrays expected by the LNA engine. If the object is a 3D array
(or `DenseNeuroVol`), a singleton fourth dimension is appended and
an attribute `lna.was_3d` is set to `TRUE`. Lower-dimensional arrays
(1D vectors and 2D matrices) are also supported and converted to 4D
with appropriate singleton dimensions. The original dimensions are 
preserved in an `lna.orig_dims` attribute.
}
\keyword{internal}
