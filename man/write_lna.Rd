\name{write_lna}
\alias{write_lna}
\title{Write data to an LNA file}
\usage{
write_lna(x, file = NULL, transforms = character(),
          transform_params = list(), mask = NULL, header = NULL,
          plugins = NULL, block_table = NULL, run_id = NULL,
          checksum = c("none", "sha256"))
}
\arguments{
  \item{x}{Numeric array or list of arrays. Each array should have at least
  three dimensions (\code{x}, \code{y}, \code{z}, optionally \code{time}). Lists
  denote multiple runs.}
  \item{file}{Path to output \code{.h5} file. If \code{NULL}, writing is
  performed in memory using the HDF5 core driver and no file is created. The
  returned result then contains \code{file = NULL}.}
  \item{transforms}{Character vector naming the transforms to apply in forward
  order.}
  \item{transform_params}{Named list of parameters for the transforms.}
  \item{mask}{Optional \code{LogicalNeuroVol} or 3D logical array used to subset
  voxels prior to compression.}
  \item{header}{Optional named list of header attributes stored under
  \code{/header}.}
  \item{plugins}{Optional named list saved to the \code{/plugins} group.}
  \item{block_table}{Optional data frame of spatial block coordinates stored
  under \code{/spatial/block_table}. Columns must contain 1-based voxel indices
  in masked space when a mask is provided.}
  \item{run_id}{Optional character vector of run identifiers. When \code{x} is a
  list these override \code{names(x)}.}
  \item{checksum}{One of \code{"none"} (default) or \code{"sha256"}. When
  \code{"sha256"} a checksum of the final file is stored in the
  \code{/lna_checksum} attribute.}
}
\value{Invisibly returns a list with elements \code{file}, \code{plan} and
\code{header} (class \code{"lna_write_result"}).}
\description{Compresses one or more fMRI runs using a sequence of transforms.
Parameter values are merged from the transform schemas, package options and the
user supplied \code{transform_params}. The optional \code{mask}, \code{header},
\code{plugins} and \code{block_table} are written to the corresponding HDF5
groups. Use \code{read_lna()} to reconstruct the data.}
\seealso{\code{\link{read_lna}}, \code{\link{validate_lna}}}
\examples{
x <- array(rnorm(64), dim = c(4,4,4,1))
write_lna(x, "ex.h5", transforms = "quant")
read_lna("ex.h5")
}
