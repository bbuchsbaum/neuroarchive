% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transform_builder.R
\name{TransformBuilder}
\alias{TransformBuilder}
\title{Transform Builder}
\description{
R6 class that provides a builder pattern for transform implementations,
reducing boilerplate code and ensuring consistency.
}
\keyword{internal}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-TransformBuilder-new}{\code{TransformBuilder$new()}}
\item \href{#method-TransformBuilder-add_dataset}{\code{TransformBuilder$add_dataset()}}
\item \href{#method-TransformBuilder-add_standard_dataset}{\code{TransformBuilder$add_standard_dataset()}}
\item \href{#method-TransformBuilder-add_report}{\code{TransformBuilder$add_report()}}
\item \href{#method-TransformBuilder-set_version}{\code{TransformBuilder$set_version()}}
\item \href{#method-TransformBuilder-set_io}{\code{TransformBuilder$set_io()}}
\item \href{#method-TransformBuilder-build}{\code{TransformBuilder$build()}}
\item \href{#method-TransformBuilder-get_base_name}{\code{TransformBuilder$get_base_name()}}
\item \href{#method-TransformBuilder-get_run_id}{\code{TransformBuilder$get_run_id()}}
\item \href{#method-TransformBuilder-params}{\code{TransformBuilder$params()}}
\item \href{#method-TransformBuilder-clone}{\code{TransformBuilder$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TransformBuilder-new"></a>}}
\if{latex}{\out{\hypertarget{method-TransformBuilder-new}{}}}
\subsection{Method \code{new()}}{
Create a new TransformBuilder
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TransformBuilder$new(handle, type, desc)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{handle}}{DataHandle object}

\item{\code{type}}{Transform type string}

\item{\code{desc}}{Transform descriptor}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TransformBuilder-add_dataset"></a>}}
\if{latex}{\out{\hypertarget{method-TransformBuilder-add_dataset}{}}}
\subsection{Method \code{add_dataset()}}{
Add a dataset to the transform
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TransformBuilder$add_dataset(path, data, role, dtype = NA_character_)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{path}}{Dataset path}

\item{\code{data}}{Data to store (or NULL if written directly to HDF5)}

\item{\code{role}}{Dataset role}

\item{\code{dtype}}{Data type string}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Self for chaining
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TransformBuilder-add_standard_dataset"></a>}}
\if{latex}{\out{\hypertarget{method-TransformBuilder-add_standard_dataset}{}}}
\subsection{Method \code{add_standard_dataset()}}{
Add a standard dataset with automatic path generation
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TransformBuilder$add_standard_dataset(
  dataset_type,
  data,
  role = NULL,
  dtype = NA_character_
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{dataset_type}}{Type of dataset (e.g., "quantized", "coefficients")}

\item{\code{data}}{Data to store}

\item{\code{role}}{Dataset role (defaults to dataset_type)}

\item{\code{dtype}}{Data type string}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Self for chaining
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TransformBuilder-add_report"></a>}}
\if{latex}{\out{\hypertarget{method-TransformBuilder-add_report}{}}}
\subsection{Method \code{add_report()}}{
Add a report dataset
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TransformBuilder$add_report(report_data, report_type = "report")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{report_data}}{Report data (will be JSON encoded and gzipped)}

\item{\code{report_type}}{Type of report (default: "report")}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Self for chaining
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TransformBuilder-set_version"></a>}}
\if{latex}{\out{\hypertarget{method-TransformBuilder-set_version}{}}}
\subsection{Method \code{set_version()}}{
Set transform version
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TransformBuilder$set_version(version = "1.0")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{version}}{Version string}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Self for chaining
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TransformBuilder-set_io"></a>}}
\if{latex}{\out{\hypertarget{method-TransformBuilder-set_io}{}}}
\subsection{Method \code{set_io()}}{
Set input/output specifications
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TransformBuilder$set_io(inputs = NULL, outputs = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{inputs}}{Character vector of input keys}

\item{\code{outputs}}{Character vector of output keys}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Self for chaining
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TransformBuilder-build"></a>}}
\if{latex}{\out{\hypertarget{method-TransformBuilder-build}{}}}
\subsection{Method \code{build()}}{
Build the transform and update handle
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TransformBuilder$build(stash_keys = character(), stash_values = list())}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{stash_keys}}{Keys to remove from stash}

\item{\code{stash_values}}{Named list of values to add to stash}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Updated handle
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TransformBuilder-get_base_name"></a>}}
\if{latex}{\out{\hypertarget{method-TransformBuilder-get_base_name}{}}}
\subsection{Method \code{get_base_name()}}{
Get the base name for this transform
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TransformBuilder$get_base_name()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
Character string
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TransformBuilder-get_run_id"></a>}}
\if{latex}{\out{\hypertarget{method-TransformBuilder-get_run_id}{}}}
\subsection{Method \code{get_run_id()}}{
Get the run ID
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TransformBuilder$get_run_id()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
Character string
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TransformBuilder-params"></a>}}
\if{latex}{\out{\hypertarget{method-TransformBuilder-params}{}}}
\subsection{Method \code{params()}}{
Get or set transform parameters
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TransformBuilder$params(params = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{params}}{Named list of parameters (if setting)}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Current parameters (if getting)
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TransformBuilder-clone"></a>}}
\if{latex}{\out{\hypertarget{method-TransformBuilder-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TransformBuilder$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
