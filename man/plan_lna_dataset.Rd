\name{plan_lna_dataset}
\alias{plan_lna_dataset}
\title{Plan a BIDS-shaped LNA exchange dataset}
\usage{
plan_lna_dataset(
  records,
  output_root,
  codec_profile = list(profile_id = "lna-independent-lossless-v1", loss_mode =
    "lossless", transforms = character(), transform_params = list()),
  fit_resources = list(),
  metadata_policy = lna_metadata_policy()
)
}
\arguments{
\item{records}{One \code{lna_source_record} or a list of source records.}

\item{output_root}{Destination dataset root. The planner does not create it.}

\item{codec_profile}{Profile ID or named list with \code{profile_id}, optional
\code{loss_mode}, \code{transforms}, \code{transform_params}, and
\code{error_contract}. Lossy profiles must declare dtype-specific maximum
absolute-error and RMSE bounds before planning succeeds.}

\item{fit_resources}{Optional shared-resource fit specifications. The fit
phase is explicit even when this list is empty.}

\item{metadata_policy}{Portable metadata, association, and privacy policy.
Pure planning permits unrecorded provenance by default.}
}
\value{
An immutable \code{lna_dataset_plan} with deterministic artifact mappings and a
\code{plan_hash}. Coercion with \code{as.data.frame()} yields an inspectable dry-run table.
}
\description{
This pure planning function accepts only neuroarchive-owned source records. It
does not parse filenames or call bidser. Every source image maps to one
\code{.lna.h5} artifact, and destination collisions fail before any filesystem
mutation.
}
\details{
Raw records retain BIDS-shaped subject, session, datatype, entity, and suffix
paths. Derivative records additionally retain their pipeline scope beneath
\code{derivatives/<pipeline>}.
Lossless profiles receive a zero-error numerical contract by default. Numerical
restoration receipts are validated separately from structural BIDS validity.
}
\section{Low-level multi-run mode}{

Dataset planning deliberately uses one source image per LNA artifact. The
existing \code{write_lna()} list-input API remains the explicit, non-default
low-level route for a multi-run-in-one-file archive.
}
