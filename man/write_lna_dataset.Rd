\name{write_lna_dataset}
\alias{write_lna_dataset}
\title{Write a transactional BIDS-shaped LNA exchange dataset}
\usage{
write_lna_dataset(
  records,
  output_root = NULL,
  source_root = ".",
  codec_profile = list(profile_id = "lna-independent-lossless-v1", loss_mode =
    "lossless", transforms = character(), transform_params = list()),
  fit_resources = list(),
  metadata_policy = lna_metadata_policy(portable = TRUE),
  overwrite = c("error", "rebuild"),
  hooks = NULL,
  owner_id = NULL,
  recover_truncated = c("error", "truncate"),
  plan_mismatch = c("error", "restart"),
  stale_after = 300,
  takeover_stale = FALSE,
  physical_root = NULL
)
}
\arguments{
\item{records}{Source records or an existing validated
\code{lna_dataset_plan}.}

\item{output_root}{Destination root, omitted when \code{records} is a plan.}

\item{source_root}{Root used to resolve source and companion paths.}

\item{codec_profile}{Codec profile forwarded to \code{plan_lna_dataset()}.}

\item{fit_resources}{Shared fit resources forwarded to the planner.}

\item{metadata_policy}{Portable metadata and privacy policy. The default
requires provenance for every resolved metadata field and excludes participant
metadata unless explicitly opted in.}

\item{overwrite}{\code{"error"} protects unreceipted destinations;
\code{"rebuild"} authorizes transactional replacement.}

\item{hooks}{Optional named integration or fault-injection hooks for source
loading, writing, validation, installation, disk checks, companion copying,
and manifest staging.}

\item{owner_id}{Optional journal writer identity.}

\item{recover_truncated, plan_mismatch, stale_after, takeover_stale}{Journal
recovery and lock policies passed to \code{open_lna_journal()}.}

\item{physical_root}{Optional physical transaction root used while building a
replacement generation. Serialized paths continue to name \code{output_root}.}
}
\value{
An \code{lna_dataset_write_result} containing the verified plan and final
manifest, written and reused record IDs, and final journal resume state.
}
\description{
Writes one LNA per source record through a unique temporary path, validates it,
atomically installs it, validates the installed artifact, and journals the
receipt. Dataset metadata and the manifest are finalized only after all record
outputs pass. Resumption is driven by the journal and output checksums.
}
\details{
The default source loader uses \code{neuroim2::read_vec()}; no bidser API is
required. Declared path-based companions are copied byte-for-byte with checksum
verification. Existing artifacts are never truncated. Explicit rebuild mode
preserves the previous file until the replacement has installed and validated.
The declaration and manifest JSON sidecar are installed before the manifest TSV,
which remains the last dataset authority artifact.
Lossy codec profiles must include a predeclared \code{lna_error_contract()}.
}
