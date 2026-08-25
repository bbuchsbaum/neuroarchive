# QPP calibration-data route

The committed QPP fixture corpus is synthetic and deterministic. Real fMRI used
to select fixed parcelwise steps or benchmark rate-distortion is calibration
data, not a test fixture and not evaluation data.

For each calibration run, record outside the package repository:

- dataset name, version, license, access date, and immutable source identifier;
- subject/run identifiers and an explicit statement that they are disjoint from
  every evaluation dataset;
- input units, mask checksum and order, parcel-map checksum, weights, and all
  preprocessing/normalization parameters and fit scopes;
- the extraction script commit, package versions, platform, and random seeds;
- SHA-256 digests for source manifests, derived parcel series, and receipts.

Do not commit restricted or large source data. Store it in an access-controlled
calibration location and commit only a non-identifying manifest template,
aggregate receipt, and the resulting externally calibrated profile parameters
when the source license permits redistribution. A profile generated from the
subject or run being encoded MUST be labelled with its actual `subject`, `run`,
or `fold` fit scope and MUST NOT be promoted to `analysis-safe`.

Calibration reports should include original and packed bytes, header overhead,
effective bits per parcel-TR, chunk coverage for representative range reads,
maximum and parcelwise reconstruction error, spectra, autocorrelation, GLM,
connectivity, single-trial, and MVPA comparisons. Profile selection is frozen on
calibration data before evaluation data are decoded.

## Reproducible synthetic court

`inst/benchmarks/qpp-rate-distortion.R` is the redistributable development
court. It simulates separate calibration and evaluation run IDs, freezes a
candidate using calibration results only, and then reports untouched evaluation
results. It measures the QPP payload, exact v1 index/step/settings components,
complete selected LNA file size, sequential and random-access runtime,
R-visible allocation, a working-set upper bound, and all required scientific
agreement estimands. Baselines state their byte boundary and loss contract;
voxel-space ratios remain labelled accounting identities unless `V/K` is
actually supplied.

Reproduce the committed receipt from the repository root with:

```r
devtools::load_all(quiet = TRUE)
source("inst/benchmarks/qpp-rate-distortion.R")
receipt <- qpp_run_rate_distortion_court(
  n_runs = 2L, n_time = 192L, n_parcel = 128L,
  alpha = c(0.00025, 0.001, 0.004),
  time_chunk = c(64L, 128L), parcel_chunk = c(64L, 128L),
  repetitions = 5L
)
qpp_write_rate_distortion_receipt(
  receipt, "inst/benchmarks/qpp-gates-v1.json"
)
qpp_write_rate_distortion_receipt(
  receipt$corpus_manifest,
  "inst/benchmarks/qpp-calibration-manifest-v1.json"
)
```

The 2026-08-25 synthetic court froze external-calibration
`alpha = 0.001`, 128-TR by 128-parcel tiles, and the raw v1 outer codec.
Untouched evaluation measured 7.65 payload bits and 8.08 total transform bits
per parcel-TR, or 3.96 times compression relative to serialized parcel
float32. The full receipt records timing, hardware, uncertainty, baselines, and
scientific metrics.

This is deliberately **not** a stable-profile receipt. The selected candidate
missed the predeclared maximum connectivity-edge error ceiling (0.008497
observed versus 0.008 allowed), and no restricted or representative external
fMRI corpus is committed here. Do not adjust that tolerance after inspecting
evaluation. A new profile iteration needs a newly versioned calibration plan
and another untouched evaluation corpus. The current QPP profile therefore
remains experimental even though the other scientific checks passed.
