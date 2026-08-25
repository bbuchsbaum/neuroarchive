# Analysis-safe QPP contract, version 1

Status: normative design decision for QPP-1, amended by remediation epic
`bd-01M0XD72K8DKD9DHJBMC6453Q4`. Implementations conform only when the exact
parcel relationship and selected-read integrity requirements below pass through
public and installed-package paths.

## 1. Decision and estimand

The named Neuroarchive profile IDs are `analysis-safe` and
`multiresolution-run`; their lifecycle states remain `experimental` and
`prototype-no-go`, respectively. Their transform type IDs are, respectively,
`parcel.mean` followed by `qpp`, and `parcel.mean` followed by
`temporal.lift53`. Planned public constructors are
`lna_profile_analysis_safe()` and `lna_profile_multiresolution_run()`; split
validation is exposed by `assert_split_compatible()`.

`analysis-safe` reconstructs a fixed parcel representation with an independently
quantized value for every parcel and TR. It does not reduce the number of time
observations. Temporal prediction changes only the byte representation and MUST
reconstruct the quantized integers exactly.

For a run, input is a finite binary64 matrix `X[T,V]`: rows are increasing TRs
and columns are active voxels in canonical mask order. Given an externally fixed
parcel assignment `c(v)` and weights `w(v)`, the scientific estimand is

```
Y[t,k] = sum(w(v) * X[t,v], c(v) = k) / sum(w(v), c(v) = k)
q[t,k] = RNE(Y[t,k] / delta[k])
Yhat[t,k] = fl64(delta[k] * q[t,k])
```

`RNE` is round-to-nearest, ties-to-even. `fl64` is correctly rounded IEEE-754
binary64 multiplication. Each `delta[k]` is a finite, strictly positive
binary64 value stored by exact bit pattern. No clipping, saturation, wraparound,
or replacement of non-finite values is permitted.

Under exact arithmetic, `abs(delta[k] * q[t,k] - Y[t,k]) <= delta[k]/2`.
The archive MUST distinguish this quantizer bound from the decoded floating-point
error. The latter is at most the quantizer error plus the actual binary64
division/rounding and multiplication error. Validation receipts MUST report the
empirical maximum absolute error against the encoder input; documentation MUST
not promise an unconditional floating-point `delta/2` bound. The integers `q`
are the exact decoded object; `Yhat` is an approximation to `Y`.

## 2. Spatial contract

- The active-mask order, volume dimensions, affine, and mask checksum are
  required geometry. A full-volume label map uses `0` for excluded voxels; its
  active entries are serialized in canonical mask order.
- Source parcel labels may be arbitrary non-negative integers. The encoder MUST
  store their deterministic remapping to dense internal IDs `1..K`; ascending
  source-label order breaks ties. Decoding MUST NOT infer a new remapping.
- Every internal parcel MUST contain at least one voxel. Weights MUST be finite
  and non-negative, and every parcel weight sum MUST be positive. Empty parcels,
  negative weights, and zero denominators fail closed.
- Unweighted means are represented by weights equal to one. Counts and, when
  applicable, weights are stored so the spatial estimand is auditable.
- Units are mandatory metadata. Percent-signal-change, centering, scaling,
  nuisance regression, and other upstream normalization are provenance-bearing
  operations, not implicit codec behavior.
- Multiple runs are encoded independently. Predictors, anchors, and chunk
  indexes never cross a run boundary.
- Every named-profile QPP descriptor binds to exactly one run-matched
  `parcel.mean@1.1.0` descriptor. The binding includes descriptor name, step
  index, run ID, immutable asset identity, map checksum, mask/support checksum,
  and label path. Descriptor order and a matching parcel count are never
  sufficient evidence.
- Spatial assets may be reused across runs only by equality of the full
  canonical asset identity. Otherwise each asset has a distinct path.

Voxel-space expansion is `Xhat[t,v] = Yhat[t,c(v)]`; it is a gather operation
and does not recover within-parcel variation.

## 3. Fit-scope and leakage contract

Every fitted or lossy operation MUST declare `fit_scope`, `lossy_support`, and
`coding_dependency`. These fields describe different things and MUST NOT be
collapsed into one generic "temporal support" field.

Allowed `fit_scope` values are `fixed`, `external-calibration`, `fold`, `run`,
`session`, and `subject`. `analysis-safe` permits only `fixed` and
`external-calibration` for parcel labels, weights, quantizer steps, and any
upstream normalization. External calibration data MUST be disjoint from the
archive's eventual evaluation observations. A step selected from the encoded
subject or run, including a maximum, MAD, variance, percentile, or target bit
rate, is not `analysis-safe`; it belongs to an explicitly scoped profile.

For `analysis-safe`:

- quantization has `lossy_support = pointwise` in parcel-by-TR space;
- parcel averaging has `lossy_support = parcel-at-same-TR` in voxel space;
- integer prediction has `coding_dependency = chunk` but is lossless;
- the machine split vocabulary records `allowed_split_unit = pointwise`; for
  this profile that means one TR observation, subject to the declared fit
  scopes of all upstream operations;
- external calibration, map, mask/support, and fitted-asset identities use the
  canonical form `sha256:` followed by exactly 64 lowercase hexadecimal digits.

Whole-run centering is not automatically safe merely because a spatial basis is
fixed. For a lossy projection `P`, reconstruction after centering can be written
`Xhat[t] = P X[t] + (I-P) mu`; if `mu` includes held-out observations and
`P != I`, training reconstructions depend on those observations. Consequently,
run-wide centering, PSC scaling, or other data-fitted normalization disqualifies
the universal profile unless its parameters are external/fixed or the
implementation proves exact cancellation for the represented estimand.

Predictor choice may inspect its current chunk because the chosen representation
is exactly invertible and cannot change any `q[t,k]`. This affects compression
rate, not reconstructed observations.

`multiresolution-run` has `allowed_split_unit = run` and run-level lossy
support. A valid fold MUST be a union
of whole encoded runs (or larger declared independence units). Splitting TRs
inside a transformed run requires recompression within each fold.

## 4. Quantizer and predictor determinism

Quantizer policies are typed, not free-form strings:

- `fixed`: literal externally specified `delta` values and units;
- `external-calibration`: values plus immutable calibration provenance;
- `fold`, `run`, `session`, or `subject`: fitted policies, forbidden by
  `analysis-safe` and valid only when the declared split contract permits them.

The encoder MUST reject a value whose rounded integer is outside the versioned
signed integer domain. It MUST NOT silently increase `delta` or clip the value.

Each independently decodable chunk stores anchors and one predictor per parcel
or declared parcel group:

- `P0`: `r[t] = q[t]`;
- `P1`: `r[t] = q[t] - q[t-1]`;
- `P2`: `r[t] = q[t] - 2*q[t-1] + q[t-2]`.

All arithmetic is checked signed integer arithmetic. A predictor that overflows
is ineligible. The canonical cost function and golden bytes are frozen by the
QPP reference-oracle ticket. Equal costs select the lowest predictor ID in the
order `P0`, `P1`, `P2`. Chunk length and parcel-grouping choices are explicit
profile parameters, never decoder heuristics.

## 5. Version and wire rules

The v1 QPP wire contract uses:

- little-endian multibyte integers and IEEE-754 values;
- two's-complement signed integers and the versioned zig-zag mapping;
- least-significant-bit first packing within each byte;
- explicit lengths, shapes, bit widths, anchors, predictor IDs, and checksums;
- independently addressable run/chunk/parcel-group payloads.

`profile_id`, `codec_version`, and `decoder_version` are mandatory. A decoder
MUST support an explicit allow-list of codec versions and fail closed on every
unknown version, including an unknown minor version; it MUST NOT guess forward
compatibility. `decoder_version` records the deterministic implementation recipe
but does not override `codec_version`. Exact field layout, numeric bounds,
checksums, and golden byte fixtures are owned by the QPP layout and oracle
tickets before any implementation is labelled v1.

The frozen inner byte codec remains `1.0.0`. The remediation adds mandatory
descriptor relationships and digest datasets without changing inner bytes, so
both `parcel.mean` and `qpp` transform descriptors advance to `1.1.0`.
`qpp@1.0.0` can be decoded only as an explicitly unvalidated legacy/generic
parcel-series archive and never validates as `analysis-safe`. Standalone
already-parcellated data uses `profile_id = "parcel-series-qpp"`; it has no
voxel-space estimand and no `parcel.mean` binding.

The file-level `required_transforms` list MUST name every non-core transform
needed to decode the archive. A missing decoder, inconsistent descriptor,
checksum failure, invalid bit width, truncated payload, integer overflow, or
out-of-range index is a hard validation/read error.

## 6. Package boundary and fmrilatent

Neuroarchive owns the stable LNA profile definitions, transform descriptors,
bitstream, reference decoder behavior, validation, random-access semantics, and
compatibility policy. The `analysis-safe` decoder MUST NOT depend on fmrilatent.
This is a portability and archive-lifecycle boundary; it is not a claim that
fmrilatent is scientifically unsafe.

fmrilatent may own and version externally fitted parcel bases, calibration
tools, scientific evaluation utilities, and an experimental lifting-wavelet
prototype. An fmrilatent-produced asset is compatible with `analysis-safe` when
its labels, weights, quantizer steps, and normalization parameters are fixed or
fit on disjoint calibration data and their provenance is serialized. A fitted
fmrilatent workflow can instead be safe for a declared run/fold boundary. It is
not universally split-safe when it fits centering, a basis, thresholds, or
quantizer parameters using observations later assigned to evaluation.

Any companion fmrilatent source change requires its own Mote issue, path
reservation, tests, and release decision. Neuroarchive MUST retain sufficient
metadata to decode archives without importing fmrilatent or locating its
original fitting environment.

## 7. Non-goals and downstream gates

This decision does not claim that QPP reduces temporal dimensionality, that a
particular effective bit rate will be achieved, or that RMSE alone establishes
scientific validity. Encryption and access-pattern side channels are out of
scope.

Before either profile is stable, downstream work MUST provide:

1. an independent integer oracle and cross-language golden fixtures;
2. a complete LNA layout, chunk index, checksums, and compatibility tests;
3. exact parcel-mean and QPP kernels with adversarial overflow/tie tests;
4. range-read capability and fail-closed validation;
5. machine-checkable split compatibility and provenance;
6. scientific gates for spectra, autocorrelation, GLM, connectivity,
   single-trial estimates, and MVPA, plus decode-rate and size benchmarks;
7. an integer CDF 5/3 specification and run-boundary tests for the optional
   multiresolution profile.

The first production release MUST call the profile experimental until all of
those gates pass on untouched evaluation data.
