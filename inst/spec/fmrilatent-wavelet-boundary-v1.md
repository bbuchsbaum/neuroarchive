# fmrilatent and multiresolution LNA: package-boundary decision

Decision date: 2026-08-25. Status: accepted design input for QPP-11; no
fmrilatent source code was changed.

## Evidence inspected

fmrilatent already contains `cdf53_time_lift` in
`src/active_pencil_wavelet.cpp`. It is a compact floating-point transform using
half/quarter updates, symmetric endpoint duplication, and recursive coarse
prefixes. It is used by the spatial wavelet-active workflow when temporal
levels are requested. It is useful implementation evidence and a scientific
prototype.

It is not the LNA wire oracle:

- its state and coefficients are binary64, not the reversible integer 5/3
  domain;
- its rounding rule is floating arithmetic rather than the frozen floor rules;
- it has no versioned temporal coefficient/band bitstream;
- it does not serialize independence units, lossy support, progressive bands,
  quantizer steps, or decoder compatibility metadata;
- no dedicated temporal golden-vector tests were found for
  `cdf53_time_lift` during this review.

Companion fmrilatent Mote item `bd-01M0WJX0QVS0ZE07EF3VCN3AG4` records this
boundary. No fmrilatent path was reserved or edited.

## Ownership decision

The production integer lifting decoder, descriptor, bitstream, validation,
random-access rules, golden vectors, and compatibility lifetime belong in
neuroarchive core. The reference decoder is small, deterministic, has no model
or binary dependency, and is required to interpret an archive long after the
fitting environment disappears. Making stable decoding depend on fmrilatent
would couple archive portability to another package's release cadence and
compiled API without buying a meaningful reduction in decoder complexity.

fmrilatent remains the natural home for:

- experimental floating wavelet and latent-representation workflows;
- fitted bases, calibration utilities, and scientific comparisons;
- optional construction of externally calibrated profile assets;
- research on thresholding or wavelet choices that may later be frozen into a
  new neuroarchive profile version.

An fmrilatent-produced asset is not inherently unsafe. It is `analysis-safe`
when labels, weights, steps, and normalization are fixed or fitted on disjoint
external data and serialized with immutable provenance. A fitted workflow may
instead be safe relative to complete run/fold units. Safety is determined by
fit scope and lossy support, not package name.

## Portability and release consequences

- `analysis-safe` continues to decode without fmrilatent.
- `multiresolution-run` should also have a core decoder if QPP-12 admits it.
- fmrilatent may export compatible assets, but an archive must retain all
  numerical parameters needed for decode.
- A future optional plugin is justified only for a materially larger learned or
  binary decoder. It must be named in `required_transforms`, version-pinned, and
  fail closed when absent; that exception does not apply to integer 5/3.
- neuroarchive and fmrilatent can release independently. Parity tests may be
  added later under the companion Mote issue, but neuroarchive golden vectors
  remain authoritative for the wire profile.

## QPP-12 admission rule

Integration proceeds only if measured sparse/bit-plane wavelet storage offers a
useful progressive or rate advantage over QPP at an explicitly acceptable
run-level fidelity cost. Full integer roundtrip alone is necessary but not a
reason to ship another stable profile. If admitted, the profile stays
experimental until the same corruption, installed-package, scientific, and
documentation gates used for QPP are satisfied.
