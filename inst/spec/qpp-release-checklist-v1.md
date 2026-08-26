# QPP release checklist v1

This checklist governs the named QPP profiles. It separates implementation
correctness from scientific admission and prevents a profile name from being
treated as stronger evidence than its receipts.

## Lifecycle decision

| Profile | Writer available | Lifecycle | Release decision |
| --- | ---: | --- | --- |
| `analysis-safe` | yes | experimental | Do not advertise as stable: untouched connectivity maximum missed its frozen ceiling. |
| `multiresolution-run` | no | prototype-no-go | Do not implement the production transform: rate margin and untouched GLM-t gates failed. |

`analysis-safe` means that every lossy amplitude operation is pointwise when
all fitted assets have fixed or disjoint-external provenance. It does not mean
that whole-run normalization, subject-fitted parcellation, or other hidden
fitted operations become safe.

## Required evidence before any stable-profile claim

- [x] Frozen wire layout, decoder version, integer domain, byte order, bit order,
  predictor cost, and canonical zero representation.
- [x] Golden byte fixtures and independent reference-decoder parity.
- [x] Exact integer round-trip tests across P0/P1/P2, tile edges, partial reads,
  and corrupt/adversarial archives.
- [x] Machine-checkable fit scope, lossy support, coding dependency, and split
  compatibility contract.
- [x] Full validation receipts with empirical maximum error and parcelwise RMSE.
- [x] Calibration/evaluation separation with an immutable benchmark manifest.
- [ ] All frozen untouched scientific ceilings pass without retuning.
- [x] Clean installed-tarball consumer using public APIs only.
- [x] Focused QPP/profile tests and the installed-package QPP vignette render.
- [ ] Full package tests and `R CMD check` are green.
- [x] The macOS, Windows, Linux release/devel/oldrel CI matrix runs both golden
  package tests and the installed consumer.
- [ ] Hosted cross-platform jobs have produced green artifacts for this change.

The unchecked scientific item blocks stable admission. Engineering gates below
may still be completed to make the experimental workflow inspectable and
reproducible.

## Remediation implementation gate

Before epic `bd-01M0XD72K8DKD9DHJBMC6453Q4` can close:

- [x] `parcel.mean@1.1.0` and `qpp@1.1.0` are emitted for the named profile;
  the inner QPP golden bytes remain codec `1.0.0`.
- [x] Every analysis-safe QPP descriptor resolves exactly one run-matched parcel
  descriptor and verifies its asset, map, mask/support, step, and path identity.
- [x] Canonical `sha256:<64 lowercase hex>` validation is shared by public and
  lower-level write paths, schemas, contracts, readers, and validators.
- [x] Two-run public-profile archives preserve independent run ownership and
  deduplicate spatial assets only by complete immutable identity.
- [x] Normal selected reads verify selected tile bytes and parcel steps before
  decode, multiplication, or cache insertion while retaining bounded range I/O.
- [x] Legacy `qpp@1.0.0`, standalone parcel-series QPP, corrupt relationships,
  and unknown versions follow the explicit compatibility matrix in the layout
  specification.
- [x] The adversarial current-source and installed-package court passes.
- [x] Local focused, integration, build/install, consumer, vignette, package
  suite, and `R CMD check` evidence is recorded separately from hosted CI.

## Evidence locations

- Format: `spec/qpp-lna-layout-v1.md`
- Leakage contract: `spec/qpp-analysis-contract-v1.md`
- Rate-distortion receipt: `benchmarks/qpp-gates-v1.json`
- Calibration manifest: `benchmarks/qpp-calibration-manifest-v1.json`
- Wavelet reference: `spec/temporal-lift53-reference-v1.md`
- Wavelet boundary decision: `spec/fmrilatent-wavelet-boundary-v1.md`
- Wavelet receipt: `benchmarks/qpp-wavelet-prototype-v1.json`
- Installed consumer: `consumer/qpp-installed-consumer.R`
- Installed adversarial court: `consumer/qpp-adversarial-consumer.R`
- Remediation certification: `benchmarks/qpp-remediation-certification-v1.json`

## Promotion rule

Promotion requires a new versioned court and checklist. Do not edit historical
receipts or widen a failed ceiling after inspecting evaluation results. A new
court must state its estimands, calibration set, untouched evaluation set,
thresholds, platform matrix, and profile version before evaluation begins.

## Local verification record (2026-08-25)

| Gate | Result | Evidence |
| --- | --- | --- |
| Focused QPP/profile/multi-run suite | PASS | 746 assertions, 0 failures/warnings/skips across profile, benchmark, kernel, Pure-R oracle, partial-read, remediation, validation, wavelet no-go, transform, and generic multi-run files. |
| Source tarball without vignettes | PASS | `R CMD build --no-manual --no-build-vignettes`. |
| Clean tarball install | PASS | Installed `neuroarchive_0.1.1.tar.gz` into an empty temporary library. |
| Installed public courts | PASS | Both installed consumers passed: two-run write, exact provenance, split compatibility, bounded selected reads, selected integrity, Pure-R integer oracle, typed mutation matrix, and wavelet no-go. |
| Current and installed QPP vignette | PASS | Both renders completed with executable checks for finite/nonempty output, full receipt, pointwise vocabulary, selected digests, partial bytes, and named lazy runs. |
| Hosted court receipt logic | PASS locally | The exact source-directory install and court block passed in an empty library; JSON assertions verified both named consumer SHA-256 values and every required scientific/engineering status field. Workflow YAML parses with independent five-entry `QPP-public-court` and `R-CMD-check` matrices; the focused court installs hard dependencies only and does not depend on optional `Suggests`. |
| Installed package suite | PASS | `R CMD check` test subprocess: 3,482 passes, 35 expected warnings, 8 documented skips, 0 failures. |
| Source-loaded package suite | ABORT, unrelated | `devtools::test()` exited 134 after 7 `experimental_api` assertions; isolated reproduction is `basis_space_embed()` under `devtools::load_all()`. |
| Normal vignette-bearing build | ABORT, unrelated | `analysis_safe_qpp.Rmd` completed; the child process exited 134 while entering the pre-existing `cookbook.Rmd`. |
| `R CMD check` on no-vignette tarball | 0 ERROR, 2 WARNING, 2 NOTE | Both warnings report intentionally absent `inst/doc`; check-time vignette rebuilding passed. Notes are CRAN incoming/fmrilatent Suggests and the legacy C++11 declaration. |
| Hosted platform matrix | PENDING | An independent `QPP-public-court` job runs both installed courts and uploads per-platform JSON receipts on macOS, Windows, Ubuntu release/devel/oldrel. The separate package-check job cannot mask the QPP result; the dirty intended slice has not been committed or pushed. |

The source-load and normal-build aborts, plus the unrun hosted matrix, prevent a
package-wide green or published claim. They do not alter the focused or
installed QPP result, and they must not be hidden or reclassified as passing.
