# Neuroarchive BIDS fMRI exchange profile 1.0

Status: experimental implementation contract for neuroarchive 0.1.1.

This profile defines how neuroarchive represents a common BIDS fMRI dataset as
a compressed, checksummed exchange dataset. It is deliberately narrower than
BIDS. BIDS discovery and inheritance belong to an adapter; neuroarchive owns
the typed compression boundary and the portable archive.

## Scope and unit of storage

Each source imaging file maps to exactly one `.lna.h5` artifact and exactly one
manifest row. Subject, session, task, acquisition, direction, run, echo, part,
space, resolution, density, description, and any additional entities remain
exact strings. A task entity identifies the acquisition task. Conditions or
trial types within a run remain content in byte-preserved `events.tsv` files.

Changing field of view or dimensions is legal between artifacts. Every record
has independent run geometry. A single LNA file containing several runs
remains a low-level `write_lna()` mode and is not the dataset profile default.

`.lna.h5` is not a recognized BIDS imaging extension. A neuroarchive exchange
tree is a BIDS-shaped envelope, not a claim of ordinary BIDS validity. The
declaration records the extension and media type under the neuroarchive
profile so downstream tools do not need to infer it from filenames.

## Ownership boundary

| Concern | Owner |
|---|---|
| Dataset discovery and BIDS entity parsing | bidser adapter |
| BIDS metadata inheritance and provenance order | bidser adapter |
| Events and derivative-pipeline discovery | bidser adapter |
| Source record schema and validation | neuroarchive |
| Grid, support, and reference identities | neuroarchive |
| LNA compression and numerical error contract | neuroarchive |
| Path planning and collision rejection | neuroarchive |
| Journal, manifest, declaration, query, update, transfer, restoration | neuroarchive |
| General BIDS validation, editing, or migration | external BIDS tooling |

The neuroarchive core must load, plan, write, query, update, transfer, and
restore source records with no bidser installation. `lna_source_records_from_bidser()`
is one optional producer of the source-record type. It requires bidser 0.5.1
or later.

## Source-record boundary

A source record contains:

- stable record identity;
- source path or URI, byte checksum, BIDS datatype, suffix, extension, scope,
  and derivative pipeline;
- exact entity values;
- fully resolved metadata;
- checksummed associations, including metadata sources and events;
- dimensions, full 4 by 4 affine, spatial units, voxel ordering, volume count,
  grid identity, and optional support identity; and
- a provenance-asserted reference identity.

Portable planning requires field-level provenance for every resolved metadata
field. Fabricated identity geometry is forbidden. Missing fMRI dimensions or
affine are hard errors.

## Spatial identity

Spatial compatibility has three levels:

1. `grid_id` is a pure digest of dimensions, affine, units, and voxel ordering.
2. `support_id` is a pure digest of the grid plus an exact logical support.
3. `reference_identifier` is a provenance assertion with kind, origin, status,
   evidence, and override history.

The first two identities are verified from stored values. Reference identity
cannot be inferred from an affine: distinct MNI templates can have identical
dimensions and voxel size. A missing reference becomes a record-unique
unknown-native assertion that matches only itself.

Shared resources require exact grid and support identity. They also require
the same asserted reference kind and identifier. Masked resources require an
explicit mask; only an all-voxel support may use the explicit full-grid mode.
All mismatches fail during planning or fit-resource preflight, before an image
artifact is written.

## Planning and two-pass execution

The execution order is:

1. inventory and preflight;
2. fit or verify shared resources;
3. write one artifact per source record;
4. validate every artifact;
5. derive and install the manifest and declaration.

Independent profiles have an empty fit pass. Shared-basis profiles first group
records by grid, support, and reference, then fit one content-addressed basis
per group. Every consumer stores the exact basis ID, checksum, filename, grid,
support, reference assertion, fit membership, member source checksums, fit
parameters, and update policy.

## Dataset authority

`lna_manifest.tsv` is the row authority. `lna_manifest.json` declares column
types, units, nullability, meanings, and schema version. `lna_dataset.json`
declares the profile, plan hash, manifest pair, extension registry, metadata
policy, fit-resource receipts, BIDS-envelope divergence, and optional
participants receipt.

The manifest contains source and LNA checksums, original extension, byte sizes,
datatype, suffix, scope, pipeline, all present entities, complete codec
profile, shared-resource provenance, associations, resolved metadata,
field-level metadata provenance, volume count, repetition time, reference
assertion, grid, support, and packaging status.

Manifest schema 1.1.0 is the current writer contract. It adds the required
`shared_basis_json` column, whose explicit no-resource value is
`{"status":"none"}`. Readers retain support for the already-produced 1.0.0
layout by inserting that value only when `shared_basis_json` is the sole
missing core column and the legacy TSV exactly matches its sidecar. Other
missing or incompatible core columns fail closed; the migration does not
weaken major-version or declared-column validation. Dataset declarations may
name either supported 1.x manifest version and must exactly agree with the
manifest they accompany.

Queries read the manifest only. Exact matching is the default. Regex and glob
matching require an explicit mode.

## Metadata and privacy

The adapter resolves BIDS inheritance before constructing a source record.
neuroarchive records the effective metadata and ordered source receipts; it
does not reimplement the BIDS inheritance algorithm.

Events TSV and JSON companions are copied byte for byte when selected by the
metadata policy. Participant data are excluded by default. Inclusion requires
an explicit column allowlist and selected participant IDs, and produces a
receipt for the subset written to the exchange archive.

## Restoration

Closed restoration writes the selected records plus their association and
metadata dependency closure. NIfTI dimensions, affine, and temporal length are
restored from per-record geometry. Values must pass the codec's declared
lossless or lossy error contract.

Partial restoration is permitted only as an explicitly marked partial output.
It must not claim to be a complete BIDS dataset. Structural validation and
numerical validation are separate gates.

## Updates

An update preview classifies the union of current and incoming records as:

- `unchanged`: source, semantics, profile, and artifact integrity agree;
- `added`: incoming record has no current identity;
- `changed`: the same record identity has a different source checksum;
- `source-missing`: current identity is absent from the incoming inventory;
- `stale`: source bytes agree but metadata, geometry, reference, path, policy,
  or codec semantics changed; or
- `corrupt`: the live artifact is missing or disagrees with manifest size or
  checksum.

`append` writes additions only and reports other differences. `refresh`
rebuilds changed, stale, and corrupt records without deleting missing sources.
`mirror` removes missing sources only when `prune = TRUE` was present in the
reviewed preview and the executor receives that preview's exact hash as
confirmation.

Execution builds a deterministic sibling generation and resumes from its
journal. The live archive is not replaced until the sibling manifest and
declaration validate. A successful swap retains the previous generation as a
recoverable sibling. Each update writes a new revision, plan hash, update
receipt, and append-only update events.

Publication also maintains an append-only parent-level transaction directory.
Each checksummed event records the reviewed update hash, base and target plan
hashes, stage and backup names, owner token, and the transition before or after
each directory rename. Recovery acquires the advisory lock before inspecting
this state and runs before reading the live root. A matching reviewed retry may
finish a validated staged generation when the old root has already been
preserved, or finish the receipt and journal when the target root is already
installed. If the stage is invalid, the valid previous generation is restored.
Missing, malformed, cross-plan, or multiply valid ambiguous combinations fail
with an update-recovery error and are not guessed through.

Lock ownership records process ID, host, acquisition time, and a random owner
token. A demonstrably live local owner cannot be overridden. A non-live or
unverifiable owner can be quarantined only after the configured stale interval
and an explicit `takeover_stale = TRUE`; the takeover and subsequent recovery
are appended to a parent-level audit journal.

The publication protocol assumes stage, live, backup, lock, and transaction
state share one filesystem and that an individual same-filesystem directory
rename is atomic. It does not claim that the two-directory swap is atomic;
the durable state machine closes that gap. Filesystems that do not provide
atomic same-filesystem rename, reliable close-to-open visibility, or durable
metadata after their normal sync policy are outside this profile. Windows and
network filesystems must establish those guarantees for their deployment;
cross-filesystem publication is unsupported.

Shared-basis update policy is explicit: append reuses a basis only when fit
membership and member source checksums are unchanged; refresh requires a refit
when either changes; repacking reuses only the exact content checksum.

## Transfer

The sender inventories every portable file with byte size and SHA-256 digest.
Operational locks, partial files, and transfer receipts are excluded. Paths
that collide after Unicode NFC normalization or case folding are rejected
before transfer.

The receiver copies with a bounded buffer. An interrupted partial is resumed
only after its complete prefix matches the sender byte for byte. A completed
file is installed only after size and checksum verification. Corrupt receiver
files may be repaired transactionally. Dataset authority files are transferred
last. The receiver cannot disable strict checksum verification. Verification
rebuilds the canonical inventory from the received bytes and reconciles it
with the manifest, declaration, plan, journal, associations, participants
receipt, and shared-basis resources; the receipt's file list is not trusted as
the definition of completeness. Its inventory hash proves internal
consistency. Callers that need sender authenticity must also provide the
trusted sender hash through `expected_inventory_hash`.

## Capability and maturity table

| Capability | Status | Boundary |
|---|---|---|
| Independent lossless per-file packaging | Implemented, experimental dataset API | One source image per LNA |
| Exact manifest query | Implemented | Does not open HDF5 artifacts |
| Closed and marked-partial restoration | Implemented | Numerical gate remains codec-specific |
| Resumable update generations | Implemented | Single coordinator with advisory lock |
| Resumable checksummed transfer | Implemented | Local transport engine with upload/download provenance labels |
| Optional bidser discovery | Implemented against bidser 0.5.1 API | Core remains bidser-free |
| Shared PCA basis fit and consumption | Implemented, advanced/experimental | Exact grid, support, reference, and source receipts required |
| General BIDS library behavior | Not a goal | Delegated to bidser and external tools |
| Image registration or template inference | Not supported | Must be supplied as provenance |
| `.lna.h5` as a recognized BIDS extension | Not supported | Divergence is declared programmatically |

## Required gates

A release candidate for this profile must pass:

- raw multi-subject/session/acquisition/multi-echo and MNI-derivative fixtures;
- geometry round trip for changed affine and changed field of view;
- exact entity query adversaries;
- lossless bit-exact and lossy tolerance oracles;
- metadata provenance and byte-preserved event checks;
- update no-op/add/change/remove/corrupt/resume cases;
- shared-basis grouping, mismatch, fit membership, and ROI reconstruction;
- truncated and corrupt transfer recovery with bounded buffers;
- bidser-absent core tests and bidser 0.5.1 adapter integration; and
- source-package tests, vignette build, and R CMD check.
