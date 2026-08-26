# LNA v2 Validation Contract

This document defines the minimum validation contract for Latent
NeuroArchive files written by `neuroarchive`. It is intentionally narrower
than the full LNA design space: a file can contain optional metadata and
transform-specific groups beyond this contract, but it is not considered
valid unless the checks below pass.

The public `validate_lna()` API is the trust boundary for this contract.
Fast reader checks such as `quick_validate_lna()` may verify only a subset of
these rules before opening an archive for normal use.

## File Identity

An LNA v2 file is an HDF5 file whose root group has these required
attributes:

- `lna_spec`: must equal `LNA R v2.0`.
- `creator`: character metadata identifying the writer.
- `required_transforms`: a character vector of transform types needed to read
  the archive.

If `lna_spec` is missing or has an unsupported value, the file is not an LNA
v2 file for this package.

## Required Groups

The root group must contain these groups:

- `/transforms`: transform descriptors and transform reports.
- `/basis`: basis payloads or an empty group when no basis is used.
- `/scans`: run-level data.

Optional groups include:

- `/header`: spatial and acquisition metadata derived from the input object or
  supplied by the user.
- `/plugins`: plugin metadata.
- `/spatial`: masks, block tables, and spatial indexing metadata.

Optional groups may be absent. Required groups must be HDF5 groups, not
datasets.

## Run Layout

Run-specific data lives below `/scans/{run_id}`. Run IDs are sanitized before
writing and normally use labels such as `run-01`, `run-02`, or user supplied
labels that pass the package sanitizer.

When a transform is fitted or executed independently for one run, its
descriptor contains a scalar `run_id` owner. Readers apply that descriptor only
while reconstructing the matching `/scans/{run_id}` group. Descriptors without
`run_id` retain the legacy archive-wide meaning used by explicitly cross-run
operations.

For archives with stored data, each descriptor-referenced run payload must be
reachable through an absolute HDF5 path. The validator should not require a
specific internal run subtree beyond the paths declared by transform
descriptors, because different transforms can use different storage layouts.

## Transform Descriptors

Transform descriptors are stored in `/transforms` as scalar UTF-8 JSON
datasets. Descriptor names must end in `.json`. Diagnostic reports such as
`*_report.json` are not transform descriptors and should be skipped by
descriptor validation unless a transform explicitly references them.

Each transform descriptor must:

- parse as a single JSON object;
- contain a string `type`;
- validate against `inst/schemas/{type}.schema.json` from `neuroarchive` or a
  loaded plugin namespace that provides the schema;
- describe all payload datasets needed by that transform using `datasets`
  entries when the transform writes payloads.

Descriptor validation is schema validation plus LNA-specific reference
validation. A descriptor that parses and matches a JSON schema can still be
invalid if it points to missing or malformed HDF5 objects.

## Dataset References

When a descriptor has a `datasets` list, each entry with a non-empty `path`
must satisfy these rules:

- `path` resolves to an existing HDF5 object;
- the object is a dataset when the descriptor entry describes data, not a
  group or other container;
- if `dims` is provided, it must match the dataset dimensions after integer
  coercion;
- if `dtype` is provided, the dataset type must be compatible with the named
  LNA dtype;
- the dataset must be readable by `hdf5r`.

Missing referenced datasets, referenced groups where datasets are expected,
dimension mismatches, incompatible dtypes, and unreadable datasets are
validation failures.

The validator should report the descriptor name and referenced path for these
failures.

## Header and Spatial Metadata

Header metadata is optional at the format level but important for reconstructing
spatially aware `neuroim2` objects. When `/header` is present, values should be
readable as HDF5 datasets or root/header attributes according to the package
helpers.

Spatial metadata under `/spatial` is optional. If present, validators may check
known objects such as block tables and masks, but absence of `/spatial` alone
is not a validation failure.

## Basis Metadata

The `/basis` group is required even when empty. Basis-producing transforms can
store basis matrices, centers, scale/offset arrays, or descriptor-specific
payloads below this group or below descriptor-declared paths. Shared-basis
archives may reference external basis files through descriptor parameters and
dataset-level metadata.

Validation must distinguish internal and external basis contracts:

- internal basis paths declared in descriptors must exist and be readable;
- external basis locators may require a lighter check unless the caller opts
  into resolving external references.

## Checksum Semantics

When checksum mode is `"sha256"`, `materialise_plan()` writes a root
`lna_checksum` attribute. The current writer computes the SHA-256 digest on the
file with `lna_checksum` temporarily represented by a 64-character zero
placeholder, then updates the attribute to the final hash.

Validation with `checksum = TRUE` must reproduce that placeholder state before
hashing and compare the result to the stored hash. Validation with
`checksum = FALSE` skips only checksum verification. It must still perform
structural, descriptor, and dataset-reference checks.

Invalid checksum format, inability to reproduce the placeholder state, and
hash mismatch are checksum validation failures.

## Strict and Non-Strict Validation

`validate_lna(file, strict = TRUE)` should abort on the first validation
failure with class `lna_error_validation`. Error messages should include a
specific location or path whenever possible.

`validate_lna(file, strict = FALSE)` should collect validation failures and
return a structured report. Until a structured report is implemented, it may
return a character vector of issue messages for backward compatibility. It
should not silently pass malformed files.

## Validation Layers

The validator should be organized as independent layers:

1. Structure: file identity, root attributes, required groups.
2. Descriptors: descriptor discovery, JSON parsing, schema validation.
3. References: referenced HDF5 paths exist, are datasets when required, and are
   readable.
4. Shapes and dtypes: descriptor hints match actual datasets.
5. Checksum: root checksum verification when requested.

Layering keeps checksum fragility from masking structural failures and makes
corrupt-file tests easier to diagnose.

## Examples of Invalid Files

These cases must fail validation:

- root `lna_spec` is missing or not `LNA R v2.0`;
- required group `/basis`, `/scans`, or `/transforms` is missing;
- a descriptor JSON dataset cannot be parsed;
- a descriptor type has no discoverable schema;
- a descriptor references `/scans/run-01/data`, but that path is absent;
- a descriptor references a group where a dataset is expected;
- descriptor `dims` do not match the referenced dataset dimensions;
- a dataset referenced by a descriptor cannot be read;
- `lna_checksum` is present and does not match under `checksum = TRUE`.

These cases are not failures by themselves:

- `/plugins` is absent;
- `/spatial` is absent;
- `/basis` exists but is empty for a pipeline that does not use basis payloads;
- checksum verification is disabled with `checksum = FALSE`.
