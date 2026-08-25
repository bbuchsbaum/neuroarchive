# LNA QPP layout and compatibility contract, version 1

Status: normative layout decision for QPP-3, amended by the QPP remediation
epic `bd-01M0XD72K8DKD9DHJBMC6453Q4`. This document specializes the LNA v2
validation contract and the analysis-safe QPP contract. Numeric fields and
offsets are zero-based unless explicitly stated otherwise.

## 1. Version identifiers and required decoders

An analysis-safe archive uses:

- root `lna_spec = "LNA R v2.0"`;
- profile ID `analysis-safe`;
- transform IDs `parcel.mean` and `qpp`;
- exact transform descriptor versions `1.1.0`;
- inner QPP byte codec version `1.0.0`;
- `required_transforms = c("parcel.mean@1.1.0", "qpp@1.1.0")`.

`required_transforms` is derived from the descriptors actually written. The
writer MUST NOT emit an always-empty value. A reader compares each identifier
with an explicit allow-list before reading tile payloads. An unsupported or
unknown exact version, including an unknown minor version, fails closed.

Adding descriptor fields, datasets, or provenance that leave existing numeric
and byte interpretation unchanged increments the minor version. Consequently,
exact parcel binding, immutable asset identity, stored tile digests, and stored
parcel-step digests advance both transform descriptors from `1.0.0` to `1.1.0`.
The inner `QPP1` frame and its golden bytes remain codec version `1.0.0`.
Changing numeric interpretation, predictor cost, inner frame bytes, index
semantics, or reconstruction increments the major version. A patch release may
clarify prose or fix an implementation to match unchanged golden bytes. A root
layout change that cannot be expressed as a new transform version requires an
`lna_spec` version change.

## 2. HDF5 object tree

For transform step `SS` and run `{run_id}`, the canonical objects are:

```text
/
├── spatial/qpp_assets/{asset_sha256}/
│   ├── parcel_labels          uint16 or uint32 [V_active]
│   ├── source_parcel_ids      uint64 [K]
│   ├── parcel_counts          uint32 [K]
│   └── parcel_weights         float64 [V_active] (optional)
├── scans/{run_id}/qpp/SS_qpp/
│   ├── payload                uint8 [N_byte]
│   ├── tile_index             packed compound [N_tile]
│   ├── parcel_steps           float64 [K]
│   ├── tile_sha256            uint8 [N_tile,32]
│   ├── parcel_step_sha256     uint8 [K,32]
│   └── quality_receipt.json   scalar UTF-8 JSON
└── transforms/
    ├── SS_parcel.mean.json
    ├── SS_qpp.json
    └── SS_qpp_report.json     scalar UTF-8 JSON (optional diagnostic)
```

There is one `payload` blob and one `tile_index` per QPP transform per run, not
one HDF5 object per tile. All descriptor paths are absolute. `parcel_labels`
lists dense internal IDs in canonical active-mask order; source-ID remapping,
counts, optional weights, geometry checksum, and mask order obey the analysis
contract. `parcel_steps` stores the exact binary64 bit patterns used by the
quantizer and is never inferred from the payload.

`asset_sha256` is the lowercase hex component of the canonical
`parcel_asset_id`. An asset identity covers map, mask/support, grid, affine,
ordering, weights, source provenance, and exact dataset contents. Assets may be
reused across runs only when that complete identity is equal; a different
identity receives a different path and is never overwritten.

The payload dataset is one-dimensional, contiguous or normally HDF5-chunked as
raw `uint8`. It MUST NOT require a non-core HDF5 filter. HDF5 chunking is an I/O
choice and does not change tile offsets, which address the logical uint8 array.

## 3. Canonical tile index

`tile_index` is an HDF5 compound dataset with a canonical 72-byte record.
Members use named fixed-width little-endian HDF5 integer types at these explicit
byte offsets:

| Member | Type | Offset | Meaning |
|---|---:|---:|---|
| `tile_id` | uint32 | 0 | Consecutive row identity, starting at zero |
| `run_ordinal` | uint32 | 4 | Run position in the descriptor's ordered run table |
| `unit_ordinal` | uint32 | 8 | Independence-unit position; equal to run ordinal for `analysis-safe` v1 |
| `time_start` | uint32 | 12 | First TR in the tile |
| `time_count` | uint32 | 16 | Number of TRs, in `1..65535` |
| `parcel_start` | uint32 | 20 | First dense parcel ID minus one |
| `parcel_count` | uint32 | 24 | Number of parcels, in `1..65535` |
| `byte_offset` | uint64 | 32 | First stored byte in `payload` |
| `byte_length` | uint64 | 40 | Stored length after outer compression |
| `uncompressed_length` | uint64 | 48 | Length of the inner `QPP1` frame |
| `inner_header_offset` | uint32 | 56 | Zero in v1 |
| `predictor_table_offset` | uint32 | 60 | Eight in v1, immediately after the fixed frame header |
| `outer_codec` | uint8 | 64 | `0=raw`, `1=gzip`, `2=lz4`, `3=zstd` |
| `flags` | uint16 | 66 | Bit 0 means SHA-256 row is present; all other v1 bits are zero |
| `digest_row` | uint32 | 68 | Row in `tile_sha256`, or `0xffffffff` when absent |

Bytes 28--31 and 65 are zero padding. The layout matches hdf5r's canonical
compound representation and is validated by member name, type, offset, and
total size rather than by a host C struct. HDF5 performs host conversion, but
the stored member types remain fixed little-endian types. Implementations using R MUST preserve
uint64 values as exact `integer64` or raw fields and MUST NOT round them through
binary64.

Rows are sorted by `(run_ordinal, unit_ordinal, time_start, parcel_start)`.
Canonical writers assign consecutive `tile_id` values and densely concatenate
stored byte ranges: the first offset is zero, each later offset is the previous
offset plus length, and the final range ends at `N_byte`. Tiles form an exact,
non-overlapping rectangular partition of `T x K`; partial final time or parcel
tiles are allowed. Duplicate, overlapping, missing, or out-of-order rectangles
are invalid.

## 4. Inner `QPP1` frame

Outer decompression yields exactly the independent-reference frame frozen by
`tests/testthat/fixtures/qpp-golden-v1.json`:

```text
offset  width  field
0       4      ASCII "QPP1"
4       2      time_count, uint16 little-endian
6       2      parcel_count, uint16 little-endian
8       ...    one parcel record followed by its packed residual bytes
```

Each parcel record is:

```text
width  field
1      predictor_id: 0=P0, 1=P1, 2=P2
1      residual_bit_width: 0..34
1      anchor_count: 0, 1, or 2 as required by predictor
1      reserved: zero
2      residual_count, uint16 little-endian
4      packed_payload_length, uint32 little-endian
4*n    signed int32 anchors, little-endian two's complement
m      zig-zag residuals, LSB-first packed, zero-padded in the final byte
```

`residual_count = time_count - anchor_count` and
`packed_payload_length = ceiling(residual_count * residual_bit_width / 8)`.
Padding bits MUST be zero. Predictor records occur in increasing dense parcel
order. The inner frame has no trailing bytes.

The descriptor freezes `integer_domain = "signed-int32-q"`,
`residual_arithmetic = "checked-exact-binary64-reference-or-wider-integer"`,
`zigzag_version = "v1"`, `byte_order = "little-endian"`,
`bit_order = "lsb-first"`, `rounding_mode = "nearest-ties-to-even"`, and
`predictor_set = c("P0", "P1", "P2")`. Production integer kernels MUST prove
equivalence to the reference oracle, including P2 values that require 34 packed
bits. The deterministic predictor cost is 32 bits per anchor plus residual
count times width; ties select P0, then P1, then P2.

`decoder_version` in the descriptor identifies the recipe that wrote the file.
It is provenance, not a request to substitute a different codec interpretation.

## 5. Outer compression and digests

Outer compression is independently selected per tile after producing the
canonical inner frame:

- `raw` stores the frame unchanged and is always supported;
- `gzip` is an RFC 1952 stream and is the portable compressed fallback;
- `lz4` and `zstd` are optional framed streams, never HDF5 plugin filters.

The canonical writer selects raw when compressed and raw stored lengths tie.
Files using LZ4 or Zstd add `outer.lz4@1.0.0` or `outer.zstd@1.0.0` to
`required_transforms`. Missing support is an error before tile allocation or
decode. A file MUST NOT silently substitute a different outer codec.

For `qpp@1.1.0`, `tile_sha256[i,]` is required and is SHA-256 over the stored
bytes identified by index row `i`, before outer decompression. Each index row
has digest flag bit 0 set and `digest_row = i - 1`; the relationship is
one-to-one. `parcel_step_sha256[k,]` is required and is SHA-256 over a domain
separator, the zero-based parcel ordinal, and the exact little-endian binary64
bit pattern stored in `parcel_steps[k]`. Selected reads fetch and verify only
the corresponding digest rows before decoding or multiplication. These digests
are independent of the optional whole-file `lna_checksum`.

## 6. Descriptor and dataset roles

The `qpp` descriptor is a JSON object with `additionalProperties` permitted for
backward-compatible provenance. These fields are required:

```json
{
  "type": "qpp",
  "version": "1.1.0",
  "profile_id": "analysis-safe",
  "decoder_version": "neuroarchive-qpp-ref-1",
  "inputs": ["parcel_series"],
  "outputs": ["quantized_parcel_series"],
  "params": {
    "matrix_order": "time-by-parcel",
    "time_count": 1,
    "parcel_count": 1,
    "time_chunk": 1,
    "parcel_chunk": 1,
    "quantizer_policy": "fixed",
    "fit_scope": "fixed",
    "lossy_support": "pointwise",
    "coding_dependency": "chunk",
    "rounding_mode": "nearest-ties-to-even",
    "integer_domain": "signed-int32-q",
    "predictor_set": ["P0", "P1", "P2"],
    "predictor_cost": "anchors32-plus-residual-width-v1",
    "predictor_tie_order": ["P0", "P1", "P2"],
    "byte_order": "little-endian",
    "bit_order": "lsb-first",
    "outer_codecs": ["raw"],
    "integrity_policy": "selected-v1",
    "parcel_binding": {
      "descriptor": "00_parcel.mean.json",
      "transform": "parcel.mean@1.1.0",
      "run_id": "run-01",
      "step_index": 0,
      "asset_id": "sha256:<64 lowercase hex>",
      "map_checksum": "sha256:<64 lowercase hex>",
      "mask_checksum": "sha256:<64 lowercase hex>",
      "labels_path": "/spatial/qpp_assets/<asset_sha256>/parcel_labels"
    }
  },
  "datasets": [
    {"role": "qpp_payload", "path": "/scans/run-01/qpp/01_qpp/payload", "dtype": "uint8", "dims": [19]},
    {"role": "qpp_tile_index", "path": "/scans/run-01/qpp/01_qpp/tile_index", "dtype": "qpp_index_v1", "dims": [1]},
    {"role": "qpp_parcel_steps", "path": "/scans/run-01/qpp/01_qpp/parcel_steps", "dtype": "float64", "dims": [1]},
    {"role": "qpp_tile_sha256", "path": "/scans/run-01/qpp/01_qpp/tile_sha256", "dtype": "uint8", "dims": [1,32]},
    {"role": "qpp_parcel_step_sha256", "path": "/scans/run-01/qpp/01_qpp/parcel_step_sha256", "dtype": "uint8", "dims": [1,32]},
    {"role": "qpp_quality_receipt", "path": "/scans/run-01/qpp/01_qpp/quality_receipt.json", "dtype": "utf8", "dims": []}
  ]
}
```

The `parcel.mean` descriptor owns the spatial label, remap, count, and optional
weight roles. Its `fit_scope` and provenance are independently validated. Every
`analysis-safe` QPP descriptor MUST contain the `parcel_binding` object above.
Validation resolves exactly one descriptor by stored descriptor name and then
requires equality of transform version, scalar run ownership, step index,
asset identity, map checksum, mask checksum, and label path. Matching only
transform type, order, or `K` is invalid; there is no first-match fallback.

The quality receipt records at least input units, exact parcel step bit patterns
or their SHA-256, maximum and parcelwise error, tile count, logical values,
payload bits, stored payload bytes, header/anchor bytes, effective bits per
parcel-TR, outer codec counts, and encoder/platform provenance.

## 7. Validation order and allocation safety

Readers and validators perform these checks before allocating a decoded tile:

1. validate LNA identity, required groups, descriptors, exact supported
   transform versions, and required-transform availability;
2. validate scalar limits `T`, `K`, chunk sizes, tile count, and configured
   maximum decoded bytes using checked multiplication;
3. validate every referenced dataset type and dimension;
4. read the bounded index and validate types, reserved values, sort order,
   rectangular coverage, dense non-overlapping byte ranges, and final blob size;
5. for a requested range, read only intersecting stored byte ranges;
6. verify each selected stored-byte digest, then outer-decompress with an explicit output
   ceiling equal to `uncompressed_length`;
7. validate inner magic, duplicate dimensions, record lengths, predictor/anchor
   consistency, zero padding, checked integer reconstruction, and no trailing
   bytes;
8. read and verify each selected parcel-step digest;
9. multiply exact `q` by the verified parcel steps.

Truncation, overflow, offset plus length wraparound, ranges beyond the blob,
overlapping stored ranges, impossible dimensions, digest mismatch, unsupported
version/codec, decompression beyond its ceiling, malformed headers, non-zero
reserved/padding bits, and decoded int32 overflow are hard errors. No failing
input may trigger allocation proportional to an untrusted product before that
product is checked against configured limits.

Public validation depth has the following fixed meaning:

- `read_lna(validate = "quick")` performs bounded version, relationship,
  extent, and index preflight; reconstruction then verifies all selected tile
  and step digests before use;
- `read_lna(validate = FALSE)` skips preflight but does not disable the current
  decoder's selected-read digest checks;
- `validate_lna(qpp = "structural")` checks the complete index, parcel binding,
  digest extents, all parcel-step digests, and receipt metadata without
  decoding payload tiles;
- `qpp = "sample"` additionally verifies and round-trips a deterministic tile
  sample;
- `qpp = "all"` verifies and round-trips every tile and is the default full
  archive court.

## 8. Minimal hand-constructed file

A minimal one-run archive for `q = matrix(-9, 1, 1)` and `delta = 1` contains
the required root groups/attributes, one active voxel labelled parcel 1, count
1, source ID 1, and the descriptors above. Its raw payload is the 19 bytes:

```text
51505031010001000005000001000100000011
```

Its sole index row is:

```text
tile_id=0 run_ordinal=0 unit_ordinal=0 time_start=0 time_count=1
parcel_start=0 parcel_count=1 byte_offset=0 byte_length=19
uncompressed_length=19 inner_header_offset=0 predictor_table_offset=8
outer_codec=0 flags=1 digest_row=0
```

The `parcel_steps` dataset is one float64 value `1`; both required digest
datasets contain one SHA-256 row. The receipt declares one logical value and 19
stored bytes, and `required_transforms` contains both exact `1.1.0` IDs. A
conforming decoder reconstructs integer `-9`, verifies the selected step, then
returns floating value `-9`.

## 9. Backward compatibility

`qpp@1.0.0` and `parcel.mean@1.0.0` remain byte-decodable for explicit legacy or
generic parcel-series recovery, but they do not satisfy the named
`analysis-safe` profile because they lack exact binding and selected-read
integrity. Normal profile validation therefore rejects an archive that claims
`profile_id = "analysis-safe"` with either legacy descriptor. A caller that
deliberately disables archive validation may decode the legacy bytes with an
`lna_warning_qpp_legacy_integrity` warning; the result cannot carry an
analysis-safe validation receipt.

Standalone already-parcellated input uses the distinct generic contract
`profile_id = "parcel-series-qpp"`. It contains no parcel binding, does not
claim voxel reconstruction, and cannot be relabelled as `analysis-safe` by a
migration. One-run and multi-run analysis-safe archives use the same `1.1.0`
contract; every run has its own binding and QPP payload. Corrupt bindings or
digests always fail before scientific values are returned.

Old readers encountering a new QPP archive must fail because
`required_transforms` names unsupported exact decoders; they must not treat the
payload as ordinary `quant`. New readers continue to validate non-QPP LNA v2
archives under their original transform contracts.

## 10. Implementation routing

`materialise_plan()` derives `required_transforms` from emitted descriptors and
records `parcel.mean@1.1.0` and `qpp@1.1.0` for the named profile. The writer
stores only terminal QPP payload, tile index, parcel steps, quality receipt, and
the parcel-map reconstruction assets; it does not materialize dense parcel
series, quantized matrices, or predictor residual matrices as intermediates.

The reader performs physical HDF5 range reads for intersecting index rows,
payload byte ranges, parcel steps, and label positions. `validate_lna()` owns
structural and corruption rejection, allocation ceilings, stored-byte digests,
exact integer round-trip evidence, and exact-version dispatch. Golden fixtures
and the Pure-R oracle define portable canonical bytes. These implementation
claims do not promote the `analysis-safe` lifecycle: scientific admission is a
separate gate recorded in `qpp-release-checklist-v1.md`.

Resource ceilings remain the limits implemented by `.qpp_validation_limits`:
at most 10,000,000 TRs, 65,535 parcels, 1,000,000 tiles, 1 TiB of declared
payload, 1 GiB of index data, and 512 MiB per decoded tile by default. A normal
range read allocates at most one selected tile plus the requested result and
reads only intersecting index rows, tile digest rows, payload ranges, selected
step values, selected step digest rows, and requested label positions.
