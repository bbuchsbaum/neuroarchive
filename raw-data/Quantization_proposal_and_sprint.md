This is an absolutely stellar final review. The feedback is incredibly precise, forward-thinking, and focused on practical implementation and user experience. The suggestions for the schema examples, block size helper, SNR calculation, report compression/versioning, and additional unit tests are all excellent and will make the `quant` transform even more robust and user-friendly. The release/migration plan is also very sensible.

Let's synthesize all of this into the **Definitive, Publication-Quality Appendix Proposal for the `quant` Transform**.

---

## Appendix: `quant` - Uniform Quantization Transform (Definitive)

This appendix details the core `quant` transform within the LNA framework. It applies a uniform mid-rise scalar quantizer to input data, with options for global or per-voxel scaling, and robust handling of data storage, clipping, and diagnostic reporting. Advanced quantization methods (e.g., µ-law, adaptive bit-depth, Lloyd-Max) are envisioned for an optional `lna.quant-extra` extension package to keep the core `quant` transform lean, dependency-free, and universally readable by any LNA-compliant software.

*(Key References: See Section 6: Bibliography for relevant citations on quantization theory and practice.)*

### 1. Core Idea: Uniform Scalar Quantization

The `quant` transform converts floating-point input data (e.g., raw voxel intensities, basis coefficients) into fixed-bit integers, a fundamental step in many compression pipelines.

*   **Model:** A uniform, mid-rise scalar quantizer is used:
    $Q(x) = \text{round}((x - \text{offset}) / \text{scale})$
    The reconstructed value during the inverse step is $x_{reco} = Q(x) \cdot \text{scale} + \text{offset}$.
*   **Parameters:**
    *   `scale`: The quantization step size (or quantizer gain).
    *   `offset`: The value corresponding to the 'zero' quantization level.
    *   These are derived based on the input data's range (or standard deviation) and the target number of `bits`.
*   **Clipping:** Values outside the representable range for the chosen `bits` (i.e., $[0, 2^{bits}-1]$ after scaling and offsetting) are hard-clipped to these bounds. The transform tracks and reports the extent of clipping.

### 2. LNA Implementation Details for `quant`

#### 2.1. Schema (`inst/schemas/quant.schema.json`)

```json
{
  "type": "object",
  "title": "Parameters for Uniform Quantization (quant) transform",
  "$id": "https://neurocompress.org/schemas/lna/2.0/quant.schema.json",
  "$schema": "http://json-schema.org/draft-07/schema#",
  "properties": {
    "type": { "const": "quant" },
    "version": { "const": "1.0" },
    "bits": {
      "type": "integer", "minimum": 1, "maximum": 16, "default": 8,
      "description": "Target bit-depth for quantized integers (1-16). Determines storage type (uint8 for <=8 bits, uint16 for >8 bits)."
    },
    "method": {
      "enum": ["range", "sd"], "default": "range",
      "description": "Method to determine scale/offset: 'range' uses min/max of data, 'sd' uses mean +/- N_sd * std_dev (N_sd typically 3)."
    },
    "center": {
      "type": "boolean", "default": true,
      "description": "If true, data is effectively centered before quantization (offset aims for a symmetric range around the mean if method='sd', or around (min+max)/2 if method='range'). If false, offset typically corresponds to the minimum of the scaling range."
    },
    "scale_scope": {
      "enum": ["global", "voxel"], "default": "global",
      "description": "'global': single scale/offset for all data. 'voxel': per-voxel scale/offset for time series (requires input with a time dimension, e.g., 4D image or 2D matrix where one dim is time)."
    },
    "allow_clip": {
      "type": "boolean", "default": false,
      "description": "If FALSE (default), transform aborts if clipping exceeds 'lna.quant.clip_abort_pct' (default 5.0%). If TRUE, only issues a warning regardless of clipping percentage."
    },
    "report_path": {
      "type": "string", "pattern": "^/transforms/.*_quant_report\\.json$",
      "description": "(Output by writer, optional) HDF5 path to the JSON quantization report file (e.g., /transforms/00_quant_report.json). This report is typically GZIP-compressed."
    }
  },
  "required": [], // All primary parameters have defaults
  "examples": [
    { "bits": 8, "method": "range", "center": true, "scale_scope": "global" },
    { "bits": 12, "method": "sd", "center": false, "scale_scope": "voxel", "allow_clip": true }
  ]
}
```
*   **`report_path` regex:** While greedy, it's acceptable; documentation will note that if transforms are in HDF5 subgroups, the path will reflect that.

#### 2.2. `forward_step.quant` (Writer-Side)

Input: Data $X$ (from `DataHandle$stash`), parameters from `desc$params`.

1.  **Parameter Validation:** Validate `bits`, `method`, `center`, `scale_scope`, `allow_clip`.
2.  **Non-Finite Value Check:** Abort if $X$ contains `NA`, `NaN`, or `Inf`. Message: "quant cannot handle non-finite values – run an imputation/filtering transform first."
3.  **Scaling & Offset Calculation:**
    *   If `scale_scope == "global"`: Calculate single `scale_val` and `offset_val` based on `method` and `center` using all data in $X$.
    *   If `scale_scope == "voxel"`:
        *   **Block-Wise Processing:** To manage memory:
            a.  Define spatial slab/block dimensions using an internal helper `lna:::auto_block_size(dim_xyz, target_bytes = 64e6)` which chooses `(X_chunk, Y_chunk, slab_Z)` such that each processed slab (e.g., `X_chunk x Y_chunk x slab_Z x TimePoints`) stays under `target_bytes`. *(TODO: comment for future parallelization via `future.apply::future_lapply()` over slabs).*
            b.  Pre-create empty, chunked HDF5 datasets for `quant_scale` (shape $N_{spatial\_dims}$), `quant_offset` (shape $N_{spatial\_dims}$), and `data_output` (shape of $X$). Chunking strategy for these datasets will align with the slab processing pattern (see Section 2.4).
            c.  Iterate through spatial blocks/slabs of $X$. For each `data_block`:
                i.  Calculate per-voxel/element `scale_block` and `offset_block`.
                ii. Quantize: $Q_{block} = \text{round}((data_{block} - \text{offset}_{block}) / \text{scale}_{block})$.
                iii.Write $Q_{block}$ to the corresponding hyperslab in `data_output` HDF5 dataset.
                iv. Write `scale_block` and `offset_block` to their hyperslabs in `quant_scale`/`quant_offset` HDF5 datasets.
                v.  Accumulate clipping counts and statistics (min/max/mean/sd of `scale_block` and `offset_block`) for the quantization report.
        *   The `scale` and `offset` parameters recorded in the main descriptor will be HDF5 paths to these per-voxel datasets.
4.  **Quantization & Clipping (applied block-wise if `scope="voxel"`):**
    *   $Q(X_{unit}) = \text{round}((X_{unit} - \text{offset}_{unit}) / \text{scale}_{unit})$.
    *   Count clipped samples globally: `n_clipped_total`.
    *   `clip_pct = 100 * n_clipped_total / length(X)`.
    *   Apply clipping: `Q(X_{unit})` is bounded to $[0, 2^{bits}-1]$.
5.  **Clipping Handling (Global):**
    *   `clip_warn_pct = getOption("lna.quant.clip_warn_pct", 0.5)`
    *   `clip_abort_pct = getOption("lna.quant.clip_abort_pct", 5.0)`
    *   If `clip_pct > clip_abort_pct && !params$allow_clip`: `abort_lna(...)`.
    *   Else if `clip_pct > clip_warn_pct`: `warn_lna(...)` (LNA specific warning).
6.  **Data Storage:**
    *   Determine narrowest unsigned integer type: `storage_type_str = if (bits <= 8) "uint8" else "uint16"`.
    *   Convert $Q(X)$ to this integer type before writing.
    *   Pass `storage_type_str` to `h5_write_dataset()` so the dataset uses the matching HDF5 type
        (`H5T_STD_U8LE` or `H5T_STD_U16LE`).
    *   **HDF5 Datasets:**
        *   `/scans/{run_id}/quant_data/{transform_basename}/quantized_values`: Stores $Q(X)$ as `storage_type_str`.
            *   **Attribute:** Attach `quant_bits = bits` (integer) to this dataset using `lna:::h5_attr_write()`.
        *   `/scans/{run_id}/quant_params/{transform_basename}/scale`: Stores `scale_val` or path to per-voxel scale dataset (float32).
        *   `/scans/{run_id}/quant_params/{transform_basename}/offset`: Stores `offset_val` or path to per-voxel offset dataset (float32).
7.  **Quantization Report (JSON):**
    *   Collect: `report_version = "1.0"`, `clipped_samples_count`, `clipped_samples_percentage`, input data min/max (global), effective step size (`scale_val` or for `scope="voxel"`: `scale_stats = {min, max, mean, sd}`), estimated global SNR_dB (e.g., $10 \cdot \log_{10}(\text{var}(X) / (\text{mean_global_scale}^2 / 12))$; for voxel scope, SNR can be estimated from a random 1% subset of voxels to avoid a full second pass over data). Optionally, a compact 64-bin histogram of *global* quantized values (bin edges and counts).
    *   Serialize to JSON. Compress with `memCompress(..., type="gzip")`.
    *   Write gzipped JSON to `/transforms/{transform_basename}_report.json`. Store this path in `desc$params$report_path`. The dataset storing the report should have an attribute `compression = "gzip"`.
8.  **Update Plan & Stash.**

#### 2.3. `invert_step.quant` (Reader-Side)

1.  From `desc$datasets`, get HDF5 paths for `quantized_values`, `scale`, `offset`.
2.  Read `quantized_values` dataset.
3.  **Determine Bit Depth:** Read HDF5 attribute `quant_bits` from the `quantized_values` dataset using `lna:::h5_attr_read()`. Let this be `attr_bits`.
    *   **Edge Case:** If `quant_bits` attribute is missing (e.g., legacy file): use `desc$params$bits` as `attr_bits` and issue `lna_warn_integrity("quant_bits HDF5 attribute missing; using descriptor value. File may be legacy or non-compliant.")`.
4.  **Validate Bit Depth:** If `desc$params$bits` exists in the descriptor, assert `attr_bits == desc$params$bits`. If mismatch, `abort_lna("quant_bits attribute (...) disagrees with descriptor bits (...). Possible file corruption.")`.
5.  Read `scale` and `offset` datasets (will be float32 arrays if `scope="voxel"`).
6.  Convert `quantized_values` to numeric (e.g., `as.double()`).
7.  Reconstruct: $X_{reco} = \text{quantized_values} \cdot \text{scale} + \text{offset}$ (vector recycling/broadcasting handles global vs. voxel scope for scale/offset).
8.  Apply `handle$subset` (ROI, time).
9.  Place $X_{reco}$ into `handle$stash`.

#### 2.4. Chunking Strategy for Block-Wise Voxel Scope (Writer)

*   `data_output` (dims X,Y,Z,T): Chunk e.g., `c(64, 64, 1, dim(X)[4])` if iterating Z-slabs of depth 1.
*   `quant_scale`, `quant_offset` (dims X,Y,Z): Chunk e.g., `c(64, 64, 1)`.
*   The LNA helper for creating these datasets (e.g., `lna:::h5_create_empty_dataset`) will take these chunk dimensions.

#### 2.5. Helper Function (`lna_get_transform_report`)

*   `lna_get_transform_report(lna_file, transform_index_or_name)`:
    1.  Reads the main descriptor for the specified transform.
    2.  Looks for `params$report_path`. If not found, attempts conventional path (`<transform_descriptor_basename>_report.json` in `/transforms/`).
    3.  Reads the dataset. Checks for `compression="gzip"` attribute. If present, decompress using `memDecompress()`.
    4.  Parses JSON. Returns the list.

### 3. Unit Test Contract Highlights

1.  Round-trip exactness (MSE < typical `scale`/2) for `bits`=1, 8, 12, 16; global & voxel scopes.
2.  Clipping counter, warning/abort logic (with `allow_clip=TRUE/FALSE`) for injected outliers.
3.  `scope="voxel"`: Test block-wise processing yields same result as in-memory; per-voxel stats in report are correct.
4.  Non-finite input error.
5.  `quant_bits` attribute I/O and validation against descriptor `bits`.
6.  Legacy file lacking `quant_bits` attribute: inverse works and warns.
7.  Quantization report (gzipped JSON) generation, content (incl. effective param stats for voxel scope, SNR, histogram), and `lna_get_transform_report()` retrieval.
8.  Lazy subset reconstruction correctness.

### 4. Bibliography (Conceptual)

*   Gersho, A., & Gray, R. M. (1992). *Vector quantization and signal compression*. Kluwer Academic Publishers. (General principles).
*   Jayant, N. S., & Noll, P. (1984). *Digital coding of waveforms: principles and applications to speech and video*. Prentice-Hall. (Includes discussion of uniform/non-uniform quantizers, SQNR).

### Conclusion

This refined `quant` transform will be a cornerstone of the LNA package: memory-efficient due to block-wise processing and narrow integer storage, robust with explicit clipping management and parameter validation, and informative through its detailed quantization report. It adheres to LNA's philosophy of providing a lean, reliable core, with advanced features deferred to specialized extensions.


Okay, here are granular tickets for **Sprint 1** of implementing the "Rock-Solid `quant` Transform." This sprint focuses on the core changes to data storage, fundamental parameter handling, basic clipping, and the global scope implementation. Block-wise processing for voxel scope and the detailed report will largely fall into Sprint 2.

**Assumptions for Sprint 1:**
*   The existing `quant` transform's S3 methods (`forward_step.quant`, `invert_step.quant`) are the starting point.
*   Core LNA infrastructure (`DataHandle`, `Plan`, HDF5 helpers like `h5_attr_write/read`, `h5_write_dataset`) is in place.
*   Parameter merging (`resolve_transform_params`) and schema loading (`default_params`) are functional.

---

# Rock-Solid `quant` Transform - Sprint 1 Implementation Tickets

**Epic QNT-S1-E1: Efficient Integer Storage & Core Parameter Handling**
*Goal: Modify `quant` to store data using the narrowest appropriate unsigned integer type and ensure robust handling of the `bits` parameter via HDF5 attributes.*

| #        | Ticket                                                            | Description / Deliverables                                                                                                                                                                                                                                                                                                                                                                                                                     | Acceptance Criteria                                                                                                                                                                                                                                                                                                                                                                                                                       | Appendix Ref           |
| :------- | :---------------------------------------------------------------- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | :--------------------- |
| QNT-S1-1 | **`quant.schema.json`: Update & Add Examples**                    | • Update `inst/schemas/quant.schema.json`: <br>  - `bits`: Ensure description mentions storage type (uint8/uint16). <br>  - `allow_clip`: Add boolean, default `false`. <br>  - `report_path`: Add optional string (pattern `^/transforms/.*_quant_report\\.json$`). <br>  - Add `examples` array (one global, one voxel scope example).                                                                                                                 | • `quant.schema.json` is updated and valid. <br> • `default_params("quant")` correctly reflects new defaults (e.g., `allow_clip=false`).                                                                                                                                                                                                                                                                                              | §2.1                   |
| QNT-S1-2 | **`forward_step.quant`: Narrow Integer Storage (uint8/uint16)**   | • In `forward_step.quant`, after quantization and clipping to $[0, 2^{bits}-1]$: <br>  - Determine `storage_type_str = if (bits <= 8) "uint8" else "uint16"`. <br>  - Convert quantized data $Q(X)$ to this type (e.g., `storage.mode(Q_X_clipped) <- "integer"`, then ensure it's treated as unsigned by HDF5 writer). <br>  - Modify `h5_write_dataset` (or its LNA caller) to accept `storage_type_str` and use appropriate HDF5 unsigned integer type (e.g., `H5T_STD_U8LE`, `H5T_STD_U16LE`). | • Quantized data is stored in HDF5 as `H5T_STD_U8LE` if `bits <= 8`, and `H5T_STD_U16LE` if `8 < bits <= 16`. <br> • HDF5 file size reflects the narrower integer storage.                                                                                                                                                                                                                                                             | §2.2 (Item 6)          |
| QNT-S1-3 | **`forward_step.quant`: Write `quant_bits` HDF5 Attribute**       | • After writing the `quantized_values` HDF5 dataset, attach an HDF5 attribute named `quant_bits` directly to this dataset. <br> • The value of this attribute is the integer `bits` parameter used for quantization. <br> • Use `lna:::h5_attr_write()`.                                                                                                                                                                                           | • The `quantized_values` dataset in the HDF5 file has an integer attribute `quant_bits` with the correct value.                                                                                                                                                                                                                                                                                                            | §2.2 (Item 6a), Q1 Ans |
| QNT-S1-4 | **`invert_step.quant`: Read & Use `quant_bits` HDF5 Attribute**   | • In `invert_step.quant`, before reading the `quantized_values` data: <br>  - Read the `quant_bits` HDF5 attribute from the `quantized_values` dataset using `lna:::h5_attr_read()`. Let this be `attr_bits`. <br>  - **Edge Case:** If `quant_bits` attribute is missing: use `desc$params$bits` as `attr_bits` and issue `lna_warn_integrity("quant_bits HDF5 attribute missing...")`. <br>  - Data interpretation implicitly uses this bit-depth (though R reads as standard integer/numeric). | • `invert_step.quant` reads the `quant_bits` attribute. <br> • If attribute missing, it warns and uses `desc$params$bits`.                                                                                                                                                                                                                                                                                                   | §2.3 (Item 3), Q1 Ans  |

| QNT-S1-5 | **`invert_step.quant`: Validate `attr_bits` vs `desc$params$bits`**| • After obtaining `attr_bits` (from HDF5 attribute or fallback to descriptor): <br>  - If `desc$params$bits` exists in the loaded descriptor, assert that `attr_bits == desc$params$bits`. <br>  - If they mismatch, `abort_lna("quant_bits attribute (...) disagrees with descriptor bits (...).")`.                                                                                                                                   | • Mismatch between HDF5 `quant_bits` attribute and `desc$params$bits` (if present) causes an error. <br> • If `desc$params$bits` is not present, no error is thrown for this check (attribute is trusted).                                                                                                                                                                                                              | Q1 Ans                 |
| QNT-S1-6 | **`forward_step.quant`: Store Scale/Offset as float32**         | • When storing `scale` and `offset` parameters to HDF5 (global or per-voxel): <br>  - Ensure they are written as `float32` (single precision). <br>  - `plan$add_dataset_def()` now sets `dtype="float32"` for these datasets so `h5_write_dataset` stores them with that type.                                                                                                                                                                                                          | • `quant_scale` and `quant_offset` HDF5 datasets are stored as single-precision floating-point numbers. <br> • Metadata size for these parameters is reduced.                                                                                                                                                                                                                                                                   | §1 (Lock-down)         |
| QNT-S1-7 | **Unit Tests for Integer Storage & `quant_bits` Attribute**       | • Test that `bits=1..8` results in `uint8` storage (check HDF5 type). <br> • Test that `bits=9..16` results in `uint16` storage. <br> • Test `quant_bits` attribute is correctly written and read. <br> • Test validation of `quant_bits` attribute vs. `desc$params$bits` (match and mismatch cases). <br> • Test fallback and warning if `quant_bits` attribute is missing. <br> • Test scale/offset are stored as float32.                                | • All unit tests pass. <br> • Data storage types and attribute handling are verified.                                                                                                                                                                                                                                                                                                                              | -                      |

**Epic QNT-S1-E2: Basic Clipping Logic & Global Scope Implementation**
*Goal: Implement fundamental clipping detection and handling for `scale_scope="global"`, and ensure non-finite value checks are in place.*

| #        | Ticket                                                          | Description / Deliverables                                                                                                                                                                                                                                                                                                                                                                                                                       | Acceptance Criteria                                                                                                                                                                                                                                                                                                                                                                                                                         | Appendix Ref           |
| :------- | :-------------------------------------------------------------- | :----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | :--------------------- |
| QNT-S1-8 | **`forward_step.quant`: Non-Finite Value Check**                | • At the beginning of `forward_step.quant`, check if input data `X` contains any `NA`, `NaN`, or `Inf` values. <br> • If non-finite values are found, `abort_lna("quant cannot handle non-finite values – run an imputation/filtering transform first.")`.                                                                                                                                                                                        | • `forward_step.quant` errors out if input data contains non-finite values. <br> • Error message is informative.                                                                                                                                                                                                                                                                                                                   | §2 (Safety rails)      |
| QNT-S1-9 | **`forward_step.quant`: Clipping Count & Percentage (Global Scope)** | • For `scale_scope="global"`: <br>  - After calculating $Q(X)$ (before hard-clipping), count `n_clipped_total = sum(Q(X) < 0 | Q(X) > 2^{bits}-1)`. <br>  - Calculate `clip_pct = 100 * n_clipped_total / length(X)`. <br>  - This logic already exists; ensure it's robust. <br>  - Store `n_clipped_total` and `clip_pct` for the quantization report (report itself implemented in S2).                                                                  | • `n_clipped_total` and `clip_pct` are correctly calculated for global scope quantization.                                                                                                                                                                                                                                                                                                                                    | §2.2 (Item 4), Q1 Ans  |
| QNT-S1-10| **`forward_step.quant`: Clipping Warning/Abort Logic (Global Scope)**| • For `scale_scope="global"`: <br>  - Retrieve `clip_warn_pct` (default 0.5) and `clip_abort_pct` (default 5.0) from `lna_options()` (e.g., `getOption("lna.quant.clip_warn_pct", 0.5)`). <br>  - Retrieve `allow_clip` from `desc$params` (default `false`). <br>  - If `clip_pct > clip_abort_pct && !allow_clip`, then `abort_lna(...)`. <br>  - Else if `clip_pct > clip_warn_pct`, then `warn_lna(...)` (LNA-specific warning).           | • Correct warning is issued if `clip_pct` exceeds `clip_warn_pct`. <br> • Correct error (abort) is thrown if `clip_pct` exceeds `clip_abort_pct` AND `allow_clip` is `FALSE`. <br> • If `allow_clip` is `TRUE`, only a warning is issued even if `clip_pct > clip_abort_pct`.                                                                                                                                               | Q_Schema 1 Ans         |
| QNT-S1-11| **`forward_step.quant`: Hard Clipping (Global Scope)**          | • For `scale_scope="global"`: <br>  - After counting clipped samples, apply hard clipping: <br>    `Q_X_clipped <- Q_X` <br>    `Q_X_clipped[Q_X_clipped < 0] <- 0L` <br>    `Q_X_clipped[Q_X_clipped > (2^bits - 1)] <- (2^bits - 1L)` <br>  - The `Q_X_clipped` is then converted to `storage_type_str` and written.                                                                                                                         | • Output quantized values are strictly within the range $[0, 2^{bits}-1]$.                                                                                                                                                                                                                                                                                                                                                   | §2.2 (Item 4)          |
| QNT-S1-12| **Unit Tests for Clipping & Global Scope**                      | • Test non-finite value check. <br> • Test clipping counter and percentage calculation for global scope. <br> • Test warning/abort logic for global scope with various `clip_pct` and `allow_clip` values. <br> • Test that hard clipping correctly limits the output range. <br> • Test round-trip fidelity for `scale_scope="global"` with `bits=1` and `bits=16`.                                                                             | • All unit tests pass. <br> • Clipping behavior and global scope quantization are verified.                                                                                                                                                                                                                                                                                                                                | §5 (Unit-test)         |

---

**Sprint 1 - Definition of Done:**

*   `quant.schema.json` is updated with `allow_clip`, `report_path`, and examples.
*   `forward_step.quant` writes quantized data as `uint8` or `uint16` to HDF5 and attaches the `quant_bits` HDF5 attribute.
*   `invert_step.quant` reads the `quant_bits` attribute (with fallback and validation) to correctly interpret data.
*   Scale and offset parameters are stored as `float32`.
*   `forward_step.quant` includes a check for non-finite input values.
*   For `scale_scope="global"`, clipping is counted, and warnings/aborts are triggered based on `clip_pct` and the `allow_clip` parameter. Output values are hard-clipped.
*   Unit tests cover all implemented functionalities for efficient storage, `quant_bits` attribute handling, non-finite checks, and global scope clipping.
*   The `quant` transform (for global scope) is more robust and storage-efficient.

Sprint 2 will then tackle the block-wise processing for `scale_scope="voxel"`, the full quantization report generation, and the associated helper functions and documentation.

Okay, here are the granular tickets for **Sprint 2** of implementing the "Rock-Solid `quant` Transform." This sprint builds directly on Sprint 1, focusing on the memory-efficient block-wise processing for `scale_scope="voxel"`, implementing the detailed quantization report, and adding related helper functions and documentation.

**Assumptions for Sprint 2:**
*   Sprint 1 deliverables are complete and stable: `quant` uses narrow integer storage, `quant_bits` HDF5 attribute, float32 for scale/offset, and basic global scope clipping logic is in place.
*   LNA HDF5 helpers can create empty, chunked datasets and write to hyperslabs (or `hdf5r` direct equivalents are used).

---

# Rock-Solid `quant` Transform - Sprint 2 Implementation Tickets

**Epic QNT-S2-E1: Block-Wise Processing for `scope="voxel"`**
*Goal: Implement memory-efficient block-wise (slab-wise) processing for `forward_step.quant` when `scale_scope="voxel"`, ensuring correct calculation and storage of per-voxel parameters.*

| #        | Ticket                                                                 | Description / Deliverables                                                                                                                                                                                                                                                                                                                                                                                                                                                    | Acceptance Criteria                                                                                                                                                                                                                                                                                                                                                                                                                                                        | Appendix Ref           |
| :------- | :--------------------------------------------------------------------- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :--------------------- |
| QNT-S2-1 | **Internal Helper: `lna:::auto_block_size()`**                           | • Implement `lna:::auto_block_size(spatial_dims, element_size_bytes, target_slab_bytes = 64e6)`. <br> • Takes spatial dimensions of the data (e.g., `dim(x)[1:3]`) and target memory for a slab. <br> • Returns a list with `slab_dims` (e.g., `c(X_chunk, Y_chunk, Z_slab_depth)`) and `iterate_slabs_info` (e.g., number of slabs along each dim). <br> • Prioritizes making `Z_slab_depth=1` if possible, then adjusts `X_chunk`, `Y_chunk`.                  | • Helper function returns sensible slab dimensions that keep `prod(slab_dims) * element_size_bytes * num_timepoints` (for `quantized_values`) and `prod(slab_dims) * param_element_size_bytes` (for scale/offset) under or near `target_slab_bytes`. <br> • Handles various input `spatial_dims`.                                                                                                                                                                       | §2 (Forward Step Impl) |
| QNT-S2-2 | **`forward_step.quant`: Pre-create HDF5 Datasets (Voxel Scope)**       | • When `scale_scope="voxel"`: <br>  - Before block processing, pre-create empty, chunked HDF5 datasets using LNA HDF5 helpers (or direct `hdf5r`): <br>    1. `quantized_values` (shape `dim(x)`, type `uint8`/`uint16` based on `bits`). <br>    2. `quant_scale` (shape `dim(x)[1:spatial_dims]`, type `float32`). <br>    3. `quant_offset` (shape `dim(x)[1:spatial_dims]`, type `float32`). <br>  - Use chunking strategy from Q_BlockProcessing_1 answer (e.g., `c(64,64,1,dim(x)[TimeDim])` for `quantized_values`, `c(64,64,1)` for params, if processing Z-slabs of depth 1). | • Correctly named, shaped, typed, and chunked HDF5 datasets are created before block processing begins. <br> • HDF5 file handles to these datasets are managed appropriately (opened for creation, then accessed for slab writes).                                                                                                                                                                                                                             | Q2 Ans                 |
| QNT-S2-3 | **`forward_step.quant`: Block-wise Loop for Voxel Scope**              | • When `scale_scope="voxel"`: <br>  - Use `lna:::auto_block_size()` to determine slab iteration. <br>  - Loop through spatial slabs of the input data `X`. <br>  - For each `data_block` (e.g., `N_vox_in_slab x Time` matrix): <br>    - Call/adapt `.quantize_voxel()` (renamed to `.quantize_voxel_block()`) to calculate `scale_block` and `offset_block` *for voxels in this block only*. <br>    - Perform quantization $Q_{block}$ for `data_block`.          | • Loop iterates correctly over all spatial data. <br> • `.quantize_voxel_block()` is called for each slab with correct data.                                                                                                                                                                                                                                                                                                                         | Q2 Ans                 |
| QNT-S2-4 | **`forward_step.quant`: Write Slabs to HDF5 (Voxel Scope)**            | • Inside the block-wise loop: <br>  - Write the slab's quantized data $Q_{block}$ to the correct hyperslab in the pre-created `quantized_values` HDF5 dataset. <br>  - Write the slab's `scale_block` to the correct hyperslab in `quant_scale` HDF5 dataset. <br>  - Write the slab's `offset_block` to the correct hyperslab in `quant_offset` HDF5 dataset. <br>  - Use LNA HDF5 helpers or `hdf5r` direct hyperslab writing (`dset_obj$write(data, hyperslab_args)`). | • Data for each processed slab is correctly written to the appropriate slice of the corresponding HDF5 datasets. <br> • No large intermediate R objects (full scale/offset arrays or full quantized array) are created in memory.                                                                                                                                                                                                      | Q2 Ans                 |
| QNT-S2-5 | **`forward_step.quant`: Accumulate Clipping & Stats (Voxel Scope)**  | • Inside the block-wise loop (or after processing each block by `.quantize_voxel_block()`): <br>  - Accumulate total `n_clipped_total` across all blocks. <br>  - For the quantization report (implemented in E2): <br>    - Accumulate min/max/sum/sum_sq of per-voxel `scale` values. <br>    - Accumulate min/max/sum/sum_sq of per-voxel `offset` values. <br>    - (Optional) Accumulate counts for global histogram of quantized values.                 | • `n_clipped_total` is correctly accumulated. <br> • Statistics for per-voxel scale/offset are correctly accumulated for later use in the report (e.g., for calculating mean/sd).                                                                                                                                                                                                                                                           | Q2 Ans, Suggestion 2   |
| QNT-S2-6 | **Unit Tests for Block-Wise Voxel Scope Processing**                 | • Test `lna:::auto_block_size()` helper with various dimensions. <br> • Test pre-creation of HDF5 datasets with correct chunking. <br> • Mock HDF5 write functions to verify slab iteration and correct hyperslab indexing. <br> • Test end-to-end `forward_step.quant` with `scale_scope="voxel"` on moderately large synthetic data: <br>   - Verify output HDF5 datasets match results from a (slower) non-block-wise in-memory calculation. <br>   - Monitor memory usage during test (qualitatively or with tools if CI supports). | • Block-wise processing for `scope="voxel"` produces identical numerical results to an in-memory approach. <br> • Memory usage remains bounded (does not scale with total N_voxels for intermediate arrays). <br> • Clipping counts are correct.                                                                                                                                                                  | -                      |

**Epic QNT-S2-E2: Quantization Report & Final Polish**
*Goal: Implement the detailed JSON quantization report, associated helper function, and finalize documentation and tests for the `quant` transform.*

| #        | Ticket                                                              | Description / Deliverables                                                                                                                                                                                                                                                                                                                                                                                                                                                   | Acceptance Criteria                                                                                                                                                                                                                                                                                                                                                                                                                                                        | Appendix Ref                 |
| :------- | :------------------------------------------------------------------ | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :--------------------------- |
| QNT-S2-7 | **`forward_step.quant`: Generate Quantization Report Content**      | • At the end of `forward_step.quant` (after all blocks if voxel-scope): <br>  - Assemble the `quant_report` list: <br>    - `report_version = "1.0"`. <br>    - `clipped_samples_count`, `clipped_samples_percentage`. <br>    - `input_data_range: [global_min_orig, global_max_orig]`. <br>    - `effective_step_size`: single `scale_val` if global; or `scale_stats = {min, max, mean, sd}` if voxel scope (calculated from accumulated stats). <br>    - `effective_offset`: single `offset_val` if global; or `offset_stats = {min, max, mean, sd}` if voxel scope. <br>    - `estimated_snr_db`: (e.g., $10 \cdot \log_{10}(\text{var}(X_{global}) / (\text{mean_scale_global}^2 / 12))$). For voxel scope, estimate from a 1% random subset of voxels. <br>    - (Optional) `histogram_quantized_values: {breaks: [...], counts: [...]}` (global, 64 bins). | • `quant_report` list is correctly populated with all specified fields. <br> • Stats for `scale_scope="voxel"` (mean/sd of scale/offset) are accurate. <br> • SNR estimate is calculated. <br> • Histogram (if implemented) has correct structure.                                                                                                                                                              | Suggestion 2, §2 (Report)    |
| QNT-S2-8 | **`forward_step.quant`: Write Quantization Report to HDF5**         | • Serialize the `quant_report` list to JSON: `json_report_str <- jsonlite::toJSON(quant_report, auto_unbox=TRUE, pretty=TRUE)`. <br> • Compress: `gzipped_report <- memCompress(json_report_str, type="gzip")`. <br> • Determine report path: `report_hdf5_path <- paste0("/transforms/", desc$descriptor_basename, "_report.json")` (where `descriptor_basename` is "00_quant" etc.). <br> • Write `gzipped_report` (raw vector) to `report_hdf5_path` as a variable-length string or opaque byte dataset. Attach HDF5 attribute `compression = "gzip"` to this dataset. <br> • Store `report_hdf5_path` in `desc$params$report_path`. Update plan. | • Gzipped JSON report is written to the conventional HDF5 path. <br> • Dataset has `compression="gzip"` attribute. <br> • `desc$params$report_path` is correctly set in the main quant descriptor.                                                                                                                                                                                                                 | §2 (Report), Suggestion 1    |
| QNT-S2-9 | **Helper Function: `lna_get_transform_report()`**                   | • Implement generic `lna_get_transform_report(lna_file, transform_index_or_name)`. <br>  - Opens LNA file, reads main descriptor for specified transform. <br>  - Looks for `params$report_path`. If missing, tries conventional path `<descriptor_basename>_report.json`. <br>  - Reads the report dataset. Checks for `compression="gzip"` attribute; if present, calls `memDecompress()`. <br>  - Parses JSON string to list. Returns list. <br> • Export this helper. | • `lna_get_transform_report()` retrieves and correctly decompress/parses the quant report. <br> • Handles missing `report_path` by trying conventional name. <br> • Errors gracefully if report not found or malformed. <br> • Specific `lna_get_quant_report(lna_file, ...)` can be a thin wrapper if desired.                                                                                                  | §2 (Report), Suggestion 1    |
| QNT-S2-10| **Full Unit Test Suite for `quant` Transform**                      | • Consolidate and expand unit tests from Sprint 1 and Sprint 2. <br> • **Cover all checklist items:** uint8/16 storage, `quant_bits` attr, clipping logic (all `allow_clip` branches, warn/abort thresholds), block-wise voxel-scope processing (vs. in-memory oracle), inverse logic (attr read, validation), non-finite checks, report generation (content and file I/O), `lna_get_transform_report()`. <br> • Include tests for `bits=1`, `bits=16`. <br> • Test legacy file reading (missing `quant_bits` attr). | • Comprehensive test suite passes. <br> • All specified behaviors and edge cases are covered. <br> • `quant` transform is demonstrably "rock-solid."                                                                                                                                                                                                                                                                       | §3 (Checklist), §5 (Tests) |
| QNT-S2-11| **Documentation: `quant` Help Page & Vignette**                   | • Update roxygen docs for `forward_step.quant`/`invert_step.quant` (if they are directly user-visible, otherwise for the DSL verb). <br> • Create/update man page for `quant()` DSL verb, detailing all params including `allow_clip`. <br> • Create a short vignette: "Choosing `bits` and Interpreting the Quantization Report". Include example `lna_get_transform_report()` output (pretty-printed) and a sample histogram plot. <br> • In DSL docs, show `allow_clip=TRUE` example. | • Documentation is clear, accurate, and comprehensive. <br> • Vignette provides practical guidance.                                                                                                                                                                                                                                                                                                                        | §6 (Docs)                  |
| QNT-S2-12| **NEWS.md Update & Deprecation Note**                             | • Add entry to `NEWS.md` for the "rock-solid quant" implementation, highlighting key improvements (storage, clipping, report). <br> • If old `method="sd"` with `center=FALSE` had asymmetric behavior, note this change/fix in NEWS as a potential (minor) behavioral difference for users of very old versions.                                                                                                                                           | • `NEWS.md` accurately reflects changes.                                                                                                                                                                                                                                                                                                                                                                             | §7 (Release)               |

---

**Sprint 2 - Definition of Done:**

*   `forward_step.quant` with `scale_scope="voxel"` uses memory-efficient block-wise processing to calculate and store per-voxel scale/offset parameters and quantized data, writing to pre-created, correctly chunked HDF5 datasets.
*   The full JSON quantization report (including clipping stats, effective parameter stats for voxel scope, SNR, optional histogram) is generated, GZIP-compressed, and stored in HDF5 per the standardized naming convention. `desc$params$report_path` is set.
*   The `lna_get_transform_report()` helper function can retrieve, decompress, and parse these reports.
*   All clipping logic (counting, warning, aborting based on `allow_clip` and thresholds) is fully implemented and tested for both global and voxel scopes.
*   A comprehensive unit test suite validates all aspects of the `quant` transform, including storage types, attribute handling, global/voxel scopes, block-wise processing correctness, clipping, report generation, and legacy file considerations.
*   Documentation (man pages, new vignette) for the `quant` transform and its report is complete.
*   `NEWS.md` is updated.
*   The `quant` transform is considered "rock-solid," dependable, efficient, and informative.

This two-sprint plan should result in a high-quality, core `quant` transform that meets all specified requirements.