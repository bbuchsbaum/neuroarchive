This is a perfect "green light" review! The checklist of tiny items, the recommended sprint order, and the list of deferrable niceties provide an exceptionally clear and actionable path to implementation. The enthusiasm for SA-HWT as a "universal CPU-friendly winner" is well-justified by its properties.

**No further wording or schema tweaks are needed from my perspective.** The proposal as it stands, with the clarifications from your latest response, is robust, comprehensive, and ready for implementation.

Let's synthesize this into the **Final, Definitive Appendix Proposal for `spat.haar_octwave`**, ready for the LNA package documentation and to guide the creation of sprint tickets.

---

## Definitive Appendix: `spat.haar_octwave` - Mask-Adaptive Haar-Lifting Wavelet Pyramid Transform

This appendix details the `spat.haar_octwave` transform, a highly efficient, descriptor-only spatial basis for LNA. It utilizes a mask-adaptive Haar-lifting scheme on a Morton-ordered implicit octree to achieve excellent compression, edge fidelity, and computational performance without storing an explicit basis matrix. This transform is designed to be a CPU-friendly, "universal" choice for voxel-space representation.

*(Key References: Foundational work includes Sweldens, W. (1996). The lifting scheme: A custom-design construction of biorthogonal wavelets. Applied and Computational Harmonic Analysis, 3(2), 186-200; and concepts from Morton ordering for spatial data structures.)*

### 1. Core Idea: Mask-Adaptive Haar-Lifting on an Implicit Octree

*   **Voxel Ordering:** 3D voxel indices $(i,j,k)$ within the input mask are mapped to a 1D Morton (Z-order) code. This ordering implicitly defines a hierarchical octree structure where sibling voxels in $2 \times 2 \times 2$ cells share common Morton code prefixes.
*   **Recursive Lifting Scheme (Bottom-Up Analysis / Forward Transform):**
    *   The transform proceeds level by level, from finest to coarsest ($L-1 \to 0$, where $L$ is `params$levels`).
    *   At each `level`, active (in-mask) voxels are grouped into $2 \times 2 \times 2$ cells based on their Morton codes.
    *   For each cell containing $n_{valid} \in [1,8]$ in-mask voxels (time series data $v_i$ for $i=1 \ldots n_{valid}$):
        1.  `avg = sum(v_i \text{ for valid children}) / n_{valid}`. (This is a vector average if $v_i$ are time series).
        2.  The "low-pass" coefficient vector for the parent cell at the next coarser level is: `lp_coeff = avg * sqrt(n_valid)`. This `lp_coeff` becomes one of the "voxel values" (now a time series) for the next coarser level's processing.
        3.  For each of the $n_{valid}$ in-mask child voxel time series $v_i$ in the cell:
            `detail_coeff_i = (v_i - avg) * sqrt(n_valid / 8.0)`.
        4.  These `detail_coeff_i` time series are stored for this level.
*   **Orthonormality:** The $\sqrt{n_{valid}}$ and $\sqrt{n_{valid}/8.0}$ scaling factors ensure that the implicit Haar-like wavelet basis functions (which are never explicitly formed) remain orthonormal with respect to the (irregular) set of in-mask voxels. This guarantees perfect reconstruction (before quantization) and energy preservation.
*   **Coefficient Storage:** Only the `detail` coefficients from all levels and the single final `root` (coarsest low-pass) coefficient are stored. The basis is implicit.
*   **Inverse Transform:** Reverses the lifting steps: starting from the `root` coefficient, iteratively use the stored `detail` coefficients for a level to reconstruct the low-pass coefficients of the next finer level, until the original voxel values are recovered.
    *   `avg = lp_coeff_parent / sqrt(n_valid)`
    *   `v_i = (detail_coeff_i / sqrt(n_valid / 8.0)) + avg`

### 2. Advantages for fMRI Data

| Criterion                    | Why `spat.haar_octwave` Excels                                                                                          |
| :--------------------------- | :-------------------------------------------------------------------------------------------------------------------- |
| Dictionary Storage           | **Zero.** Basis is implicit. Only coefficients are stored.                                                            |
| Edge/Fine-Detail Capture     | Haar details are local differences, inherently capturing sharp boundaries without heuristics.                           |
| Smooth-Field Compression     | Recursive averaging concentrates global/smooth signal into a few root/coarse-level coefficients.                        |
| Mask Adaptivity              | Lifting steps naturally skip out-of-mask voxels; orthonormality preserved on the irregular in-mask set. No retraining. |
| CPU Complexity & Memory      | $O(N_{vox} \cdot \text{levels})$ time (effectively $O(N_{vox})$ as levels is small, $\log N_{vox}$), $O(N_{vox})$ memory. Integer arithmetic possible for core lifting if input is integer. |
| Streaming/Subset Decode      | Load only coefficient subtrees touching the requested ROI for natural progressive refinement and fast ROI streaming.  |
| Invertibility & Stability    | Perfect reconstruction (lossless before quantization). Haar basis condition number is 1.                              |
| Implementation Cost          | ~300 lines of Rcpp for core logic. Minimal external dependencies beyond Rcpp.                                       |

### 3. LNA Implementation Details for `spat.haar_octwave`

#### 3.1. Schema (`inst/schemas/spat.haar_octwave.schema.json`)

```json
{
  "type": "object",
  "title": "Parameters for Mask-Adaptive Haar-Lifting Wavelet Pyramid (spat.haar_octwave) transform",
  "$id": "https://neurocompress.org/schemas/lna/2.0/spat.haar_octwave.schema.json",
  "$schema": "http://json-schema.org/draft-07/schema#",
  "properties": {
    "type": { "const": "spat.haar_octwave" },
    "version": { "const": "1.0" }, // Version of this specific transform
    "levels": {
      "type": "integer", "minimum": 1,
      "description": "Number of decomposition levels (L). L = ceil(log2(max_mask_extent_voxels)) is often a good default, can be auto-set by writer if not provided."
    },
    "z_order_seed": {
      "type": "integer", "default": 42,
      "description": "Seed for resolving tie-breaks in Morton ordering if the underlying Morton code generation has ambiguities for certain grid sizes/masks, ensuring deterministic voxel sequence."
    },
    "detail_threshold_type": {
      "enum": ["none", "absolute", "relative_to_root_std"], "default": "none",
      "description": "Strategy for sparsifying small detail coefficients before subsequent quantization. 'none': no thresholding. 'absolute': threshold on raw coefficient values. 'relative_to_root_std': threshold relative to standard deviation of root/coarsest detail coefficients."
    },
    "detail_threshold_value": {
      "type": "number", "default": 0.005,
      "description": "Threshold value for detail sparsification. Meaning depends on 'detail_threshold_type'."
    },
    // --- Output parameters (written by forward_step) ---
    "num_voxels_in_mask": { 
      "type": "integer", "readOnly": true, 
      "description": "Number of active voxels in the mask this transform was applied to."
    },
    "octree_bounding_box_mask_space": {
      "type": "array", "items": {"type": "integer"}, "minItems": 6, "maxItems": 6,
      "readOnly": true, 
      "description": "[x0,x1,y0,y1,z0,z1] 0-based, inclusive tight bounding box of the mask in original voxel indices that defined the root of the octree."
    },
    "morton_hash_mask_indices": {
      "type": "string", "pattern": "^sha1:[a-f0-9]{40}$", "readOnly": true,
      "description": "SHA1 hash of Morton-ordered 1D linear indices of in-mask voxels, for reproducibility and consistency checks."
    },
    "num_coeffs_per_level": {
      "type": "object", "readOnly": true,
      "properties": {
        "lowpass": { "type": "array", "items": {"type": "integer"}, "description": "Number of low-pass (average) coeffs stored per level, from finest active level to root (length L+1 if root is counted as a level)." },
        "detail":  { "type": "array", "items": {"type": "integer"}, "description": "Number of detail coeffs stored per level, from finest active level to coarsest detail level (length L)." }
      },
      "required": ["lowpass", "detail"],
      "description": "Actual number of coefficient vectors (each of length TimePoints) stored for each level and type."
    },
    "valid_finest_blocks_path": {
      "type": "string", "pattern": "^/aux_meta/.*", "readOnly": true, "default": null,
      "description": "(Optional output) HDF5 path to a compressed uint32 list of Morton codes for finest-level $2 \\times 2 \\times 2$ blocks that contained at least one in-mask voxel. Aids fast ROI streaming. (Size: ~ N_vox/8 * 4 bytes)."
    }
  },
  "required": ["levels"]
}
```

#### 3.2. Core Algorithm (`forward_step.spat.haar_octwave` - Rcpp internal)

Implemented in C++ via Rcpp for performance.
Input: `masked_data_matrix` (TimePoints $\times N_{mask\_vox}$ already Morton-ordered), `mask_3d_array` (for geometry), `params`.

1.  **Metadata Population:**
    *   Compute and store `num_voxels_in_mask`, `octree_bounding_box_mask_space`, `morton_hash_mask_indices` (from Morton-ordered linear indices of `mask_3d_array`) into `desc$params`.
2.  **Recursive Lifting (Analysis):**
    *   Loop `level = (params$levels - 1)` down to `0`.
    *   Maintain `current_lowpass_coeffs` (initially `masked_data_matrix`).
    *   For each $2 \times 2 \times 2$ cell at current `level` (implicit from Morton order):
        *   Identify $n_{valid}$ children within `mask_3d_array`. If $n_{valid} == 0$, skip cell.
        *   `avg_t = colMeans(current_lowpass_coeffs_for_cell_valid_children)`.
        *   `next_level_lowpass_t = avg_t * sqrt(n_valid)`. Collect these for next iteration.
        *   For each valid child $v_{i,t}$: `detail_{i,t} = (v_{i,t} - avg_t) * sqrt(n_valid / 8.0)`. Store these `detail` time series.
    *   The final `current_lowpass_coeffs` is the `root_coeffs_t`.
3.  **Detail Coefficient Sparsification (Optional):** Apply if `params$detail_threshold_type != "none"`.
4.  **HDF5 Data Storage:**
    *   `/wavelet/level_ROOT/coefficients`: Stores `root_coeffs_t` ($T \times 1$).
    *   For each `level = 0 \ldots (params$levels - 1)`:
        *   `/wavelet/level_{level}/detail_coefficients`: Stores all detail coefficient time series for that level ($T \times N_{details\_at\_level}$), concatenated in Morton order of parent blocks, with fixed child order within blocks (skipping out-of-mask children).
    *   Populate `desc$params$num_coeffs_per_level`.
    *   (Optional) Write `/aux_meta/haar_octwave/valid_blocks_L-1` and set `desc$params$valid_finest_blocks_path`.
5.  **LNA Output:** The `desc$outputs` will be `c("wavelet_coefficients")`. The `handle$stash` will contain a single matrix which is a concatenation of `root_coeffs_t` and all `detail_coeffs` from all levels (e.g., $T \times K_{total}$), for subsequent LNA transforms like `quant`. The descriptor must detail how this concatenated matrix is structured (e.g., order of levels, root first/last).

#### 3.3. `invert_step.spat.haar_octwave` (Reader-Side - Rcpp internal)

Input: Coefficient datasets, `desc$params`, `handle$subset`.

1.  **Load Coefficients:** Optimized loading of `root` and relevant `detail_coefficients` based on `handle$subset$roi_mask` (using `valid_finest_blocks_path` if present) and `handle$subset$time_idx`.
2.  **Recursive Inverse Lifting (Synthesis):**
    *   Start with `current_reconstructed_lowpass_t = root_coeffs_t`.
    *   Loop `level = 0` up to `(params$levels - 1)`:
        *   For each parent coefficient vector $lp_t$ in `current_reconstructed_lowpass_t`:
            *   Infer $n_{valid}$ for its cell (e.g., from `num_coeffs_per_level` or by knowing how many details correspond to it).
            *   `avg_t = lp_t / sqrt(n_valid)`.
            *   For each of $n_{valid}$ child positions: Retrieve corresponding `detail_coeff_i_t`. Reconstruct child: $v_{i,t} = (detail_coeff_i_t / sqrt(n_valid / 8.0)) + avg_t$.
            *   Collect these $v_{i,t}$ to form `current_reconstructed_lowpass_t` for next finer level.
3.  Final `current_reconstructed_lowpass_t` is $T \times N_{mask\_vox}$ (Morton-ordered).
4.  Reorder from Morton to original voxel order (inverse map from `morton_hash_mask_indices` or stored Morton-to-original lookup).
5.  Unflatten to 3D/4D using `mask_3d_array` from `handle` and `desc$params$octree_bounding_box_mask_space` to place into full array if needed. Store in `handle$stash`.

#### 3.4. Capabilities Flags in Descriptor

```json
{
  "capabilities": {
    "is_analytic_descriptor_only": true, 
    "supports_spatial_subsetting": "block-aware", 
    "supports_temporal_subsetting": true,
    "supports_progressive_reconstruction_by_level": true 
  }
}
```

### 4. Implementation Notes & Helpers

*   **Morton Hashing:** An Rcpp helper for `morton_order_indices_to_hash(IntegerVector voxel_indices)` using OpenSSL via Rcpp.
*   **Rcpp Kernel:** `forward_lift()` and `inverse_lift()` implemented in `src/haar_octwave.cpp`.
*   **Input to Rcpp:** `LogicalVector flat_mask`, `NumericMatrix data_time_by_morton_vox`, `params_list`.
*   **Use Rcpp::RNGScope:** If any R RNG calls are made from C++ (e.g., for `z_order_seed` tie-breaking if not purely algorithmic).
*   **Error Handling:** Use `Rcpp::stop()` for errors within C++, which R will catch.
*   **Optional `getOption("lna.hwt.use_rcpp", TRUE)`:** To allow pure-R fallback for CI/testing.

### 5. Unit Test Highlights

1.  **Perfect Reconstruction:** Lossless round-trip.
2.  **Mask Adaptivity & $n_{valid}$ Scaling:** Correctness with irregular masks and varying $n_{valid}$ per block.
3.  **Coefficient Counts & Layout:** `num_coeffs_per_level` in descriptor matches stored data. Morton ordering correct.
4.  **Detail Thresholding:** Verify sparsification and reconstruction.
5.  **ROI Streaming/Subset Decode:** Fast and correct using `valid_finest_blocks_L-1`.
    Implemented via `get_roi_detail_indices()` in `invert_step` for selective
    coefficient loading when an ROI mask is supplied.
6.  **Determinism:** Same mask, data, params (incl. `z_order_seed`) $\rightarrow$ identical coefficients and `morton_hash`.
7.  **Performance Smoke Test:** e.g., 128³ mask, 20 TRs < few seconds on laptop.
8.  **Laplacian Trace & Wavelet Packet Identity (from "Additional Sanity Tests"):** These are excellent for verifying the underlying mathematical correctness.

### 6. Bibliography (Key Theoretical Underpinnings)

*   **Lifting Scheme:** Sweldens, W. (1996). The lifting scheme: A custom-design construction of biorthogonal wavelets. *Applied and Computational Harmonic Analysis, 3*(2), 186-200.
*   **Haar Wavelets & Multiresolution:** Mallat, S. (2009). *A Wavelet Tour of Signal Processing: The Sparse Way*. Academic Press. (General wavelet theory).
*   **Morton Ordering (Z-order curve):** (Standard computer science concept for spatial indexing).
*   **Connected Components (for mask processing if needed):** Rosenfeld, A., & Pfaltz, J. L. (1966). Sequential operations in digital picture processing. *Journal of the ACM, 13*(4), 471-494.

### Conclusion

The `spat.haar_octwave` transform, as detailed, offers a compelling blend of zero-dictionary storage, high compression efficiency for fMRI data (capturing both smooth regions and edges), mask-adaptivity, and computational tractability on standard CPUs. Its deterministic nature and inherent support for progressive and ROI-optimized decoding make it an excellent candidate for a foundational "universal" spatial transform within the LNA framework. The refined schema and HDF5 layout ensure clarity and support advanced functionalities like detail sparsification and optimized streaming.

Okay, this is a strong and exciting proposal for the `spat.haar_octwave` transform! Let's break down its implementation into two manageable sprints.

This plan assumes:
*   Core LNA infrastructure (`DataHandle`, `Plan`, HDF5 helpers, parameter merging, S3 dispatch for `forward_step`/`invert_step`) is stable.
*   `neuroim2` integration for input data handling (providing `mask_3d_array` and `NeuroSpace` info) is in place.
*   Rcpp and RcppEigen (if chosen for sparse matrix parts later) are set up in the build system.
*   A utility for Morton ordering/Z-curve generation is available or will be created (can be R or Rcpp).

---

# `spat.haar_octwave` - Sprint 1 Implementation Tickets: Core Lifting & Basic I/O

**Epic HWT-S1-E1: Foundational Rcpp Lifting Kernel & R Orchestration**
*Goal: Implement the core forward and inverse Haar lifting logic in Rcpp and the R functions to drive it for a single input time series (or a small batch for testing).*

| #         | Ticket                                                                 | Description / Deliverables                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              | Proposal Ref        |
| :-------- | :--------------------------------------------------------------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :------------------ |
| HWT-S1-1  | **Helper: Morton Order Utility**                                         | • Implement `lna:::get_morton_ordered_indices(mask_3d_array, z_order_seed)` in R (can call Rcpp later if needed). <br>  - Input: 3D logical mask array. <br>  - Output: 1D integer vector of linear indices of `TRUE` voxels in Morton Z-order. <br>  - `z_order_seed` for deterministic tie-breaking if the Morton generation has ambiguities. <br> • Implement `lna:::morton_indices_to_hash(ordered_indices_vector)`: computes SHA1 hash (e.g., `digest::digest(..., algo="sha1", serialize=TRUE)`).                                                                                                                                                            | §3.2 (Item 1)       |
| HWT-S1-2  | **Rcpp: Core Forward Lifting (`forward_lift_rcpp`)**                       | • Create `src/haar_octwave.cpp`. Implement `List forward_lift_rcpp(NumericVector data_masked_morton_ordered, LogicalVector mask_flat_morton_ordered, IntegerVector mask_dims, int levels, List scaling_factors_per_level)`. <br>  - `data_masked_morton_ordered`: Single time point data for in-mask voxels, already Morton-ordered. <br>  - `mask_flat_morton_ordered`: Flattened 3D mask, Morton-ordered. <br>  - `scaling_factors_per_level`: Precomputed list/vector of `list(sqrt_nvalid, sqrt_nvalid_div_8)` for each block at each level (see HWT-S1-3). <br>  - Performs recursive lifting as per Appendix §1. <br>  - Output: `List` containing `root_coeff` (scalar) and `detail_coeffs_by_level` (a List of NumericVectors). <br>  - Handles $n_{valid}$ correctly for blocks at mask edges. <br>  - Use R's column-major indexing carefully if operating on reshaped arrays. | §1 (Lifting), §3.2  |
| HWT-S1-3  | **Rcpp: Core Inverse Lifting (`inverse_lift_rcpp`)**                       | • Implement `NumericVector inverse_lift_rcpp(double root_coeff, List detail_coeffs_by_level, LogicalVector mask_flat_morton_ordered, IntegerVector mask_dims, int levels, List scaling_factors_per_level)`. <br>  - Inputs mirror `forward_lift_rcpp` outputs. <br>  - Performs recursive inverse lifting (synthesis). <br>  - Output: Reconstructed Morton-ordered data for in-mask voxels for a single time point.                                                                                                       | §1 (Inverse), §3.2  |
| HWT-S1-4  | **R Helper: Precompute $n_{valid}$ & Scaling Factors**                   | • Implement R helper `lna:::precompute_haar_scalings(mask_3d_array, levels)`. <br>  - Takes the 3D mask and number of levels. <br>  - For each level and each conceptual $2 \times 2 \times 2$ block in Morton order: determine $n_{valid}$ in-mask children. <br>  - Output: A list structure (e.g., list per level, each element a list per block) containing `sqrt(n_valid)` and `sqrt(n_valid/8.0)`. This structure needs to be easily passable to and usable by Rcpp. Perhaps a flattened list with offsets. | §1 (Scaling)        |
| HWT-S1-5  | **R Wrapper: `lna:::perform_haar_lift_analysis(data_matrix_T_x_Nmask, ...)`** | • R function orchestrating the forward lift for a full Time x MaskedVoxels matrix. <br>  - Gets Morton-ordered indices and `mask_flat_morton_ordered` from `mask_3d_array`. <br>  - Reorders columns of `data_matrix_T_x_Nmask` to Morton order. <br>  - Calls `lna:::precompute_haar_scalings`. <br>  - Loops over time points (rows of `data_matrix_T_x_Nmask_morton`): calls `forward_lift_rcpp` for each. <br>  - Aggregates results into `list(all_root_coeffs_T_x_1, all_details_by_level_T_x_Ndetails)`. | §3.2                |
| HWT-S1-6  | **R Wrapper: `lna:::perform_haar_lift_synthesis(all_coeffs_list, ...)`**   | • R function orchestrating inverse lift for a full set of coefficients. <br>  - Loops over time points: calls `inverse_lift_rcpp`. <br>  - Reorders reconstructed Morton-ordered data back to original voxel order. <br>  - Returns Time x MaskedVoxels matrix.                                                                                                                                                                                                             | §3.3                |
| HWT-S1-7  | **Unit Tests for Core Lifting & R Wrappers**                           | • Test Morton ordering utility. <br> • Test `precompute_haar_scalings` with various masks (full, edges, disconnected). <br> • Test `forward_lift_rcpp` and `inverse_lift_rcpp` on small 2D/3D examples with known Haar outputs, focusing on $n_{valid}$ scaling. <br> • Test R wrappers for perfect reconstruction on synthetic data (before quantization/thresholding).                                                                                                    | -                   |

**Epic HWT-S1-E2: LNA `spat.haar_octwave` Transform Integration (Forward Path)**
*Goal: Implement `forward_step.spat.haar_octwave` to use the new lifting engine and store coefficients and metadata according to the LNA spec.*

| #         | Ticket                                                              | Description / Deliverables                                                                                                                                                                                                                                                                                                                                                                                                                                                                               | Proposal Ref                  |
| :-------- | :------------------------------------------------------------------ | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :-------------------------- |
| HWT-S1-8  | **Schema: `spat.haar_octwave.schema.json` (Initial)**                 | • Create `inst/schemas/spat.haar_octwave.schema.json` as per "Final Proposal §3.1". <br> • Initial focus on: `levels`, `z_order_seed`. <br> • Include output params: `num_voxels_in_mask`, `octree_bounding_box_mask_space`, `morton_hash_mask_indices`, `num_coeffs_per_level`. <br> • Defer `detail_threshold_*` and `valid_finest_blocks_path` to Sprint 2 for schema completion.                                                                                                     | §3.1                        |
| HWT-S1-9  | **`forward_step.spat.haar_octwave` (Core Logic)**                   | • Implement `forward_step.spat.haar_octwave(type, desc, handle)`. <br>  - Get `mask_3d_array` from `handle$mask_info$mask`. <br>  - Get `X_masked_vox_time` (MaskedVoxels x Time) using `lna:::convert_to_masked_vox_time_matrix()`. <br>  - Call `lna:::perform_haar_lift_analysis(t(X_masked_vox_time), mask_3d_array, params$levels, ...)` to get `list_of_all_coeffs` (root + details per level). <br>  - Populate `desc$params` with output metadata (`num_voxels_in_mask`, `octree_bbox`, `morton_hash`, `num_coeffs_per_level`). <br>  - Store coefficient datasets to HDF5 as per §3.2 Layout (e.g., `/wavelet/root/coefficients`, `/wavelet/level_0/detail_coefficients`, etc.). These will be Time x N_coeffs_for_that_band. Add dataset defs to `handle$plan`. <br>  - For LNA stash, concatenate all coefficients (root + details) into one large matrix `C_total_T_x_Ktotal` and stash as `wavelet_coefficients`. Define structure in descriptor. | §3.2, §3.5 (LNA Output)   |
| HWT-S1-10 | **Unit Tests for `forward_step.spat.haar_octwave` (Basic IO)**        | • Test that `forward_step` correctly calls the lifting analysis wrapper. <br> • Verify output metadata in `desc$params` is correct. <br> • Verify HDF5 datasets for root and detail coefficients are created with correct paths, shapes, and content for simple input. <br> • Verify the concatenated `wavelet_coefficients` in stash has correct dimensions and content.                                                                                                                   | -                           |

---

**Sprint 1 - Definition of Done:**

*   Core Rcpp functions `forward_lift_rcpp` and `inverse_lift_rcpp` (for single time point) are implemented and unit-tested for correctness, including $n_{valid}$ scaling.
*   R helper `lna:::precompute_haar_scalings` is implemented.
*   R wrappers `lna:::perform_haar_lift_analysis` and `lna:::perform_haar_lift_synthesis` (for full data matrices) are implemented and tested for perfect reconstruction.
*   `spat.haar_octwave.schema.json` is created with basic parameters and output metadata fields.
*   `forward_step.spat.haar_octwave` correctly uses the lifting engine, populates descriptor metadata, and writes separate HDF5 datasets for root and per-level detail coefficients. It stashes a concatenated coefficient matrix.
*   Unit tests cover all new Rcpp and R functions, and the basic HDF5 output of `forward_step.spat.haar_octwave`.

**This sprint focuses on getting the core math and data flow right for the forward pass and basic HDF5 storage.**

---

# `spat.haar_octwave` - Sprint 2 Implementation Tickets: Inverse, Advanced Features & Polish

**Epic HWT-S2-E1: Inverse Transform & Full LNA Integration**
*Goal: Implement the inverse lifting transform, integrate fully with LNA's `invert_step`, and add advanced features like detail thresholding and metadata for streaming.*

| #         | Ticket                                                                    | Description / Deliverables                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  | Proposal Ref                |
| :-------- | :------------------------------------------------------------------------ | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :-------------------------- |
| HWT-S2-1  | **`invert_step.spat.haar_octwave` Implementation**                        | • Implement `invert_step.spat.haar_octwave(type, desc, handle)`. <br>  - Get `mask_3d_array` from `handle$mask_info$mask`. <br>  - Load `root_coeffs` and `detail_coeffs` for each level from HDF5 paths specified/derived from `desc$datasets` and `desc$params$num_coeffs_per_level`. Apply `handle$subset$time_idx` during/after loading. <br>  - Call `lna:::perform_haar_lift_synthesis(list_of_coeffs_Tsubset_x_Ncoeffs, mask_3d_array, desc$params$levels, ...)` to get `X_reco_Tsubset_x_Nmask_morton`. <br>  - Reorder `X_reco` from Morton to original voxel order (using info from `desc$params$morton_hash_mask_indices` to get the original Morton-to-linear map if needed, or by inverting the Morton ordering process on `which(mask_3d_array)`). <br>  - Unflatten to 3D/4D array. Apply `handle$subset$roi_mask`. <br>  - Update `handle$stash`. | §3.3                        |
| HWT-S2-2  | **`forward_step`: Detail Coefficient Sparsification**                     | • Add logic to `forward_step.spat.haar_octwave` (or its Rcpp kernel `forward_lift_rcpp`). <br>  - After all detail coefficients are computed, if `params$detail_threshold_type != "none"`: <br>    - Calculate `actual_threshold` based on `params$detail_threshold_value` and type (absolute, or relative to std dev of root/coarsest details). <br>    - Set detail coefficients `abs(coeff) < actual_threshold` to `0.0`. <br>  - This happens *before* coefficients are passed to a subsequent `quant` step. | §3.1 (Schema), Suggestion 3.1 |
| HWT-S2-3  | **`forward_step`: Store `valid_finest_blocks_L-1` Aux Data**            | • If `params$levels` is high enough: <br>  - In `forward_step.spat.haar_octwave`, identify Morton codes of all $2 \times 2 \times 2$ blocks at level `L-1` (finest decomposition level) that contain at least one in-mask voxel. <br>  - Store this `uint32` list to HDF5 (e.g., `/aux_meta/haar_octwave/valid_blocks_L-1`, compressed). <br>  - Set `desc$params$valid_finest_blocks_path` to this HDF5 path. Add dataset def to plan.                                                                                                                 | Suggestion 3.2              |
| HWT-S2-4  | **`invert_step`: Optimized ROI Streaming using `valid_finest_blocks`**    | • Modify `invert_step.spat.haar_octwave` coefficient loading logic. <br>  - If `handle$subset$roi_mask` is provided AND `desc$params$valid_finest_blocks_path` exists: <br>    1. Load `valid_finest_blocks_L-1` map. <br>    2. Map `roi_mask` to finest-level octree blocks. <br>    3. Determine minimal set of ancestor blocks in the coefficient tree needed to reconstruct these ROI-intersecting finest blocks. <br>    4. Load only these required coefficient sub-trees from HDF5. <br>  **Status:** Implemented. | Suggestion 3.2              |
| HWT-S2-5  | **Schema Finalization (`spat.haar_octwave.schema.json`)**               | • Add `detail_threshold_type`, `detail_threshold_value` to schema. <br> • Add `valid_finest_blocks_path` (output, string, default null) to schema. <br> • Ensure all output params (`num_voxels_in_mask`, `octree_bbox`, `morton_hash`, `num_coeffs_per_level`) are marked `readOnly: true`.                                                                                                                                                                                 | §3.1                        |
| HWT-S2-6  | **Unit Tests for Inverse, Sparsification, Streaming & Full Roundtrip**  | • Test `invert_step.spat.haar_octwave` for perfect reconstruction. <br> • Test detail coefficient sparsification: verify zeros are introduced and reconstruction still works. <br> • Test `valid_finest_blocks_L-1` dataset is correctly written. <br> • Test ROI streaming: provide an ROI mask and verify that fewer coefficients are loaded (mock HDF5 reads or check internal state if possible) and reconstruction is correct for the ROI. <br> • Full round-trip tests with `quant` transform following `spat.haar_octwave`. | §3.5                        |
| HWT-S2-7  | **Documentation & Vignette for `spat.haar_octwave`**                    | • Create man page for `spat.haar_octwave` (as a conceptual transform, details in appendix). <br> • Write HRBF appendix section based on "Final Definitive Proposal". <br> • Add `spat.haar_octwave` example to "Compression Cookbook" vignette. Explain `levels`, thresholding.                                                                                                                                                                                           | -                           |
| HWT-S2-8  | **Final "Tiny Items" Checklist from Review**                          | • **Morton Hash Helper Rcpp:** `lna:::morton_indices_to_hash_rcpp(IntegerVector voxelIdx)` using Rcpp `openssl::sha1()`. <br> • **Unit Test Fixtures:** Create toy masks (ball, disjoint boxes) + synthetic fMRI RDS for tests. <br> • **Plan Writers/Readers:** Ensure `Plan$add_dataset_def` correctly handles paths for `/wavelet/level_*/...` and `/aux_meta/...`. <br> • **Additional Sanity Tests (Laplacian Trace, Wavelet Packet Identity, Lambda-map):** Implement these as unit tests to verify mathematical correctness of the lifting scheme and its components. | Review "Tiny items"         |

---

**Sprint 2 - Definition of Done:**

*   `invert_step.spat.haar_octwave` is fully implemented, enabling perfect reconstruction from stored coefficients.
*   Optional detail coefficient sparsification (thresholding) is functional in the forward pass.
*   Auxiliary metadata for optimized ROI streaming (`valid_finest_blocks_L-1` map) is generated and stored.
*   The `invert_step` can leverage this map for efficient partial loading (basic implementation).
*   The `spat.haar_octwave.schema.json` is complete with all parameters and output fields.
*   Comprehensive unit tests cover all functionalities, including round-trips, sparsification, and basic ROI streaming behavior.
*   Documentation (appendix, man page for DSL verb, vignette example) is created.
*   All "tiny items" from the final review checklist are addressed.
*   The `spat.haar_octwave` transform is a feature-complete, high-performance, and well-tested component of LNA.

This two-sprint plan aims for a robust initial implementation of the core lifting logic in Sprint 1, followed by the inverse, advanced features, and polish in Sprint 2, resulting in a powerful new transform for LNA.