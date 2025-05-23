You're absolutely right to push for this level of detail and foresight! Making `spat.oct_pca` (or any complex transform) "stand-alone" as much as possible is key for modularity, reusability, and easier maintenance, even if it initially lives within the main `lna` package.

This review focuses on ensuring the Octree PCA proposal is self-contained in its logic and data structures, minimizing dependencies on LNA-specific `DataHandle` internals beyond the defined input/output contract (primarily data matrices and parameter lists).

**Review of "Ready for Implementation" Octree PCA Proposal (with "Stand-Alone" Lens):**

The previous iteration was already very strong. The key aspects for stand-alone capability are:

1.  **Input/Output Data:** The transform takes `X (dense_mat; shape Voxels × Time)` and outputs `block_table`, local bases, and `coeffs_blocks`. These are standard data structures (arrays, tables).
2.  **Parameters:** Defined via a JSON schema, which is portable.
3.  **Core Logic (`build_octree`, `fast_svd`, blending):** These are algorithmic components that can be implemented independently of LNA's `DataHandle` or `Plan` objects during their execution.
4.  **HDF5 Interaction:** The HDF5 layout is well-defined. A stand-alone module would need its own HDF5 writing/reading helpers if not using LNA's, but the *structure* it writes to is the important part for LNA compatibility.

**Addressing Previous Questions/Suggestions with "Stand-Alone" in Mind:**

*   **Q_Schema_Refine_1 (Default for `k_local` vs. `split_var_thresh`):**
    *   **Stand-alone perspective:** The *algorithm* for Octree PCA needs a clear rule. The adaptive `k_star` logic (`k_star = min(which(cum_energy >= 1 - p$split_var_thresh), p$k_local)`) is internal to the algorithm.
    *   **LNA Schema:** It's fine for `k_local` to have a default (e.g., 8 or a higher cap like 64). If a user wants `split_var_thresh` to dominate, they can set `k_local` to a large enough value. The schema should document this interaction: "`k_local` acts as an upper cap on the number of components chosen adaptively by `split_var_thresh`."
    *   **Decision:** Keep `k_local` with a default in the schema. The core algorithm uses both.

*   **Q_Global_Stage_Refine_1 (Data for Global PCA):**
    *   **Stand-alone perspective:** The Octree PCA algorithm needs to know what data to use for its optional global stage.
    *   **LNA Schema:** Adding `global_pca_input: {"enum": ["block_means", "full_residuals"], "default": "block_means"}` to the schema *is* the correct way to make this choice explicit and configurable.
    *   **HDF5 Layout Reconciliation:**
        *   If `global_pca_input: "block_means"`:
            *   The "basis" is effectively on the block means. Input to PCA: `Mean_Activations_Per_Block x Time`. Output global basis: `Number_of_Leaf_Blocks x k_global`. Output global coefficients: `Time x k_global`. The HDF5 path `/spat/oct_pca/global_basis/matrix` would store this `Number_of_Leaf_Blocks x k_global` matrix. Reconstruction would involve multiplying global coeffs by this basis to get reconstructed *block mean time series*, which then need to be appropriately added back (this is more complex than adding back a voxel-space residual).
        *   If `global_pca_input: "full_residuals"`:
            *   Input to PCA: `Full_Residual_Matrix (Voxels_in_Mask x Time)`. Output global basis: `Voxels_in_Mask x k_global`. Output global coefficients: `Time x k_global`. Path `/spat/oct_pca/global_basis/matrix` stores `Voxels_in_Mask x k_global`. This is simpler for reconstruction.
    *   **Decision:** The HDF5 layout `/spat/oct_pca/global_basis/matrix` currently implies a voxel-space basis. To keep it simple and align with the "full residual matrix" being more direct for reconstruction, **let's specify that if `k_global > 0`, the global PCA is performed on the full voxel-space residuals $R = X - \hat{X}_{\text{block-recon}}$**. The option for "block_means" could be a future variant if a strong use case emerges, but it complicates reconstruction blending. This simplifies the schema (no `global_pca_input` needed for v1.0).

*   **Suggestion: Explicit `block_id_format_string`:**
    *   **Stand-alone perspective:** The core algorithm will generate integer block IDs. How these are mapped to HDF5 group names is an I/O detail.
    *   **LNA Schema:** This is a good idea for consistent HDF5 path generation. Add `params$block_id_hdf5_format_string: {"type": "string", "default": "%06d"}`. The LNA S3 methods would use this when creating HDF5 groups like `/basis/oct_pca/blocks/{formatted_block_id}/`.

*   **Suggestion: Store "Effective Overlap" in `block_table`:**
    *   **Stand-alone perspective:** The core algorithm computes with the *requested* overlap, clipped by mask/image boundaries. The blending weights (from `weight_window`) should implicitly handle the actual geometry.
    *   **Decision:** Defer this. While potentially useful for diagnostics or very advanced custom blending, it adds complexity to `block_table` for marginal initial gain. The `weight_window` itself is the key for reconstruction.

**Finalizing the "Stand-Alone" Design Aspects:**

*   The core Octree PCA algorithm (`build_octree`, local PCA fitting, coefficient calculation, global PCA on residuals, reconstruction blending) should be implementable as a set of R functions that take data matrices, mask arrays, and a list of parameters (matching the schema) as input, and return lists/arrays/tables of results.
*   LNA's `forward_step.spat.oct_pca` and `invert_step.spat.oct_pca` S3 methods would then be relatively thin wrappers around this core algorithmic toolkit. They would handle:
    *   Getting `dense_mat` and `mask` from `handle`.
    *   Passing `desc$params` to the core algorithm.
    *   Using LNA's HDF5 helpers (`h5_write_dataset`, `h5_create_empty_dataset`, creating groups using the `block_id_hdf5_format_string`) to store the results according to the defined LNA HDF5 layout.
    *   Updating `handle$plan` with dataset definitions.
    *   Managing the `handle$stash` for data flow to/from other LNA transforms.

This separation makes the core Octree PCA logic reusable outside LNA and testable independently.

---

## Definitive Proposal: `spat.oct_pca` - Octree Spatial PCA Transform (Appendix for LNA)

This appendix details the `spat.oct_pca` transform, which applies Principal Component Analysis (PCA) within spatially distinct, hierarchically organized blocks (an octree structure) of the input volume. This allows the transform to adapt to local data complexity, potentially offering better compression than global PCA when signals are spatially heterogeneous.

### 1. High-Level Concept

Traditional PCA applied to fMRI data learns global spatial modes. Octree PCA, by contrast, subdivides the brain volume (defined by a mask) into smaller cubic regions using an octree. PCA is then performed independently within each of these "leaf" cubes (optionally considering an overlap region with neighbors).

*   **Adaptivity:** The octree recursion stops based on criteria like maximum depth, minimum voxels per block, or local explained variance, allowing shallow trees in homogeneous regions and deeper trees in complex areas.
*   **Locality:** Each local PCA captures local spatial "texture" and correlations.
*   **Overlap & Blending:** To mitigate block-edge artifacts during reconstruction, cubes can overlap, and a blending window (e.g., triangular, Hann) is used to smoothly combine contributions from multiple overlapping cubes.
*   **Global Stage (Optional):** A final global PCA can be applied to the residuals after local reconstructions, or to the mean activity of leaf blocks, to capture any remaining large-scale variance. (For v1.0, global PCA on full voxel-space residuals is specified).

### 2. Transform Name and Version

*   `type`: `"spat.oct_pca"`
*   `version`: `"1.0"`

### 3. Schema (`inst/schemas/spat.oct_pca.schema.json`)

```json
{
  "$id": "https://neurocompress.org/schemas/lna/2.0/spat.oct_pca.schema.json",
  "type": "object",
  "title": "Parameters for Octree Spatial PCA (spat.oct_pca) transform",
  "$schema": "http://json-schema.org/draft-07/schema#",
  "properties": {
    "type":              { "const": "spat.oct_pca" },
    "version":           { "const": "1.0" },
    "max_depth":         { 
      "type": "integer", "minimum": 0, "default": 3,
      "description": "Maximum octree subdivision depth (root block is level 0)." 
    },
    "min_block_vox":     { 
      "type": "integer", "minimum": 8, "default": 512,
      "description": "Stop splitting a block if its (interior) voxel count falls below this." 
    },
    "split_var_thresh":  { 
      "type": "number",  "minimum": 0, "maximum": 1.0, "default": 0.01,
      "description": "Stop splitting if local PCA (using up to 'k_local' components) explains >= (1 - split_var_thresh) of variance. Also used for adaptive k_local_actual selection per block." 
    },
    "k_local":           { 
      "type": "integer", "minimum": 1, "default": 8,
      "description": "Maximum number of local PCA components to keep per leaf block. Actual number may be less if adaptive k (via split_var_thresh) is used." 
    },
    "svd_method":        { 
      "enum": ["irlba", "randomized", "svd"], "default": "irlba",
      "description": "SVD method for local PCA (e.g., irlba for speed, svd for base R)." 
    },
    "center":            { 
      "type": "boolean", "default": true,
      "description": "If true, center data (subtract voxel-wise means) within each block before local PCA." 
    },
    "fit_on_halo":       {
      "type": "boolean", "default": true,
      "description": "If true, local PCA for a block is fit on data including its overlap halo. If false, only on interior voxels (overlap still used for coefficient calculation if present)."
    },
    "overlap":           { 
      "type": "integer", "minimum": 0, "default": 2,
      "description": "Number of voxels of half-overlap on each face of a block (e.g., overlap=2 means extend 2 voxels out)." 
    },
    "blend_method":      { 
      "enum": ["window", "iterative", "rect"], "default": "window",
      "description": "'window': use precomputed weight_window. 'iterative': apply iterative refinement. 'rect': simple averaging (boxcar window)."
    },
    "blend_window_type": { 
      "enum": ["triangular", "hann"], "default": "triangular",
      "description": "Type of blending window if blend_method is 'window' or 'iterative'." 
    },
    "iter_blend_passes": { 
      "type": "integer", "minimum": 1, "maximum": 5, "default": 2,
      "description": "Number of passes for iterative blending (if blend_method='iterative')." 
    },
    "k_global":          { 
      "type": "integer", "minimum": 0, "default": 0,
      "description": "If >0, perform a global PCA on full voxel-space residuals after local reconstructions. Stores a global basis and global coefficients." 
    },
    "block_id_hdf5_format_string": {
      "type": "string", "default": "%06d",
      "description": "Printf-style format string for converting integer block_ids to HDF5 group names (e.g., for /basis/oct_pca/blocks/{formatted_block_id}/)."
    },
    "chunk_vox_fitting": { 
      "type": "integer", "minimum": 4096, "default": 32768,
      "description": "Approximate number of voxels (within a block's extended bounds) to process in RAM per SVD minibatch if incremental SVD is needed for very large blocks/time series." 
    },
    "seed":              { 
      "type": "integer", "default": 1,
      "description": "RNG seed for any stochastic parts of SVD (e.g., randomized SVD) to ensure reproducibility."
    }
  },
  "required": ["max_depth", "min_block_vox", "split_var_thresh", "k_local"]
}
```

### 4. HDF5 Layout

*   **/spat/oct_pca/block_table** (HDF5 Compound Dataset or structured array)
    *   Fields:
        *   `block_id`: uint32 (0-based unique ID for this leaf block)
        *   `parent_id`: int32 (-1 for root, else parent's `block_id`)
        *   `level`: uint8 (octree depth, root=0)
        *   `x0, x1, y0, y1, z0, z1`: uint16 (inclusive voxel bounds of the block's *interior/owned* region within the input mask space)
        *   `k_local_actual`: uint8 (actual number of PCA components kept for this block, adaptively chosen)
        *   `coeff_offset`: uint32 (0-based starting column index for this block's coefficients within `/scans/.../coeff_blocks`)
*   **(Optional) /spat/oct_pca/global_basis/**
    *   `matrix`: float32 array [`N_mask_voxels`, `k_global`] (Global PCA basis vectors for residuals)
    *   `mean`: float32 array [`N_mask_voxels`] (Mean of residuals if centered before global PCA)
*   **/basis/oct_pca/blocks/{formatted_block_id}/** (Group for each leaf block, `{formatted_block_id}` from `sprintf(params$block_id_hdf5_format_string, block_id)`)
    *   `matrix`: float32 array [`N_voxels_in_block_extended_bounds`, `k_local_actual`] (Local PCA basis vectors $U_{block}$)
    *   `mean`: float32 array [`N_voxels_in_block_extended_bounds`] (Voxel-wise means if `params$center=TRUE`)
    *   `(Optional) weight_window`: float32 array matching block interior dimensions (e.g., $X_b \times Y_b \times Z_b$) for blending if not 'rect'.
*   **/scans/{run_id}/embedding/oct_pca/**
    *   `(Optional) coeff_global`: float32 array [`TimePoints`, `k_global`]
    *   `coeff_blocks`: float32 array [`TimePoints`, $\sum k_{local\_actual}$] (Concatenated local PCA coefficients)
*   **(Optional) /aux_meta/oct_pca/coding_index** (uint32 array, if needed for complex variable k_local mapping, though `block_table`'s `coeff_offset` and `k_local_actual` might suffice).

### 5. Algorithms

#### 5.1. `forward_step.spat.oct_pca` (Writer-Side)

Input: `X_input` (Voxels_in_Mask $\times$ TimePoints) from `handle$stash$dense_mat`, `mask_array` from `handle$mask_info`.

1.  **Build Octree (`lna:::build_octree_for_pca`)**:
    *   Input: `mask_array`, `params$max_depth`, `params$min_block_vox`, `params$overlap`.
    *   Output: `octree_structure` (list of leaf node definitions: `leaf$id`, `leaf$parent_id`, `leaf$level`, `leaf$interior_bounds_mask_space`, `leaf$extended_bounds_mask_space` (interior + overlap, clipped to main mask)).
    *   Recursive splitting: Start with bounding box of `mask_array`. Split along largest dimension.
    *   Stop criteria for a node:
        *   `current_level == params$max_depth`.
        *   `N_voxels_in_node_interior < params$min_block_vox`.
        *   (If `params$split_var_thresh < 1.0`): Preliminary PCA on node's *interior* data $X_{node\_interior}$ using `params$k_local` (max). If `explained_variance >= (1 - params$split_var_thresh)`, stop.
2.  **Fit Local PCA & Collect Coefficients:**
    *   Initialize `block_table_data = []`, `all_coeff_blocks_list = []`, `current_coeff_offset = 0`.
    *   `set.seed(params$seed)` (for SVD reproducibility if `svd_method` is stochastic).
    *   For each `leaf` in `octree_structure$leaves`:
        a.  Extract $X_{leaf\_extended}$ (data from `leaf$extended_bounds_mask_space`, shape $N_{vox\_extended} \times T$). If `params$fit_on_halo == FALSE`, use $X_{leaf\_interior}$ instead for SVD fitting, but still calculate coefficients for the extended region.
        b.  If $X_{leaf\_extended}$ is too large for memory (`N_vox_extended * T > params$chunk_vox_fitting * T_approx_element_size`), use an incremental SVD approach or process time chunks.
        c.  If `params$center`, compute `mean_vector_leaf = rowMeans(X_{leaf_extended})`. $X_c = X_{leaf_extended} - mean_vector_leaf$. Else $X_c = X_{leaf_extended}$.
        d.  SVD: `svd_result = perform_svd(X_c, k_max = params$k_local + 4, method = params$svd_method)`. (e.g., `irlba::irlba`, `rsvd::rsvd`, or `base::svd`).
        e.  **Adaptive `k_local_actual`:**
            *   `cumulative_energy = cumsum(svd_result$d^2) / sum(svd_result$d^2)`.
            *   `k_star = min(which(cumulative_energy >= (1 - params$split_var_thresh)))`.
            *   `k_local_actual = min(k_star, params$k_local)`.
        f.  Basis $U_{leaf} = svd\_result\$u[, 1:k_{local\_actual}, drop=FALSE]$ ($N_{vox\_extended} \times k_{local\_actual}$).
        g.  Coefficients $C_{leaf} = \text{crossprod}(U_{leaf}, X_c)$ ($k_{local\_actual} \times T$). Transpose to $T \times k_{local\_actual}$ for storage convention.
        h.  Add $C_{leaf}$ (as $T \times k_{local\_actual}$) to `all_coeff_blocks_list`.
        i.  Store $U_{leaf}$ and `mean_vector_leaf` to HDF5: `/basis/oct_pca/blocks/{formatted_leaf_id}/matrix` and `/mean`.
        j.  If `params$blend_method == "window"`, pre-compute and store `weight_window` for `leaf$interior_bounds` (e.g., based on `params$blend_window_type`) to `/basis/oct_pca/blocks/{formatted_leaf_id}/weight_window`.
        k.  Append row to `block_table_data`: `leaf$id`, ..., `leaf$interior_bounds`, `k_local_actual`, `current_coeff_offset`.
        l.  `current_coeff_offset = current_coeff_offset + k_local_actual`.
3.  **Concatenate & Write Coefficients:**
    *   `coeff_blocks_matrix = do.call(cbind, all_coeff_blocks_list)` ($T \times \sum k_{local\_actual}$).
    *   Write `coeff_blocks_matrix` to `/scans/{run_id}/embedding/oct_pca/coeff_blocks` (chunked appropriately).
    *   Write `block_table_data` to `/spat/oct_pca/block_table`.
4.  **Optional Global PCA on Residuals:**
    *   If `params$k_global > 0`:
        a.  Reconstruct $\hat{X}_{block\_recon}$ by applying inverse of local PCAs (using $U_{leaf}, C_{leaf}, mean_{leaf}$) for all *interior* voxels of each leaf and blending/summing where necessary (this is a partial reconstruction only on owned voxels, no complex blending yet).
        b.  Residuals $R = X_{input} - \hat{X}_{block\_recon}$ (Voxels\_in\_Mask $\times T$).
        c.  PCA on $R$: `global_svd = perform_svd(R, k = params$k_global)`.
        d.  Store global basis $U_{global}$ (Voxels\_in\_Mask $\times k_{global}$) and $mean_{global}$ to `/spat/oct_pca/global_basis/`.
        e.  Store global coefficients $C_{global}$ ($T \times k_{global}$) to `/scans/{run_id}/embedding/oct_pca/coeff_global`.
5.  **Update Plan & Stash:**
    *   Add all dataset definitions to `handle$plan`.
    *   `handle$update_stash(keys="dense_mat", new_values=list())` (coefficients are read from HDF5 by next step).

#### 5.2. `invert_step.spat.oct_pca` (Reader-Side)

Input: `handle$subset` (ROI mask, time indices), parameters from `desc$params`.

1.  Load `/spat/oct_pca/block_table`.
2.  Identify `intersecting_blocks` based on `handle$subset$roi_mask` (if any) and `block_table` (interior bounds).
3.  Initialize $X_{reconstructed}$ (e.g., `neuroim2::NeuroVol` or array of zeros, $N_{vox\_roi} \times T_{subset}$) and `sum_weights_array` (for blending, $N_{vox\_roi}$).
4.  For each `block_id` in `intersecting_blocks`:
    a.  Read `k_local_actual` and `coeff_offset` from `block_table` row for `block_id`.
    b.  Load local basis $U_{block}$ from `/basis/oct_pca/blocks/{formatted_block_id}/matrix`.
    c.  Load local `mean_vector_block` from `/basis/oct_pca/blocks/{formatted_block_id}/mean`.
    d.  Load relevant columns $C_{block}$ from `/scans/{run_id}/embedding/oct_pca/coeff_blocks` using `coeff_offset` and `k_local_actual`, subsetting by `handle$subset$time_idx`. ($T_{subset} \times k_{local\_actual}$).
    e.  Reconstruct data for this block's *extended bounds*: $\hat{X}_{block\_extended} = C_{block} \cdot U_{block}^\top + mean\_vector\_block$ ($T_{subset} \times N_{vox\_extended}$).
    f.  Determine voxels `vox_in_roi_and_block_extended` that are both in the user's ROI and this block's extended region.
    g.  Load/compute blending `weights` for `vox_in_roi_and_block_extended` (from stored `weight_window` or generated based on `params$blend_window_type` and distance from block center/edges, considering `params$overlap`).
    h.  Add to global reconstruction: $X_{reconstructed}[\text{vox_indices_in_output}, :] \leftarrow X_{reconstructed}[\text{vox_indices_in_output}, :] + \hat{X}_{block\_extended}[:, \text{vox_indices_in_block_extended}] \cdot \text{weights}$.
    i.  `sum_weights_array[\text{vox_indices_in_output}] \leftarrow sum_weights_array[\text{vox_indices_in_output}] + \text{weights}$.
5.  **Normalize Blended Regions:** $X_{reconstructed} \leftarrow X_{reconstructed} / \text{sum_weights_array}$ (where `sum_weights_array > 0`).
6.  **Iterative Blending (if `params$blend_method == "iterative"`):**
    *   For `p` in `1...params$iter_blend_passes`:
        *   For each `block_id` in `intersecting_blocks`:
            *   Extract current estimate for its extended region $X^{(p-1)}_{block\_extended}$ from $X_{reconstructed}$.
            *   Refine coefficients: $C'_{block} = (X^{(p-1)}_{block\_extended} - mean\_vector\_block) \cdot U_{block}$.
            *   Reconstruct with refined coeffs: $\hat{Y}_{block\_extended} = C'_{block} \cdot U_{block}^\top + mean\_vector\_block$.
            *   Store $\hat{Y}_{block\_extended}$ temporarily.
        *   Re-blend all $\hat{Y}_{block\_extended}$ contributions to get updated $X_{reconstructed}$.
7.  **Add Global Component:** If `params$k_global > 0`:
    *   Load $U_{global}$ and $mean_{global}$. Load $C_{global}$ (subset by time).
    *   $\hat{X}_{global} = C_{global} \cdot U_{global}^\top + mean_{global}$.
    *   $X_{reconstructed} \leftarrow X_{reconstructed} + \hat{X}_{global}$ (subset $U_{global}$ by ROI if needed).
8.  Return $X_{reconstructed}$ (reshaped to 4D array if original input was such) into `handle$stash`.

### 6. Storage Size Estimate & Advantages

(As in previous proposal - still valid, <20MiB raw for example, good gzip potential). Advantages include adaptive rank, overlap blending, and coefficient slab structure.

### 7. Implementation Notes

*   Core algorithmic parts (`build_octree_for_pca`, local SVD loop, blending) should be modular R functions.
*   LNA S3 methods (`forward_step.spat.oct_pca`, `invert_step.spat.oct_pca`) will wrap these, handle `DataHandle`/`Plan` interaction, and use LNA HDF5 helpers (including the new `lna:::h5_create_empty_dataset` for pre-allocation in block-wise writes).
*   `neuroim2` used for mask handling (`igraph::components` can operate on `as.array(mask_neurovol)`), coordinate transformations if needed for block definitions in world space, and potentially for `NeuroVol` representations of intermediate data if beneficial. Input `X` is assumed to be voxels $\times$ time matrix extracted from `NeuroVec`.

This detailed proposal should provide a solid blueprint for implementing a powerful and flexible Octree PCA transform that can operate stand-alone algorithmically while integrating cleanly into LNA.