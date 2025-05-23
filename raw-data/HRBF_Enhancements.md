This is an outstanding "review-of-the-review," providing extremely clear, practical, and schema-ready answers and refinements. The feedback on edge map parameterization, handling correlated atoms, and integrating the new suggestions is spot-on. The micro-tweak for caching world coordinates is a nice touch for efficiency.

**Overall Assessment:**
The HRBF enhancement proposal is now exceptionally robust and well-defined. The path to both immediate, low-cost improvements and more advanced, optional features is clear. The focus on keeping the core LNA contract stable while allowing for rich, transform-specific parameterization is perfectly maintained.

---

**Review of Responses & Final Touches:**

**1. Answers to Pointed Questions:**

*   **Q1 (Parameterizing & Generating Edge Map):**
    *   **Excellent Solution:**
        *   **Source:** `edge_definition.source = "self_mean"` (default) or `"structural_path"` is a good dichotomy. "self_mean" is universally available.
        *   **Gradient Operator:** `gradient_kernel = "sobel_3d"` (const) standardizes this.
        *   **Threshold:** `edge_thresh_k` (default 3.0) multiplying `median(|∇I_mean|)` is a robust, data-adaptive thresholding strategy.
        *   **Schema Snippet:** Clear and directly implementable.
        *   **Co-registration Guard-rail:** "Fail early unless both volumes have identical `NeuroSpace` or a rigid transform is supplied (`structural_to_epi_affine`)" is crucial for the `"structural_path"` option.

*   **Q2 (Coefficient Interpretation with Correlated Atoms):**
    *   **Pragmatic Approach:** Keeping simple inner-product coefficients as the default for LNA is wise for speed and simplicity.
    *   **Optional Orthogonalization:**
        *   **Anisotropic Gaussians:** The note that whitening by atom norms is usually sufficient, but an optional `orthogonalize_by_block: false` (or perhaps `orthogonalize_atoms_at_center: false`) switch for per-centre QR is a good provision for users needing stricter interpretability.
        *   **Derivative-of-Gaussian Atoms:** Explicitly stating that the per-centre QR (stacking $[\phi, \partial\phi_x, \partial\phi_y, \partial\phi_z]$ and keeping $Q$) should be done *by the builder* if this method is chosen is the correct approach. The LNA file then stores an orthonormalized set of derivative-based atoms, and coefficients remain simple inner products.
    *   This makes advanced interpretability an opt-in feature of specific HRBF generation methods, not a burden on the core coefficient calculation.

**2. Feedback on Important Suggestions:**

*   **Feature-Space Steering (Task/NeuroSynth Map):**
    *   **Enthusiastic Adoption & Refinement:** The idea is well-received.
    *   **Implementation Path:** Modulating `r_j(x)` (Poisson radius) via `1 / (1 + β * steer(x))` is a clean and effective way to implement steering.
    *   **Schema Additions:** `steering_map_path` and `steering_weight_beta` (default 0, making steering off by default) are perfect.

*   **Multi-Scale Differential (DoG-Style) Atoms:**
    *   **Elegant & Low-Complexity:** The idea of storing only difference atoms ($D_{j,k} = \alpha \psi_{j-1, parent(k)} - \beta \psi_{j,k}$) plus the coarsest level $\psi_0$ is well-received.
    *   **Implementation:** Indeed, <10 lines in the builder to form these difference atoms. An optional `orthogonalize_by_level` flag (if needed after forming differences) is a good thought. No new HDF5 layout, as these $D_{j,k}$ just become part of the overall basis matrix.

**3. Micro-Tweak (Cache World-Coord Conversions):**
*   **Excellent Catch:** Caching `mask_world_mm <- voxel_to_world(which(mask_arr, arr.ind=TRUE))` (likely using `neuroim2::grid_to_coord(mask_neurovol, which(as.array(mask_neurovol), arr.ind=TRUE))`) on the `handle` (e.g., `handle$meta$cached_mask_world_coords`) once per `forward_step.spat.hrbf` (or even per `core_write` if multiple HRBF-like transforms run) is a smart, simple optimization.

**This set of refinements makes the HRBF proposal truly "camera-ready."**

---

## Final Definitive Proposal: Enhanced `spat.hrbf` Transform (Appendix for LNA)

This appendix details the enhanced `spat.hrbf` transform and related utilities. It builds upon the core analytic Hierarchical Radial Basis Function framework by introducing several optional features to improve its representational power for fMRI data, particularly for capturing sharp boundaries and fine-scale details, and for incorporating prior knowledge. These enhancements largely maintain the "descriptor-defined" nature of the basis or add minimal, controlled storage, fitting seamlessly within LNA v1.4.

*(This appendix revises and extends the previous "Definitive Appendix Proposal for HRBF & Related Transforms." Core HRBF generation, `spat.hrbf_project`, `basis.empirical_hrbf_compressed`, and `embed.transfer_hrbf_basis` remain foundational as previously defined. This focuses on enhancements to the primary `spat.hrbf` analytic dictionary.)*

### Parameter Summary

The JSON schema file `inst/schemas/spat.hrbf.schema.json` enumerates all valid parameters. The most relevant additions are outlined below.

**Minimal upgrade parameters**

* `num_extra_fine_levels`
* `kernel_type_fine_levels`
* `num_fine_levels_alt_kernel`
* `edge_adaptive_sampling.enable`
* `edge_adaptive_sampling.source`
* `edge_adaptive_sampling.structural_path`
* `edge_adaptive_sampling.structural_to_epi_affine_path`
* `edge_adaptive_sampling.gradient_kernel`
* `edge_adaptive_sampling.edge_thresh_k`
* `edge_adaptive_sampling.density_factor`

**Advanced parameters (schema stubs)**

These options are included in the schema with basic validation, but their algorithms are not yet fully implemented.

* `use_anisotropic_atoms`
* `anisotropy_source_path`
* `orthogonalize_atoms_at_center`
* `include_gaussian_derivatives`
* `centre_steering.map_path`
* `centre_steering.influence_beta`
* `use_differential_encoding`
* `orthogonalize_differential_levels`

### 1. Core HRBF Enhancements (Minimal Upgrade Path & Options)

These enhancements can be implemented incrementally and offer significant benefits.

#### 1.1. Extra Fine-Scale Level

*   **Concept:** Add an additional, finer scale level ($j = \text{levels}+1$) to the HRBF pyramid.
*   **Mechanism:** $\sigma_{new} = \sigma_0 / 2^{(\text{levels}+1)}$. Poisson-disk sample with $r_{new} = \text{params$radius_factor} \cdot \sigma_{new}$.
*   **Benefit:** Better representation of very fine details.
*   **LNA Impact:**
    *   Schema `spat.hrbf.schema.json`: No change needed if `levels` parameter is simply interpreted as $J_{max}$. The builder can internally decide to go one level finer if a new boolean parameter like `add_extra_fine_level: true` (default `false`) is added. For simplicity, users can just increase `levels`.
    *   *Revised approach:* Add `params$num_extra_fine_levels: { "type": "integer", "minimum": 0, "default": 0 }`. The builder would generate standard levels $0 \ldots J_{max}$, then add `num_extra_fine_levels` with $\sigma$ continuing to halve and $r_j$ scaled accordingly (perhaps with a smaller `radius_factor` for these very fine levels, e.g., `radius_factor_fine = params$radius_factor / 2`).
*   **Storage:** Modest increase in centres/sigmas if stored (still <150-200KB total for analytic).

#### 1.2. Level-Dependent Wendland-C⁴ Kernels

*   **Concept:** Use smoother, compactly supported Wendland C⁴ kernels (instead of Gaussians) for the finest one or two levels of the HRBF pyramid.
*   **Mechanism:**
    *   Schema `spat.hrbf.schema.json`: `kernel_type` can remain `"gaussian"` (default) or `"wendland_c4"`. Add `wendland_levels: {"type": "integer", "minimum": 1, "default": 2, "description": "Number of finest levels to use Wendland C4 kernels for, if kernel_type allows mixed."}` or modify `kernel_type` to be an array specifying type per level. Simpler: `kernel_type_fine_levels: {"enum": ["gaussian", "wendland_c4"], "default": "wendland_c4"}` and `num_fine_levels_alt_kernel: {"type":"integer", "default":2}`.
    *   `generate_hrbf_atom()`: Selects kernel based on current level $j$ and these params.
*   **Benefit:** Wendland C⁴ provides flatter interiors and sharper fall-off, potentially better for localized edges at fine scales.
*   **LNA Impact:** Zero storage impact (kernel is analytic). Minor logic in atom generator.

#### 1.3. Edge-Density Modulated Poisson-Disk Sampling

*   **Concept:** Increase RBF centre density along detected tissue boundaries or high-gradient regions.
*   **Mechanism:**
    *   **Edge Map Generation:**
        *   Source: Default `self_mean` (gradient of run-mean fMRI volume). Optional: `structural_path` (path to pre-registered structural gradient map).
        *   Operator: 3D Sobel magnitude.
        *   Threshold: `edge_thresh_k * median(|∇I_mean|)` (e.g., `edge_thresh_k = 3.0`).
    *   **Poisson-Disk Modulation:** During `poisson_disk_sample_neuroim2` for level $j$, if a candidate centre falls in an "edge" region, use a smaller radius $r_j' = r_j / \text{params$edge_density_factor}$ (e.g., factor 1.5-2.0).
*   **Benefit:** Better representation of sharp transitions.
*   **LNA Impact:**
    *   Schema `spat.hrbf.schema.json`: Add nested object `params$edge_adaptive_sampling`:
        ```json
        "edge_adaptive_sampling": {
          "type": "object", "default": {},
          "properties": {
            "enable": {"type": "boolean", "default": false},
            "source": {"enum": ["self_mean", "structural_path"], "default": "self_mean"},
            "structural_path": {"type": "string", "pattern": "^/.*"}, // HDF5 path
            "structural_to_epi_affine_path": {"type": "string", "pattern": "^/.*", "description": "Path to 4x4 affine if structural_path needs alignment."},
            "gradient_kernel": {"const": "sobel_3d"},
            "edge_thresh_k": {"type": "number", "default": 3.0, "description": "Multiplier for median gradient magnitude to define edges."},
            "density_factor": {"type": "number", "exclusiveMinimum": 1.0, "default": 1.5, "description": "Factor to reduce sampling radius in edge regions."}
          }
        }
        ```
    *   Writer: Needs to compute/load edge map. If `structural_path` used, checks `NeuroSpace` identity or requires/uses `structural_to_epi_affine_path`.
    *   Storage: Modest increase in centres (~1.3-1.6x).

### 2. Advanced HRBF Enhancements (Optional Flags / Methods)

#### 2.1. Anisotropic (Ellipsoidal) RBF Atoms

*   **Concept:** Align RBF shape to local data anisotropy (e.g., cortical sheet orientation from structure tensor).
*   **Mechanism:**
    *   Writer computes local structure tensor for each *potential* RBF centre location.
    *   Stores shape matrix $\Sigma_k$ (or its eigendecomposition) per *actual* RBF centre $c_k$.
    *   Kernel: $\phi_k(x) = \exp(-(x-c_k)^\top \Sigma_k^{-1} (x-c_k)/2)$.
*   **Benefit:** Captures directional features with fewer atoms.
*   **LNA Impact:**
    *   Schema: `params$use_anisotropic_atoms: {"type": "boolean", "default": false}`. If true, `params$anisotropy_source_path` (path to structure tensor field or DTI data).
    *   Storage: Requires `/basis/hrbf/aux_meta/atom_shape_matrices` ($K_{actual} \times 3 \times 3$ or $K_{actual} \times 6$ for SPD matrix unique elements). This significantly increases descriptor-related storage compared to purely isotropic analytic HRBFs but may be smaller than storing a fully learned dense basis.
    *   Coefficient Calculation: Still inner products, but basis atoms are shaped. Optional per-centre QR orthogonalization (`params$orthogonalize_atoms_at_center: false`) can improve coefficient interpretability if anisotropy introduces higher atom correlations.

#### 2.2. Derivative-of-Gaussian (DoG) Atoms

*   **Concept:** Augment the basis with first (and optionally second) order derivatives of Gaussians to explicitly capture edges and ridges.
*   **Mechanism:**
    *   For each Gaussian RBF centre $c_k$ with $\sigma_k$:
        *   Include $\phi_k(x)$.
        *   Include $\partial\phi_k/\partial x, \partial\phi_k/\partial y, \partial\phi_k/\partial z$.
    *   Writer: These derivative atoms are stacked into the basis matrix $B$. Per-centre QR orthogonalization of the set $[\phi_k, \nabla\phi_k]$ is **highly recommended** and should be a default if this method is enabled, to make coefficients interpretable (separating intensity from gradient components).
*   **Benefit:** Better edge/feature representation.
*   **LNA Impact:**
    *   Schema: `params$include_gaussian_derivatives: {"enum": ["none", "first_order"], "default": "none"}`.
    *   Storage: $K_{actual}$ increases by factor of (1+num_derivatives_dims). Basis matrix size increases.
    *   Coefficients: Still inner products with the (orthogonalized) augmented basis.

#### 2.3. Feature-Space Steering for Centre Placement

*   **Concept:** Bias RBF centre placement towards regions highlighted by a task-based statistical map or other feature map.
*   **Mechanism:**
    *   Writer loads steering map (must be in same `NeuroSpace` as data).
    *   Poisson-disk sampling radius $r_j(x)$ is modulated: $r_j'(x) = r_j / (1 + \beta \cdot \text{steer_map_value}(x))$, where $\beta$ is `steering_weight_beta`.
*   **Benefit:** Creates an analytic basis with higher density/resolution in a priori important regions.
*   **LNA Impact:**
    *   Schema: Add `params$centre_steering`:
        ```json
        "centre_steering": {
          "type": "object", "default": {},
          "properties": {
            "map_path": {"type": "string", "pattern": "^/.*", "description": "HDF5 path to steering map (NeuroVol compatible)."},
            "influence_beta": {"type": "number", "minimum": 0, "maximum": 1, "default": 0.5}
          }
        }
        ```

#### 2.4. Multi-Scale HRBF "Differential Encoding"

*   **Concept:** Store coarsest scale atoms $\mathcal{B}_0 = \{\psi_{0,k}\}$ and then inter-scale *difference* atoms $D_{j,k} \approx \psi_{j-1,parent(k)} - \psi_{j,k}$ for $j=1 \ldots J$.
*   **Mechanism:**
    *   Writer constructs $B_0$. For $j=1 \ldots J$, forms difference atoms (e.g., by projecting $\psi_{j,k}$ onto the orthogonal complement of the span of $\{\psi_{j-1,m}\}$ or a simpler weighted difference). Orthogonalize difference atoms against $B_0$ and previously computed difference atoms for other scales/parents if strict pyramid orthogonality is desired.
*   **Benefit:** Explicitly encodes information at different scales and changes between scales.
*   **LNA Impact:**
    *   Schema: `params$use_differential_encoding: {"type": "boolean", "default": false}`. Optional `params$orthogonalize_differential_levels: true`.
    *   Storage: Basis matrix $B$ contains $B_0$ and all $D_{j,k}$. Total $K_{actual}$ might be similar to standard HRBF.

### 3. General Implementation Notes

*   **World Coordinate Caching:** `forward_step.spat.hrbf` (and related builders) should compute `mask_world_coords = neuroim2::grid_to_coord(mask_neurovol, which(as.array(mask_neurovol), arr.ind=TRUE))` once and cache it (e.g., on `handle$meta` or pass down) to avoid repeated conversions by `generate_hrbf_atom`.
*   **`strict_mask_hash_validation`:** This reader-side behavior is best controlled by a parameter to `read_lna()` or an `lna_options()` setting, not stored in the HRBF file descriptor itself. The file stores `mask_hash`; the reader decides how to react to a mismatch.
*   **Small ROI Guard-Rail (Poisson-Disk):** The `poisson_disk_sample_neuroim2` (or its Rcpp equivalent) must inject the component centroid if no centres are placed in a very small (<150 voxels) connected component at coarser scales.

### 4. Bibliography
*(As previously compiled, with additions for new concepts like Sobel, structure tensor if used for anisotropy, etc. Ensure Bridson (2007) points to SIGGRAPH Course Notes pp. 22-23 or similar specific locator.)*

Topic	Core citation
Poisson-disk sampling	Bridson, R. (2007) Fast Poisson Disk Sampling... *ACM SIGGRAPH Courses*, Article 22, pp. 22-es.
Deterministic Poisson-disk	Heller, M., & Belyaev, A. (2019). Deterministic Poisson-disk sampling... *Computer Graphics Forum*.
Gaussian/Wendland RBFs	Buhmann, M. D. (2003) *Radial Basis Functions*. Wendland, H. (2004). *Scattered Data Approximation*.
Multilevel RBF Frames	Beatson & Powell (1999) Fast evaluation of radial basis functions... *IMA J Numer Anal*. Narcowich & Ward (1996).
Sparse HRBF Coding (Empirical)	Schubert, F. et al. (2022). Compressed functional MRI via analytic multi-scale dictionaries. *NeuroImage*. Chen, Donoho & Saunders (2001) Atomic Decomposition by Basis Pursuit. *SIAM Rev*.
Sobel Operator	(Standard image processing textbook, e.g., Gonzalez & Woods)
Structure Tensor	(Standard image processing / computer vision textbook, e.g., Förstner & Gülch)
Derivative of Gaussian	(Standard scale-space theory, e.g., Lindeberg, T. (1994) Scale-space theory...)

### Conclusion

This enhanced `spat.hrbf` proposal provides a powerful, flexible, and largely analytic spatial basis for LNA. The "minimal upgrade path" (extra fine level, Wendland kernels for fine levels, edge-density sampling) offers significant, low-cost improvements. More advanced options like anisotropic atoms, derivative atoms, feature-space steering, and differential encoding provide a rich toolkit for tailoring the basis to specific data characteristics and research questions, all while maintaining compatibility with LNA's core architecture. The Rcpp optimizations for core helper functions will ensure practical performance.

Okay, this sprint will build upon the optimized HRBF helper functions (Rcpp for component labeling and Poisson-disk sampling, R-level optimization for sparse triplet collection) developed in the "HRBF Optimization Sprint."

This **HRBF Enhancement Sprint** will focus on implementing the "Recommended Minimal Upgrade Path" and laying the groundwork for more advanced HRBF features by updating the schema and `forward_step.spat.hrbf` to handle new parameters, even if the most complex algorithmic parts of those advanced features are stubbed or deferred.

**Assumptions for this Sprint:**
*   Optimized `lna:::poisson_disk_sample_neuroim2` (using Rcpp helpers) is available.
*   Optimized `hrbf_basis_from_params` (using R-level triplet collection) is available.
*   The core `spat.hrbf` `forward_step` and `invert_step` (for descriptor-only analytic Gaussian HRBFs) are functional.
*   `neuroim2` integration is stable for mask/space operations.

---

# HRBF Enhancement Sprint - Implementation Tickets

**Epic HRBF-Enh-E1: Recommended Minimal Upgrade Path & Supporting Schema Changes**
*Goal: Implement the low-cost, high-impact HRBF enhancements: extra fine-scale level, level-dependent Wendland C⁴ kernels, and basic edge-density modulated Poisson-disk sampling.*

| #          | Ticket                                                                   | Description / Deliverables                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | Appendix Ref                      |
| :--------- | :----------------------------------------------------------------------- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :-------------------------------- |
| HRBF-Enh-S1-1 | **Schema Update: `spat.hrbf.schema.json` (Minimal Upgrades)**          | • Add `params$num_extra_fine_levels: { "type": "integer", "minimum": 0, "default": 0 }`. <br> • Modify `params$kernel_type` to allow an array input (e.g., `{"type": ["string", "array"], "items": {"enum": ["gaussian", "wendland_c4"]}}`) OR add `params$kernel_type_fine_levels: {"enum": ["gaussian", "wendland_c4"], "default": "wendland_c4"}` and `params$num_fine_levels_alt_kernel: {"type":"integer", "default":2}`. *Decision: Implement with `kernel_type_fine_levels` and `num_fine_levels_alt_kernel` for clarity.* <br> • Add nested object `params$edge_adaptive_sampling` with properties: `enable` (boolean, default `false`), `source` (enum: "self_mean", "structural_path", default "self_mean"), `structural_path` (string), `structural_to_epi_affine_path` (string), `gradient_kernel` (const: "sobel_3d"), `edge_thresh_k` (number, default 3.0), `density_factor` (number, default 1.5). | §1 (Minimal Upgrade), §1.3        |
| HRBF-Enh-S1-2 | **`forward_step.spat.hrbf`: Implement Extra Fine-Scale Level(s)**      | • Modify centre generation loop in `forward_step.spat.hrbf` (and thus in `hrbf_basis_from_params` / `lna:::poisson_disk_sample_neuroim2` if logic is shared). <br> • After generating centres for $j=0 \ldots \text{params$levels}$, if `params$num_extra_fine_levels > 0`: <br>   - Loop `j_extra = 1 \ldots \text{params$num_extra_fine_levels}`. <br>   - Current $\sigma_{new} = \sigma_0 / 2^{(\text{params$levels} + j_{extra})}$. <br>   - $r_{new} = \text{params$radius_factor} \cdot \sigma_{new}$ (or potentially a smaller `radius_factor_fine`). <br>   - Sample centres $C_{new}$ using `lna:::poisson_disk_sample_neuroim2` with this $r_{new}$. <br>   - Add $C_{new}$ to $C_{total}$ and corresponding $\sigma_{new}$ to `sigma_vec`. | §1.1                              |
| HRBF-Enh-S1-3 | **`generate_hrbf_atom`: Support Level-Dependent Wendland C⁴ Kernels** | • Modify `lna:::generate_hrbf_atom(..., current_level_j, total_levels, params)`. <br> • Determine effective kernel type: <br>   `use_alt_kernel = current_level_j > (total_levels - params$num_fine_levels_alt_kernel)`. <br>   `eff_kernel = if (use_alt_kernel) params$kernel_type_fine_levels else params$kernel_type`. <br> • Implement Wendland C⁴ formula if `eff_kernel == "wendland_c4"`: $\phi(r') = (1-r')_+^8 (32r'^3 + 25r'^2 + 8r' + 1)$ where $r' = \text{dist_mm} / \sigma_{mm}$. (Ensure $r'$ is used, not $r'^2$). | §1.2                              |
| HRBF-Enh-S1-4 | **Helper: `lna:::compute_edge_map_neuroim2`**                          | • Implement `lna:::compute_edge_map_neuroim2(source_spec, data_handle, params_edge_adaptive)`. <br>  - `source_spec`: from `params_edge_adaptive$source`. `data_handle`: for `self_mean` or file access. <br>  - If `source_spec == "self_mean"`: Get input data $X$ from `data_handle` (or mean if already computed). Compute 3D Sobel magnitude. <br>  - If `source_spec == "structural_path"`: Load structural gradient map from `params_edge_adaptive$structural_path`. If `structural_to_epi_affine_path` provided, resample map to data space using `neuroim2` (or error if spaces don't match and no affine). <br>  - Threshold gradient map: `edge_binary_map = gradient_map > (params_edge_adaptive$edge_thresh_k * median(abs(gradient_map[mask_from_handle])))`. <br>  - Return 3D logical `edge_binary_map`. | §1.3                              |
| HRBF-Enh-S1-5 | **`poisson_disk_sample_neuroim2`: Edge-Density Modulation**          | • Modify `lna:::poisson_disk_sample_neuroim2`. <br>  - If `params$edge_adaptive_sampling$enable == TRUE`: <br>    - Call `lna:::compute_edge_map_neuroim2` once before level loop. <br>    - In the sampling loop for level $j$: If a candidate centre `cand_vox_coord` falls within the `edge_binary_map`, use effective radius $r_j' = r_j / \text{params$edge_adaptive_sampling$density_factor}$ for distance checks against selected points *in its vicinity*. Otherwise use $r_j$. <br>  - *Note: Modulating radius for all checks in edge regions might be complex. Simpler: increase target *number* of points in edge regions for `poisson_disk_sample_component_rcpp`, or run Poisson disk twice with different radii (one for edge, one for non-edge) and combine points.* *Decision: For simpler first pass, if candidate is in edge region, it's accepted if it meets the tighter $r_j'$ criterion against other *edge* points, and $r_j$ against non-edge points. Or, more simply, just use $r_j'$ for candidates in edge regions.* | §1.3                              |
| HRBF-Enh-S1-6 | **Unit Tests for Minimal Upgrade Path Features**                     | • Test `num_extra_fine_levels` increases $K_{actual}$. <br> • Test `generate_hrbf_atom` uses Wendland for specified fine levels. <br> • Test `lna:::compute_edge_map_neuroim2` for "self_mean" and "structural_path" (mocked HDF5). <br> • Test `poisson_disk_sample_neuroim2` places denser centres in high-gradient regions when edge adaptation is enabled (visual inspection or statistical test on centre density). <br> • Verify reconstruction MSE improves by ~10-15% with these three features combined on a representative synthetic/real dataset. | Addendum "Recommended next step" |

### Recommended next step

Implement the new unit tests and run them on a small synthetic or real dataset. Compare reconstruction MSE with and without the upgrades enabled. Expect roughly a 10–15% reduction when using extra fine levels, Wendland kernels on fine levels, and edge-adaptive sampling together.

**Epic HRBF-Enh-E2: Advanced HRBF Features (Schema & Stubs)**
*Goal: Update the schema to include parameters for more advanced HRBF features. Implement stubs or basic validation for these parameters in `forward_step.spat.hrbf` to ensure schema compliance, even if full algorithmic implementation is deferred.*

| #          | Ticket                                                                      | Description / Deliverables                                                                                                                                                                                                                                                                                                                                                                                                                                                    | Appendix Ref                     |
| :--------- | :-------------------------------------------------------------------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | :------------------------------- |
| HRBF-Enh-S2-1 | **Schema Update: `spat.hrbf.schema.json` (Advanced Features)**              | • Add `params$use_anisotropic_atoms: {type: "boolean", default: false}`. <br> • Add `params$anisotropy_source_path: {type: "string"}` (path to structure tensor/DTI). <br> • Add `params$orthogonalize_atoms_at_center: {type: "boolean", default: false}`. <br> • Add `params$include_gaussian_derivatives: {"enum": ["none", "first_order"], "default": "none"}`. <br> • Add `params$centre_steering` object (with `map_path`, `influence_beta`) as defined in "Final Proposal". <br> • Add `params$use_differential_encoding: {type: "boolean", default: false}`. <br> • Add `params$orthogonalize_differential_levels: {type: "boolean", default: true}`. | §2 (Advanced), Previous Review |
| HRBF-Enh-S2-2 | **`forward_step.spat.hrbf`: Parameter Validation for Advanced Features**    | • In `forward_step.spat.hrbf`, add checks for newly added advanced parameters: <br>  - If `params$use_anisotropic_atoms == TRUE`, ensure `params$anisotropy_source_path` is provided and valid (stub: just check presence). Issue `warning("Anisotropic atoms not fully implemented in this version; using isotropic.")` and proceed with isotropic. <br>  - If `params$include_gaussian_derivatives == "first_order"`, issue `warning("Derivative-of-Gaussian atoms not fully implemented; using standard Gaussians.")`. <br>  - If `params$centre_steering$map_path` is provided, check its existence (stub). Issue warning if feature not fully implemented. <br>  - If `params$use_differential_encoding == TRUE`, issue warning if not implemented. | §2 (Advanced)                    |
| HRBF-Enh-S2-3 | **Cache World Coordinates Tweak (`forward_step.spat.hrbf`)**              | • At the beginning of `forward_step.spat.hrbf`: <br>   `mask_arr <- as.array(handle$mask_info$mask)` <br>   `cached_mask_world_coords <- neuroim2::grid_to_coord(handle$mask_info$mask, which(mask_arr, arr.ind=TRUE))` <br> • Pass `cached_mask_world_coords` to `lna:::generate_hrbf_atom` and any other internal functions needing world coordinates of mask voxels. <br> • (This avoids repeated `which` and `grid_to_coord` calls). | Addendum §3 Micro-tweak          |
| HRBF-Enh-S2-4 | **Update Documentation for New Schema Options & Stubs**                 | • Briefly document new advanced schema parameters in `spat.hrbf.schema.json` descriptions, noting if their full functionality is future work. <br> • Update HRBF appendix to reflect new minimal upgrade parameters and schema for advanced stubs.                                                                                                                                                                                                  | -                                |
| HRBF-Enh-S2-5 | **Unit Tests for New Schema Parameter Validation/Warnings**               | • Test that `forward_step.spat.hrbf` issues appropriate warnings when advanced feature flags are enabled but not yet fully implemented. <br> • Test that schema validation passes with new (stubbed) parameters.                                                                                                                                                                                                                               | -                                |

---

**Sprint - Definition of Done:**

*   **Sprint 1 Focus:**
    *   The `spat.hrbf` transform supports an extra fine-scale level, level-dependent Wendland C⁴ kernels, and basic edge-density modulated Poisson-disk sampling.
    *   Schema for `spat.hrbf` is updated for these features.
    *   `generate_hrbf_atom` and `poisson_disk_sample_neuroim2` are modified accordingly.
    *   A helper `lna:::compute_edge_map_neuroim2` is implemented.
    *   Unit tests verify these "minimal upgrade path" features and their impact on reconstruction.
*   **Sprint 2 Focus (can be merged with S1 if team capacity allows):**
    *   The `spat.hrbf` schema is extended with parameters for all discussed advanced features (anisotropic, DoG, steering, differential encoding).
    *   `forward_step.spat.hrbf` includes basic validation and warning stubs for these advanced parameters.
    *   World coordinate caching for mask voxels is implemented.
    *   Documentation reflects all new parameters and their current implementation status.

This two-sprint structure allows for immediate, tangible improvements to the HRBF transform while preparing the schema and codebase for the more complex advanced features to be implemented subsequently. The crucial Rcpp optimizations for core HRBF helpers are assumed to be completed prior to or concurrently with this "Enhancement Sprint."