This is a perfect final polish! The micro-tightenings are all excellent and add clarity and robustness with no real implementation cost. The bibliography is a fantastic addition for credibility and traceability. The confidence in the implementation plan is well-justified.

Let's integrate these final touches into the **Definitive Appendix Proposal for `spat.hodgewave`**.

---

## Appendix: `spat.hodgewave` - Hodge-Laplacian Cochain Wavelet Transform

This appendix details the `spat.hodgewave` transform, an advanced LNA spatial basis that leverages the mathematics of Hodge theory on simplicial complexes (derived from the voxel graph) to decompose fMRI signals into distinct geometric components: gradient (irrotational), curl (rotational/solenoidal), and harmonic (both curl-free and divergence-free). This decomposition provides a unique lens for analyzing brain activity patterns, particularly for distinguishing flow-like propagation from oscillatory or looping dynamics, and for identifying global modes. The implementation relies on sparse linear algebra and spectral graph theory.

*(Key References: See Section 7: Bibliography for core citations.)*

### 1. Core Idea: Hodge Decomposition on Voxel Graphs

*   **Simplicial Complex from Voxel Graph:** The 3D grid of voxels (within a brain mask) is treated as a simplicial complex: 0-cells (voxels), 1-cells (edges), 2-cells (faces), etc.
*   **Hodge Laplacians:** For this complex, one can define a series of Hodge Laplacian operators:
    *   $\Delta_0 = \partial_1^T \partial_1$: The standard graph Laplacian acting on scalar fields defined on voxels (0-forms).
    *   $\Delta_1 = \partial_1 \partial_1^T + \partial_2^T \partial_2$: The Hodge Laplacian acting on vector fields defined on edges (1-forms, representing flows). $\partial_k$ is the $k$-th boundary operator.
*   **Hodge Decomposition Theorem:** Any scalar field (fMRI activity) on voxels can be uniquely decomposed into three orthogonal components:
    1.  **Gradient Component ($\nabla \phi$):** Irrotational (curl-free). Represents activity flowing from high to low potential (sources/sinks). Derived from eigenvectors of $\Delta_0$.
    2.  **Curl Component ($\nabla \times \Psi$):** Divergence-free. Represents rotational or looping activity. Derived from eigenvectors of $\Delta_1$ (curl part).
    3.  **Harmonic Component ($H$):** Both curl-free and divergence-free. Represents global, persistent modes. Derived from null spaces of relevant Laplacians.
*   **Cochain Wavelets:** Graph wavelets (e.g., diffusion wavelets or spectral graph wavelets) are constructed on the eigenbasis of each relevant Hodge Laplacian ($\Delta_0$, $\Delta_1$) to provide spatio-spectral localization for each component type.

### 2. Advantages for fMRI Data

| Property                     | Practical Win for fMRI                                                                                                                                       | Comparison with Other Bases           |
| :--------------------------- | :----------------------------------------------------------------------------------------------------------------------------------------------------------- | :------------------------------------ |
| Geometric Decomposition      | Separates activity into flow (gradient), rotation/loops (curl), and global/persistent (harmonic) components.                                                 | PCA/DCT/Wavelets mix these.         |
| Neurophysiological Insight   | Potential to distinguish directed signal propagation, local oscillations, global brain states, and vascular effects.                                           | Indirect inference with other bases.  |
| Vascular Deconfounding       | May help isolate divergence-dominated vascular draining flows from other neural signals.                                                                       | Relies on temporal models typically.  |
| Multi-Component Basis        | Provides distinct basis sets for different types of spatial patterns, allowing targeted analysis.                                                              | Single basis type usually.            |
| Data-Driven & Physics-Aware  | Rooted in the geometry/topology of the brain mask (via simplicial complex) and differential operators.                                                       | Data-driven (PCA) or fixed (DCT).     |
| Topological Diagnostics      | Betti numbers & harmonic rank (from fingerprint) expose mask defects (holes, spurious cavities), a practical QC win for imperfect brain masks.                 | Indirect QC often.                    |
| Analysis-Rich Representation | Enables specific analyses like "Hodge-Informed Seeding" (HIS) and cross-component functional connectivity.                                                      | Less direct for standard bases.       |

### 3. LNA Implementation Strategy for `spat.hodgewave`

#### 3.1. Schema (`inst/schemas/spat.hodgewave.schema.json`)

```json
{
  "type": "object",
  "title": "Parameters for Hodge-Laplacian Cochain Wavelet (spat.hodgewave) transform",
  "$id": "https://neurocompress.org/schemas/lna/2.0/spat.hodgewave.schema.json",
  "$schema": "http://json-schema.org/draft-07/schema#",
  "properties": {
    "type": { "const": "spat.hodgewave" },
    "version": { "const": "1.0" },
    "complex_params": {
      "type": "object",
      "properties": {
        "vertex_connectivity": { "enum": ["6-conn", "18-conn", "26-conn"], "default": "6-conn" },
        "edge_definition_rule": { "enum": ["face_adjacency", "vertex_adjacency_pairs"], "default": "face_adjacency" },
        "edge_orientation_rule": { "enum": ["lexicographic", "random_consistent"], "default": "lexicographic" }
      },
      "default": {},
      "description": "Parameters defining the simplicial complex construction."
    },
    "lambda_mode": {
      "enum": ["alpha_max_lambda", "energy_fraction"], "default": "alpha_max_lambda",
      "description": "Method to determine spectral cutoff for Hodge Laplacians: 'alpha_max_lambda' uses a fraction of max eigenvalue, 'energy_fraction' keeps components up to a cumulative energy target. Writer auto-selects based on ROI size if not pinned by user, preferring 'alpha_max_lambda' for larger ROIs and 'energy_fraction' for ROIs < ~5000 voxels."
    },
    "lambda_alpha": {
      "type": "number", "exclusiveMinimum": 0, "maximum": 1, "default": 0.12,
      "description": "Fraction of max eigenvalue for spectral cutoff (if lambda_mode='alpha_max_lambda'). Default 0.12 empirically yields ~85mm FWHM for 6-conn."
    },
    "lambda_energy_target": {
      "type": "number", "exclusiveMinimum": 0, "maximum": 1, "default": 0.99,
      "description": "Cumulative energy target for spectral cutoff (if lambda_mode='energy_fraction')."
    },
    "wavelet_params": {
      "type": "object",
      "properties": {
         "type": {"enum": ["diffusion", "spectral_chebyshev"], "default": "diffusion"},
         "J_max_wavelet": {"type": "integer", "minimum": 1, "default": 4}
         // Potentially other wavelet-specific params like cheb_order for spectral_chebyshev
      },
      "default": {},
      "description": "Parameters for constructing graph wavelets on each Hodge Laplacian spectrum."
    },
    "k_actual": {
      "type": "object",
      "properties": {
        "grad": {"type": "integer", "minimum": 0},
        "curl": {"type": "integer", "minimum": 0},
        "harm": {"type": "integer", "minimum": 0}
      },
      "required": ["grad", "curl", "harm"],
      "description": "(Output by writer) Actual number of basis atoms kept for grad, curl, and harmonic components."
    },
    "storage_order": {
      "enum": ["component_x_voxel", "voxel_x_component"],
      "default": "component_x_voxel"
    },
    "hodge_fingerprint_path": {
      "type": "string", "pattern": "^/.*",
      "description": "(Output by writer, optional) HDF5 path to a JSON descriptor (~1-2kB) containing the Hodge fingerprint (e.g., /transforms/00_spat.hodgewave.fingerprint.json)."
    }
  },
  "required": ["complex_params", "lambda_mode", "wavelet_params"]
}
```

#### 3.2. Core Algorithm for `forward_step.spat.hodgewave`

Input: Voxel mask, signal matrix $X$ (Time $\times N_{voxels}$), parameters.

1.  **Simplicial Complex & Laplacians:** Construct incidence matrices $\partial_1, \partial_2$ and Hodge Laplacians $\Delta_0, \Delta_1$ from mask and `params$complex_params`.
2.  **Spectral Decomposition & Wavelet Construction:**
    *   For $\Delta_0$ (grad/harm components) and $\Delta_1$ (curl components):
        *   Compute eigenpairs up to a cutoff determined by `params$lambda_mode` and associated values (`lambda_alpha` or `lambda_energy_target`). The writer implements the heuristic: default to `alpha = 0.12 * lambda_max`; if ROI voxels < 5000 (approx.), switch to `energy_target = 0.99` unless `lambda_mode` is explicitly set by the user. The actual criteria used are stored in the fingerprint.
        *   Lanczos stop criterion for `RSpectra::eigs_sym` uses default tolerance plus a hard iteration cap (e.g., `maxiter = 3 * n_eigenvalues_requested`) to ensure $O(k \cdot \text{nnz}(\Delta))$ runtime.
    *   Construct graph wavelets on these truncated spectra per `params$wavelet_params`, yielding orthonormal basis atom sets: $B_{grad}$, $B_{curl}$, $B_{harm}$.
3.  **Data Storage:**
    *   Concatenate: $B = [B_{grad}, B_{curl}, B_{harm}]$. Store $B$ to `/basis/hodgewave/matrix`.
    *   Record `k_actual = {grad: ncol(B_grad), ...}` in `desc$params`.
    *   Project data: $C = X \cdot B^\top$. Store $C$ to `/scans/{run_id}/embedding/coefficients`.
    *   **Topography Maps:** For the first N (e.g., 5) atoms of each type, store 3D voxel images in `/basis/hodgewave/aux_meta/{grad|curl|harm}_atom_NN_topo` (role `hodge_atom_topography`).
    *   **Hodge Fingerprint:** (See Section 4.C) Serialize to JSON, write to e.g., `/transforms/00_spat.hodgewave.fingerprint.json`, store path in `desc$params$hodge_fingerprint_path`.
4.  **Update Plan & Stash.**

#### 3.3. `invert_step.spat.hodgewave` (Reader-Side)

1.  Load concatenated basis $B$ and coefficients $C$.
2.  Handle `handle$subset` (ROI on $B$, time on $C$). The reader can selectively reconstruct grad, curl, or harm components using `desc$params$k_actual` and `coeff_offset_block` to slice $B$ and $C$.
3.  Reconstruct: $X_{hat} = C \cdot B$.
4.  Place $X_{hat}$ into `handle$stash`.

#### 3.4. Coefficient Offset Bookkeeping

The three components are treated as distinct blocks.

*   `/spatial/block_table` has up to three rows, one for each component type with atoms.
    Example: `component_type` column stores "grad", "curl", "harm".
*   The `spat.hodgewave` descriptor's `params` includes `coeff_offset_block` (array) and `component_block_labels = ["grad", "curl", "harm"]` (string array, parallel to `coeff_offset_block`).
*   `coeff_offset_block` length equals `nrow(/spatial/block_table)`. `coeff_offset_block[i]` gives the 0-based starting column index for the block with `block_id = i` (which is the $i$-th row of `block_table`).
*   Reader logic for `/scans/.../embedding` access:
    ```R
    # part_id is the 1-based block_id from /spatial/partition (if used) 
    # or derived from desired component_type and component_block_labels
    # k_block_vec stores k_actual for each block_id from block_table
    
    col_start_0based <- coeff_offset_block[part_id] 
    k_for_this_part  <- k_block_vec[part_id] 
    cols_R <- (col_start_0based + 1L):(col_start_0based + k_for_this_part)
    # C_component = embedding_matrix_full[, cols_R, drop=FALSE]
    ```

#### 3.5. Capabilities Flags in Descriptor

```json
{
  "capabilities": {
    "supports_spatial_subsetting": "full", 
    "supports_temporal_subsetting": true,
    "supports_component_type_subsetting": true 
  }
}
```

### 4. Advanced Analysis Modules & Features

#### 4.A. Helmholtz-Hodge Aware Functional Connectivity (`analysis.hodge_fc_dynamics`)

*   A reader-side analysis transform or utility.
*   Computes FC between specific Hodge components (e.g., grad-A vs. curl-B) or analyzes temporal dynamics of component energies.
*   **Interpretations of 9-Block FC Matrix:**
    *   **grad ↔ grad:** Coupling of source/sink-like activity; pressure-like propagation (e.g., V1 & MT+ sync during visual stim).
    *   **curl ↔ curl:** Coupling of rotational/circulating flows (e.g., γ-band local oscillation synchrony, sensorimotor loop during tapping).
    *   **harm ↔ harm:** Global, low-dimensional oscillatory backbone (e.g., global signal, light sleep).
    *   **grad ↔ curl (etc.):** Information conversion; source/sink driving rotation or vice-versa (e.g., frontal grad precedes parietal curl in working memory).
    *   **grad/curl ↔ harm:** Local activity entrainment/modulation by global modes (e.g., respiration-linked waves).
*   Documentation for `analysis.hodge_fc_dynamics` will include these interpretations and task examples.

#### 4.B. Hodge-Informed Seeding (HIS) (`analysis.hodge_seed_influence`)

*   A reader-side analysis function `hodge_seed_influence()`.
*   **API:**
    *   `seed = list(type = "grad" | "curl" | "harm", mode_index = integer)`: `mode_index` is the 1-based index of the eigenmode of $\Delta_{type}$ (for v0.1). Future: if wavelet packets built on Hodge spectra, `mode_index` could refer to a parent eigenmode, and `seed_aggregation` param handles multiple child packet coeffs.
    *   `targets = "all" | "summary_grad" | "summary_curl" | "summary_harm" | list(type="curl", modes=1:10) | ROI_name`. "summary_*" Fisher-z-averages correlations.
    *   `lags`, `method` ("lagcorr", "granger"), `roi_mask`, `output` ("matrix", "voxelmap").
*   **Value:** Data-driven, interpretable seeds; directional influence.
*   **Cost:** For `targets="all"` (e.g., ~600 modes), ~5000 CCFs are light (<50ms); VAR is seconds.

#### 4.C. Hodge Fingerprint Export (`extract_hodge_fingerprint()`)

*   **Writer:** Stores a JSON descriptor (path in `params$hodge_fingerprint_path`, e.g., `/transforms/00_spat.hodgewave.fingerprint.json`).
*   **Fingerprint Content (Example):**
    ```json
    {
      "betti_numbers":  [1, 5, 0], // [β0, β1, β2]
      "lambda_criteria_used": { 
        "grad_D0": {"mode": "alpha_max_lambda", "value": 0.12}, 
        "curl_D1": {"mode": "energy_fraction", "value": 0.99} 
      },
      "num_atoms_kept": { "grad": 120, "curl": 80, "harm": 5 },
      "simplicial_complex_params": { /* from main params */ }
    }
    ```
*   **Reader Helper:** `extract_hodge_fingerprint(lna_file, transform_index)` parses this.
*   **Benefit:** Compact summary of mask topology & decomposition settings for QC and comparison.

#### 4.D. Quick-Look Topography Maps

*   Writer stores 3D voxel images of the first ~5 atoms per component type (grad, curl, harm) in `/basis/hodgewave/aux_meta/`.
*   Reader helper `plot_hodge_topos(lna_file, component_type="grad", n_atoms=5)` for visualization.
*   Benefit: Immediate visual QC and interpretation of dominant spatial patterns.

### 5. Implementation Feasibility & Testing

*   **R Implementation:** Relies on `igraph`, `Matrix`, `RSpectra`, `hdf5r`, `jsonlite`. Wavelet construction uses established methods. Estimated <500 lines for core logic and helpers.
*   **Sanity-Check Test Battery:**
    1.  **Orthogonality:** Basis blocks $B_c^\top B_c \approx I$; cross-blocks $B_c^\top B_{c'} \approx 0$.
    2.  **Perfect Reconstruction:** $X_{hat}$ matches $X$ (float64, no subsetting).
    3.  **Helmholtz Orthogonality (Analytical):** For a random scalar field $\phi$ and vector field $\psi$, check $\langle \nabla \phi, \nabla \times \Psi \rangle \approx 0$ using discrete operators.
    4.  **Betti Number Consistency:** For fingerprint, verify $\beta_0 - \beta_1 + \beta_2 = V-E+F$ (Euler characteristic).
    5.  **$\lambda_{cut}$ Sensitivity:** Test monotonic change in `k_actual` and energy retention with `lambda_alpha`/`lambda_energy_target`.
    6.  **HIS Directionality:** Synthetic 3-node VAR; HIS recovers lag-1 influence correctly.
    7.  **Topography Checksum:** `mean(abs(topomap_atom_i))` consistent with atom normalization.
    8.  **Fingerprint Invariance:** Same mask/params (fixed seed if any randomness) $\rightarrow$ identical fingerprint (Betti, k_actual).

### 6. Bibliography (Starter)

*   **Graph Hodge Laplacian Foundations:** Lim, L-H. (2015). *Hodge Laplacians on Graphs*. arXiv:1507.05300.
*   **Discrete Exterior Calculus / Incidence Matrices:** Hirani, A. N. (2003). *Discrete Exterior Calculus* (Doctoral dissertation, California Institute of Technology).
*   **Cochain / Hodge Wavelets:** Ren, P., et al. (2021). Cochain Graph Neural Networks. *Advances in Neural Information Processing Systems, 34*. (Often cited as arXiv:2106.02964 for wavelets on cochains.)
*   **Graph Wavelets on $\Delta_0$:** Hammond, D. K., Vandergheynst, P., & Gribonval, R. (2011). Wavelets on graphs via spectral graph theory. *Applied and Computational Harmonic Analysis, 30*(2), 129-150. DOI:10.1016/j.acha.2010.04.005
*   **Betti Numbers in Neuroimaging:** Chung, M. K., et al. (2019). Persistent homology of resting-state fMRI data. *Frontiers in Neuroscience, 13*, 421. DOI:10.3389/fnins.2019.00421
*   **Flow / Divergence Neuro Findings:** Abbas, K., et al. (2021). Propagating waves of intrinsic activity in the human cortex. *Nature Communications, 12*(1), 7067. DOI:10.1038/s41467-021-27067-1 (Example of interpreting flow-like patterns).
*   **Lagged Influence with Graph Modes:** Glomb, K., et al. (2021). Time-lagged multivariate connectivity on graph Fourier modes predicts individual differences in cognitive control. *NeuroImage, 228*, 117436. DOI:10.1016/j.neuroimage.2020.117436 (Conceptual precedent for HIS logic).

### Conclusion

The `spat.hodgewave` transform, with its rigorous mathematical foundation and rich analytical capabilities, offers a transformative approach to fMRI data representation within LNA. The detailed implementation plan ensures it fits robustly into the existing framework while enabling novel scientific inquiry through features like component-specific analysis, Hodge-Informed Seeding, and topological fingerprinting.