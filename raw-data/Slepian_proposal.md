Okay, this is excellent, actionable feedback! The "patch-style summary" is particularly helpful. I'll incorporate these changes to produce a revised Appendix for the LNA package documentation.

The core of the proposal remains sound and well-integrated. The tweaks focus on clarity, schema robustness, and fully leveraging the capabilities.

---

## Appendix: `spat.slepian` - Slepian Function Spatial Basis Transform (Revised)

This appendix details the `spat.slepian` transform, a sophisticated spatial basis for LNA that utilizes Discrete Prolate Spheroidal Sequences (DPSS), commonly known as Slepian functions. It offers excellent energy concentration within a chosen spatial region (e.g., whole brain mask, anatomical ROI, or an octree block) while maintaining smoothness and near-orthogonality, making it a powerful tool for fMRI data compression and analysis.

### Core Idea & Properties

Slepian functions are a set of bandlimited functions that are maximally concentrated in a specific region (spatial or temporal). For fMRI volumes, they offer several advantages:

1.  **Energy Concentration:** Each Slepian eigenvector is optimally concentrated within a chosen spatial support (e.g., a brain mask or a sub-block). This means more signal variance is packed into the leading components compared to standard DCT or unlocalized Fourier modes, leading to efficient compression. *Bandwidth W is expressed in voxel indices (half the number of retained Fourier bins per axis); for example, at 2mm isotropic resolution, a W of 1 voxel corresponds to a spatial frequency cutoff of 0.25 mm⁻¹.*
2.  **Smoothness and Locality:** Like low-order spherical harmonics, Slepians are inherently smooth, which is beneficial for quantization (reducing ringing artifacts). Crucially, they decay rapidly outside their defined support region, allowing a single Slepian basis (or a set of block-wise bases) to effectively serve both whole-brain analyses and fine-grained investigations within specific ROIs.
3.  **Orthonormality within Support:** The discrete Slepian functions are orthogonal over the set of voxels defining their spatial support. This simplifies the forward projection (encoding) and inverse reconstruction (decoding) to single matrix multiplications, fitting perfectly into the LNA `basis` and `embed` transform contract.
4.  **Tunable "Shannon Number" & Concentration:** The discrete Shannon number ($N_d = 2 W_x W_y W_z \cdot V_{vox}$, where $W$'s are half-bandwidths in samples per dimension and $V_{vox}$ is the number of voxels in the support) provides a theoretical guideline for the number of Slepian components needed to capture a very high percentage (e.g., ≥98%) of the signal energy within the support. This allows for principled selection of the number of components (`k`) to keep.
5.  **Eigenvalues as Concentration Metrics:** The eigenvalues ($\lambda_i$) associated with each Slepian component directly quantify the fraction of that component's energy residing *inside* the target support region Ω. This property is unique among common basis sets and enables advanced analytical capabilities (see Section C below).

### LNA Implementation Strategy for `spat.slepian`

#### 1. Transform Definition (`spat.slepian.schema.json`)

A new transform type, `spat.slepian`, is introduced. Its schema allows for flexible configuration:

```json
{
  "type": "object",
  "title": "Parameters for spat.slepian transform",
  "$id": "...", "$schema": "...",
  "properties": {
    "type": { "const": "spat.slepian" },
    "version": { "const": "1.0" },
    "method": {
      "enum": ["global", "global_plus_roi_residual", "block_octree"],
      "default": "global",
      "description": "Slepian fitting strategy: 'global' for whole mask/input, 'global_plus_roi_residual' for global + ROI-specific, 'block_octree' for per-octree-block."
    },
    "region_params": {
      "type": "object",
      "properties": {
        "bandwidth": { "type": "number", "exclusiveMinimum": 0, "description": "Discrete half-bandwidth W (voxel indices, e.g., samples per dimension)." },
        "k": { "type": "integer", "minimum": 1, "description": "Number of Slepian components to keep." },
        "target_energy": { "type": "number", "exclusiveMinimum": 0, "maximum": 1, "description": "Target energy fraction to retain (alternative to k)." },
        "bandwidth_selection": {
          "enum": ["fixed", "auto_energy"], "default": "fixed",
          "description": "'fixed': use specified bandwidth. 'auto_energy': select bandwidth to meet target_energy with up to k components."
        },
        "solver_params": {
          "type": "object",
          "properties": { "max_iter": {"type": "integer"}, "ncv": {"type": "integer"} },
          "additionalProperties": true, "maxProperties": 6,
          "description": "Parameters for RSpectra::eigs_sym (e.g., tol, max_iter, ncv)."
        }
      },
      "oneOf": [
        { "properties": { "bandwidth_selection": { "const": "fixed" } },  "required": ["bandwidth", "k"] },
        { "properties": { "bandwidth_selection": { "const": "auto_energy" } }, "required": ["target_energy"] }
      ]
    },
    "roi_residual_params": {
      "type": "object",
      "properties": {
        "mask_path": { "type": "string", "description": "HDF5 path to the ROI mask dataset." },
        "bandwidth": { "type": "number", "exclusiveMinimum": 0, "description": "Bandwidth W for ROI residual." },
        "k": { "type": "integer", "minimum": 1, "description": "Components for ROI residual." },
        "target_energy": { "type": "number", "exclusiveMinimum": 0, "maximum": 1 },
        "bandwidth_selection": { "enum": ["fixed", "auto_energy"], "default": "fixed" }
        // Similar oneOf for ROI bandwidth/k/target_energy
      },
      "required": ["mask_path"] // Plus conditional requirements based on bandwidth_selection
    },
    "storage_order": {
      "enum": ["component_x_voxel", "voxel_x_component"],
      "default": "component_x_voxel"
    },
    "eigenvalues": { "type": "array", "items": {"type": "number"}, "maxItems": 4096, "description": "(Output by writer) Eigenvalues of computed Slepian components." },
    "k_actual": { "type": "integer", "description": "(Output by writer) Actual number of components kept for the primary region." },
    "k_roi_actual": { "type": "integer", "description": "(Output by writer) Actual number of components kept for the ROI residual (if applicable)." },
    "bandwidth_actual": { "type": "number", "description": "(Output by writer) Actual bandwidth used if auto-selected for primary region." },
    "bandwidth_roi_actual": { "type": "number", "description": "(Output by writer) Actual bandwidth used if auto-selected for ROI residual." },
    "bandwidth_tested": { "type": "array", "items": {"type": "number"}, "description": "(Output by writer) Grid of W values tested if auto-selected."}
  },
  "required": ["region_params"],
  "allOf": [
    {
      "if": { "properties": { "method": { "const": "global_plus_roi_residual" } } },
      "then": { "required": ["roi_residual_params"] }
    }
  ]
}
```

*   **Key Parameters (Refined):**
    *   `method`: Explicitly includes `"global_plus_roi_residual"`.
    *   `bandwidth`: Type `number` with `exclusiveMinimum: 0`.
    *   `solver_params`: `additionalProperties: true` retained but `maxProperties: 6` added for sanity.
    *   `eigenvalues`: `maxItems: 4096` added.
    *   Added `k_roi_actual`, `bandwidth_roi_actual` for clarity when `method="global_plus_roi_residual"`.
*   The `lna:::default_params("spat.slepian")` helper would use schema defaults.

#### 2. Slepian Construction (`forward_step.spat.slepian`)

The core of the Slepian construction involves solving a generalized eigenproblem.

*   **Mathematical Formulation:**
    The discrete Slepian functions are eigenvectors of the symmetric concentration operator $\mathbf{C} = \mathbf{P A P}$, where `A` is typically a spatial adjacency matrix (e.g., 6- or 26-connectivity, default 6-conn.) and `P` is a diagonal masking matrix. After restricting to the in-mask subspace Ω, this becomes the sparse matrix $\tilde A$ solved below. For an arbitrary spatial support (mask) Ω within a larger N-voxel domain:
    1.  Define a spatial adjacency matrix `A` (e.g., 6-connectivity for voxels is a good default).
    2.  Define a diagonal masking matrix `P` where `P_ii = 1` if voxel `i` is in Ω, and `0` otherwise.
    3.  The problem is to solve the generalized eigenproblem `Av = λPv`. *(Note: Explicit formation of `L⁻¹` in `S = PL⁻¹P` is avoided; the problem is solved directly.)*

*   **Numerical Solution within LNA `forward_step`:**
    1.  **Mask Subspace:** Identify the set of in-mask voxel indices Ω from the input mask. Let $N_\Omega$ be the number of voxels in Ω.
    2.  **Compressed Adjacency/Laplacian:** Extract the sub-matrix $\tilde A = A[\Omega,\Omega]$ (size $N_\Omega \times N_\Omega$).
    3.  **Standard Eigenproblem:** The generalized problem `Av = λPv` restricted to Ω reduces to an ordinary symmetric eigenproblem on this subspace: $\tilde A \tilde v = \lambda \tilde v$.
    4.  **Eigen-Solver:** Use `RSpectra::eigs_sym()` with a matrix-free operator for $\tilde A \tilde v$.
        ```R
        op_fun <- function(x, args_list) { drop(args_list$tilde_A_masked %*% x) }
        eig_result <- RSpectra::eigs_sym(op_fun, k = k_desired, n = N_omega,
                                         A = list(tilde_A_masked = tilde_A),
                                         opts = params$solver_params %||% list(tol = 1e-6)) # Use solver_params
        ```
        Track convergence and consider adaptive `ncv` or `max_iter` based on `solver_params`.
    5.  **Padding:** Eigenvectors $\tilde v$ are padded with zeros for out-of-mask voxels to form $v$.
    6.  Keep top `k_actual` eigenvectors; store `eig_result$values[1:k_actual]` in `desc$params$eigenvalues`.

*   **Adaptive Bandwidth (`bandwidth_selection="auto_energy"`):**
    *   Iterate over a candidate grid of `W` values (e.g., `candidate_W <- 2^(seq(log2(W_min), log2(W_max), length.out=6))`). For each `W`, determine `keep` components for `params$target_energy`.
    *   Choose best `W`. Record `params$bandwidth_actual`, `params$k_actual`, `params$bandwidth_tested`.

*   **Hybrid Global + ROI Residual (`method="global_plus_roi_residual"`):**
    1.  Fit global Slepian basis (`global_basis`, `k_global` components) on the whole mask. Store at `/basis/slepian_global/matrix`.
    2.  Project input data: `global_coeffs = data %*% global_basis`.
    3.  Reconstruct global part: `data_global_recon = global_coeffs %*% t(global_basis)`. *Projecting onto the orthonormal global basis before computing the residual effectively removes variance captured by global components, helping to ensure that the subsequent ROI-specific basis captures orthogonal information relative to the global basis within that ROI.*
    4.  Compute ROI residual: `data_roi_residual = data_roi_original - data_global_recon_roi`.
    5.  Fit ROI Slepian basis (`roi_residual_basis`, `k_roi` components) on `data_roi_residual` using ROI mask. Store at `/basis/blocks/slepian_roi_{roi_name}/matrix`.
    6.  Project `data_roi_residual`: `roi_residual_coeffs = data_roi_residual %*% roi_residual_basis`.
    7.  **Bookkeeping:**
        *   Descriptor lists both basis datasets.
        *   Coefficients `global_coeffs` and `roi_residual_coeffs` are concatenated column-wise into `/scans/{run_id}/embedding/coefficients`.
        *   `/spatial/block_table` gets a new row for the ROI (e.g., `block_id = B+1`). `coeff_offset_block` (in the `spat.slepian` or `embed` descriptor) for this ROI block is `last_total_offset + k_global` (assuming global coeffs are first). The value of `coeff_offset_block` is monotonically increasing.
        *   `/spatial/partition` dataset marks voxels in the ROI as belonging to `block_id = B+1`. *Writer must ensure that after adding the ROI row to `block_table` and defining its `block_id`, the `/spatial/partition` entries for all voxels within that ROI are updated to this new `block_id`, overriding any previous octree block assignment ("ROI wins" rule).*
    *   **Reader Logic:** Sum contributions: `Recon_voxel = GlobalBasis_voxel_col ⋅ GlobalCoeffs + (if partition[v]==ROI_ID then ROI_Basis_voxel_col ⋅ ROI_Coeffs else 0)`.

*   **Data Storage:** As previously described, with eigenvalues stored in the descriptor.

#### 3. Inverse Step (`invert_step.spat.slepian`)

Analogous to `invert_step.basis`, loading appropriate basis/bases and coefficient slices based on `desc` (including `method` for hybrid reconstruction) and `handle$subset`.

#### 4. Capabilities Flags in Descriptor

```json
"capabilities": {
  "supports_spatial_subsetting": "block-aware", // Or "full" if global only
  "supports_temporal_subsetting": true,
  "supports_roi_query": true // If method supports specific ROI optimization
}
```
*   `supports_roi_query: true` indicates that if a read request's mask exactly matches a stored ROI mask (e.g., from `roi_residual_params$mask_path`), the reader can optimize by potentially skipping the global basis reconstruction if only ROI-specific information is desired and the method allows separable reconstruction.

### C. Innovative Extra Capability: Leakage-Bounded Functional Connectivity

The Slepian eigenvalues ($\lambda_i$) quantify component energy concentration within support Ω. This enables bounding inter-ROI leakage in FC estimates.

*   **Concept: "Leakage-Aware Correlation"**
    *   For two ROIs, Ω₁ and Ω₂, compressed with Slepian bases, with coefficient time series $c^{(1)}_i(t)$, $c^{(2)}_j(t)$ and eigenvalues $\lambda^{(1)}_i$, $\lambda^{(2)}_j$.
    *   A leakage-corrected Pearson correlation $\rho^\star$:
        $begin:math:display$
        \rho^\star = \frac{\sum_i \lambda_i^{(1)}\lambda_i^{(2)}\,\mathrm{Cov}(c^{(1)}_i,c^{(2)}_i)}
                           {\sqrt{\sum_i (\lambda_i^{(1)})^2 \mathrm{Var}(c^{(1)}_i)}\,
                            \sqrt{\sum_i (\lambda_i^{(2)})^2 \mathrm{Var}(c^{(2)}_i)}}.
        $end:math:display$
        This discounts contributions from poorly localized components.

*   **LNA Implementation:**
    *   **Writer (`forward_step.spat.slepian`):** Stores `eigenvalues` in `desc$params`.
    *   **New Analysis Transform (Optional, Reader-Side or Persisted):**
        ```json
        {
          "type": "analysis.leakage_fc",
          "version": "0.1",
          "params": {
            "roi_pair_descriptors": [
              { "coeffs_path": "/scans/run-01/embedding/coefficients", "slepian_desc_path": "/transforms/00_spat.slepian.json", "roi_name_in_slepian_desc": "primary" }, // or "roi_dlpfc" if hybrid
              { "coeffs_path": "/scans/run-01/embedding/coefficients", "slepian_desc_path": "/transforms/01_spat.slepian.json", "roi_name_in_slepian_desc": "primary" }  // Example if second ROI has its own Slepian transform
            ],
            "component_indices_roi1": [1, 2, 5], // Optional: subset of components to use
            "component_indices_roi2": [1, 3, 4]
          },
          "datasets": [ // Optional: if persisting results
             { "path": "/analysis/leakage_fc/results/dlpfc_precuneus_corr", "role": "scalar_result" }
          ],
          "inputs": ["coefficients_roi1_key", "coefficients_roi2_key"], // Stash keys for coefficients
          "outputs": ["rho_star_value"] // Stash key for output
        }
        ```
        *   **Inputs:** Takes coefficient time series (from stash or HDF5) and paths to the relevant `spat.slepian` descriptors to fetch `eigenvalues`.
        *   **Logic:** Computes $\rho^\star$.
        *   **Outputs:** Corrected correlation(s) to stash. Optionally writes results to a new dataset if persistence is desired.

*   **Benefit:** Analytically justifiable leakage correction without complex deconvolution, unique to Slepian bases.

### Summary for `spat.slepian` Implementation

*   **Robust Construction:** Solves `Av = λPv` via a subspace eigenproblem $\tilde A \tilde v = \lambda \tilde v$ on the mask support Ω.
*   **Hybrid Global + ROI Residuals:** Supported via `method="global_plus_roi_residual"`, with clear HDF5 layout and `block_table`/`partition` updates. Writers must ensure ROI voxels in `/spatial/partition` are correctly assigned to the ROI block ID, taking precedence.
*   **Adaptive Bandwidth:** Supported via `bandwidth_selection="auto_energy"`.
*   **Schema:** Updated for new methods, parameter constraints (`exclusiveMinimum`, `maxItems`), and output fields.
*   **Capabilities:** Includes `supports_roi_query` flag.
*   **LNA Blueprint Impact:** Minimal. Requires new schema, `lna:::default_params` update, and new unit tests.
*   **Advanced Analysis:** Slepian eigenvalues enable novel reader-side analyses like leakage-bounded FC, which can be formalized as an `analysis.leakage_fc` transform.

The `spat.slepian` transform, with these refinements, is a powerful and versatile addition, enhancing LNA's compression and analytical capabilities significantly while maintaining compatibility with the v1.4 framework.