## Appendix: `spat.hkwp` - Heat-Kernel Wavelet Packet Spatial Transform

This appendix details the `spat.hkwp` transform, an advanced spatial basis for LNA rooted in spectral graph theory and information-theoretic wavelet packet analysis. It constructs an adaptive, multi-resolution basis derived from the diffusion process (heat kernel) on the underlying voxel graph, offering excellent spatio-spectral localization and energy compaction. This transform requires only sparse linear algebra (primarily Lanczos iteration via `RSpectra`) and fits seamlessly into the LNA v1.4 framework.

*(Key References: Coifman, R. R., & Maggioni, M. (2006). Diffusion wavelets. Applied and Computational Harmonic Analysis, 21(1), 53-94; Hammond, D. K., Vandergheynst, P., & Gribonval, R. (2011). Wavelets on graphs via spectral graph theory. Applied and Computational Harmonic Analysis, 30(2), 129-150; Harte, D. (2019). Numerical evaluation of the matrix exponential. Journal of Computational and Applied Mathematics, 354, 470-491 – for Chebyshev approximation of matrix exponential; Coifman, R. R., & Wickerhauser, M. V. (1992). Entropy-based algorithms for best basis selection. IEEE Transactions on Information Theory, 38(2), 713-718; Donoho, D. L. (1995). De-noising by soft-thresholding. IEEE Transactions on Information Theory, 41(3), 613-627 – for principles behind information-based splitting.)*

### 1. Core Idea: Heat-Kernel Wavelet Packets

*   **Heat Kernel Smoothing:** Convolving a brain volume with the heat kernel $e^{-t\Delta}$ (where $\Delta$ is the graph Laplacian of the voxel connectivity and $t$ is a diffusion time-scale) produces a smoothed version of the data. The extent of smoothing is controlled by $t$.
*   **Dyadic Diffusion Wavelets:** By iterating this idea at dyadically spaced time-scales ($T_0 = I, T_1 = e^{-2^{-1}\Delta}, T_2 = e^{-2^{-2}\Delta}, \ldots$), a multi-resolution family of operators is formed. The differences $W_j = T_{j-1} - T_j$ define graph diffusion wavelets. These wavelets:
    *   Are orthogonal across scales.
    *   Localize simultaneously in space (decaying like a Gaussian with voxel distance) and frequency (each scale $j$ acts as a dyadic band-pass filter).
    *   Can be made to form an orthonormal basis or a tight frame, allowing signal reconstruction.
*   **Adaptive Wavelet Packets:** Instead of a fixed wavelet tree, the HKWP transform adaptively decides whether to further decompose (split) or keep a sub-band at each scale based on an information criterion (e.g., entropy or BIC of the projected data in that sub-band relative to a `lambda_threshold`). This yields a basis optimized for the specific fMRI data, allocating more basis atoms to complex regions and fewer to smooth regions.
*   **Orthonormal Basis (Post-Pruning):** After the adaptive packet pruning, the collection of chosen leaf packet atoms (basis vectors) is explicitly re-orthogonalized (e.g., via thin QR decomposition or Gram-Schmidt) to ensure they form a truly orthonormal basis. This guarantees that their "eigenvalues" (operator norms) are all 1.0, simplifying subsequent quantization.

### 2. Advantages for fMRI Data

| Property                     | Practical Win for fMRI                                                                                             | Comparison with Other Bases          |
| :--------------------------- | :----------------------------------------------------------------------------------------------------------------- | :----------------------------------- |
| Physics-Informed             | Intrinsically models haemodynamic blurring via heat diffusion on the grey-matter graph.                            | PCA/DCT ignore anatomy/spatial coupling. |
| Dyadic Multiresolution       | Natural for whole-brain context and ROI zoom-ins (coarse view from early scales, refined with later scales).         | Octree is spatial-only multiscale.   |
| Wavelet-Packet Adaptivity    | Data-driven allocation of basis functions: more where data is complex (e.g., vascular borders), less where smooth. | Slepian/Fixed Wavelets use fixed `k`.  |
| Sparse Transform Operators   | Each $W_j$ is sparse, making forward/inverse operations efficient ($O(k N_s)$, $N_s \ll N_{voxels}$).               | Graph-Fourier bases are dense.       |
| Orthonormal Basis Atoms      | All final basis atoms have norm 1, simplifying quantization (uniform dynamic range for coefficients).              | Block-PCA needs per-component scales. |
| Analysis-Rich Representation | Atoms tagged by scale and (optionally) centroid/bounding box enable multiscale connectivity, lag analysis, etc.      | Less direct for PCA/DCT.             |

### 3. LNA Implementation Strategy for `spat.hkwp`

#### 3.1. Schema (`inst/schemas/spat.hkwp.schema.json`)

```json
{
  "type": "object",
  "title": "Parameters for Heat-Kernel Wavelet Packet (spat.hkwp) transform",
  "$id": "https://neurocompress.org/schemas/lna/2.0/spat.hkwp.schema.json",
  "$schema": "http://json-schema.org/draft-07/schema#",
  "properties": {
    "type": { "const": "spat.hkwp" },
    "version": { "const": "1.0" },
    "J_max": {
      "type": "integer", "minimum": 1,
      "description": "Deepest dyadic scale explored (e.g., 2^-J_max relates to smallest diffusion time)."
    },
    "split_rule": {
      "enum": ["entropy", "bic"], "default": "entropy",
      "description": "Criterion for splitting a wavelet packet node."
    },
    "lambda_threshold": {
      "oneOf": [
        { "type": "number", "exclusiveMinimum": 0 },
        { "type": "array", "items": { "type": "number", "exclusiveMinimum": 0 }, "description": "Array of length J_max + 1, for scales 0 to J_max." },
        { "type": "string", "pattern": "^/.*", "description": "HDF5 path to a 3D lambda map (voxel-wise thresholds)." }
      ],
      "default": 0.15,
      "description": "Threshold for the split rule. If a single number (including the default 0.15), it is recycled to length J_max+1 and applied to all scales. If an array, it must have J_max+1 elements for per-scale thresholds. If an HDF5 path, provides a spatially varying lambda map."
    },
    "cheb_order": {
      "type": "integer", "minimum": 3, "maximum": 30, "default": 6,
      "description": "Order of Chebyshev polynomial approximation for heat kernel operator. Higher orders increase accuracy but also computation time."
    },
    "lap_spec_hash": {
      "type": "string", "pattern": "^[a-f0-9]{40}$",
      "description": "(Output by writer, optional) SHA-1 hash of the first ~32 eigenvalues of the graph Laplacian. Used for quick equality check of the underlying geometry if attempting to reuse a precomputed basis across different LNA write operations or runs."
    },
    "k_actual": {
      "type": "integer",
      "description": "(Output by writer) Actual number of wavelet packet atoms kept in the final orthonormal basis."
    },
    "storage_order": {
      "enum": ["component_x_voxel", "voxel_x_component"],
      "default": "component_x_voxel"
    },
    "eigenvalues": {
      "type": "array", "items": {"type":"number", "const": 1.0},
      "maxItems": 8192,
      "description": "(Output by writer, optional) Operator norm for each atom; all will be 1.0 if the final basis is re-orthogonalized as recommended."
    },
    "hkwp_tree_structure_path": {
      "type": "string", "pattern": "^/.*",
      "description": "(Output by writer, optional) HDF5 path to a ~1-2kB JSON descriptor containing the HKWP tree structure for fingerprinting (e.g., /transforms/00_spat.hkwp.tree.json)."
    }
  },
  "required": ["J_max"]
}
```

#### 3.2. Adjacency, Laplacian, and Heat Kernel Operators

1.  **Graph Construction & Masking (Writer Side):**
    *   From the input mask (e.g., `handle$mask_info$mask`), extract the set of in-mask voxel indices Ω.
    *   Build the sparse adjacency matrix `A` for the full N-voxel domain (e.g., 6-connectivity).
    *   Extract the sub-matrix $\tilde A = A[\Omega,\Omega]$ (an $N_\Omega \times N_\Omega$ symmetric sparse matrix).
    *   Compute the graph Laplacian for this sub-graph: $\tilde \Delta = \tilde D - \tilde A$, where $\tilde D$ is the degree matrix for $\tilde A$.
    *   If the mask Ω is disconnected, this process is applied to each connected component $\Omega_c$ independently, yielding $\tilde \Delta_c$. The final basis atoms are concatenated.
2.  **Heat Kernel Application $e^{-t\tilde{\Delta}}x$ via Chebyshev Approximation:**
    *   (Algorithm Crib Sheet B-1)
    *   Symbols:
        *   $\tilde{\Delta}$: Graph Laplacian of the *in-mask subgraph* $\Omega$.
        *   $\lambda_{min}, \lambda_{max}$: Extreme eigenvalues of $\tilde{\Delta}$ (estimate using power iterations on $\tilde{\Delta}$). $\lambda_{min}$ is often 0 for connected graphs.
        *   $s = (\lambda_{max} - \lambda_{min})/2$; $c = (\lambda_{max} + \lambda_{min})/2$.
        *   $\tilde{\Delta}' = (\tilde{\Delta} - cI)/s$ (scaled Laplacian with eigenvalues in $[-1,1]$).
    *   Approximation: $e^{-t\tilde{\Delta}} x \approx \sum_{m=0}^{M} \beta_m(t) T_m(\tilde{\Delta}') x$.
    *   Recursion for $T_m(\tilde{\Delta}')x$:
        *   $T_0x
*   Recursion for $T_m(\tilde{\Delta}')x$:
    *   $T_0x = x$
    *   $T_1x = \tilde{\Delta}' x$
    *   $T_{m+1}x = 2 \tilde{\Delta}' (T_m x) - T_{m-1}x$ for $m \ge 1$.
*   Coefficients $\beta_m(t)$: (Using the formulation involving modified Bessel functions $I_m$, suitable for $e^{-\tau X}$ where $X$ has spectrum in $[-1,1]$ and $\tau = ts$).
    *   Let $\alpha = ts$.
    *   $\beta_0(t) = I_0(\alpha) e^{-\alpha} e^{-tc}$
    *   $\beta_m(t) = 2 I_m(\alpha) e^{-\alpha} e^{-tc}$ for $m \ge 1$.
    *   *(Note: The $e^{-tc}$ term accounts for the shift back from $\tilde{\Delta}'$ to $\tilde{\Delta}$ within the sum. Alternatively, apply $e^{-tc}$ to the final sum).*
*   $M$: Order of approximation, e.g., `params$cheb_order` (default 6). Typically $M \approx \lceil \max(10, 5 t \lambda_{max}) \rceil$ for desired accuracy, where $\lambda_{max}$ is for the original $\tilde{\Delta}$.

#### 3.3. `forward_step.spat.hkwp` (Writer-Side Algorithm Sketch)

Input: Graph Laplacian $\tilde{\Delta}$ (for in-mask voxels Ω), signal matrix $X$ (Time $\times N_\Omega$), parameters from `desc$params`.

1.  **Initialize:** `J = params$J_max`. Load/prepare `lambda_threshold` vector (`lambda_vec`) of length `J+1` based on scalar, array, or map input (scalar default `0.15` is repeated for all scales). `T_prev_operator_on_X = X`. `kept_wavelet_atoms_raw = []` (list to store atom matrices). `atom_info_gains = []`. `atom_scales = []`.
    *   (Optional) Compute `lap_spec_hash` from first ~32 eigenvalues of $\tilde{\Delta}$ (via `RSpectra::eigs_sym(tilde_Delta, k=32, which="SM")`) and store in `desc$params`.
2.  **Adaptive Wavelet Packet Tree Construction (Iterative over scales $j=0 \ldots J$ and nodes at each scale):**
    *   Maintain a queue/list of active nodes (initially one node: the full data $X$ at scale $j=0$, associated with $T_0=I$). Each node contains {projected data $X_{node}$, parent operator $T_{parent}$, scale $j_{node}$, spatial support/centroid (optional)}.
    *   For each node in the queue at scale $j_{node}$:
        a.  If $j_{node} > J$, add its associated operator/basis to `kept_wavelet_atoms_raw` and continue.
        b.  `t_current = 2^{-j_{node}}` (or other dyadic scaling for diffusion time).
        c.  Compute current diffusion operator's action $T_{current\_node\_action\_on\_X} = e^{-t_{current}\tilde{\Delta}} X_{node}$ (via Chebyshev approximation on $X_{node}$).
        d.  Form band-pass (wavelet) operator's action $W_{node\_action\_on\_X} = T_{parent\_node\_action\_on\_X} - T_{current\_node\_action\_on\_X}$.
        e.  Compute information metric $I(W_{node\_action\_on\_X})$ (e.g., entropy of normalized column energies, see Crib Sheet B-2).
        f.  `current_lambda = lambda_vec[j_{node}+1]` (or lookup from map based on node centroid).
        g.  If $I(W_{node\_action\_on\_X}) > current\_lambda$ AND $j_{node} < J$:
            *   Node is split. Add two children to queue for processing at scale $j_{node}+1$:
                *   Child 1 (low-pass): Data $T_{current\_node\_action\_on\_X}$, Parent Op $T_{current}$.
                *   Child 2 (band-pass/detail): Data $W_{node\_action\_on\_X}$, Parent Op $W_{current}$. *(This part needs care: the children for splitting a packet are usually derived from further filtering operations, not just passing $T$ and $W$ actions as "data" to the next level. A more standard packet decomposition would use $W_j$ to project $X$ onto two orthogonal subspaces using $W_j$ and $T_j$ as projectors, then recurse on those subspaces.)*
            *   **Simplified/Corrected Packet Logic:** A common wavelet packet approach involves, at each node, projecting the current data onto two orthogonal subspaces (e.g., using $W_j$ and $T_j$ or similar decomposition filters derived for that node) and then deciding whether to keep the parent node's projection or the children's projections.
            *   Let's assume for simplicity that if a node is *not* split, the operator $W_j$ (or its effective action that generated the current node's data) becomes a "leaf" operator.
        h.  Else (node is a leaf):
            *   The operator that produced $X_{node}$ (which is a specific combination of $W_k$'s from previous scales) is a leaf atom/group of atoms. Add this operator (or its matrix form if small enough) to `kept_wavelet_atoms_raw`.
            *   Record `info_gain` (e.g., `parent_entropy - current_node_entropy`) and `scale_index = j_{node}`.
3.  **Construct & Orthonormalize Basis:**
    *   Form preliminary basis matrix $B_{raw}$ ( $N_\Omega \times K_{raw}$) by taking all columns from all matrices in `kept_wavelet_atoms_raw`.
    *   Run a thin QR decomposition: `qr_decomp = qr(B_raw)` (if $B_{raw}$ stored as atoms $\times$ voxels, then `qr(t(B_raw))`).
    *   Extract the orthonormal basis $B$ ($N_\Omega \times k_{actual}$), ensuring columns have norm 1. $k_{actual} = \text{ncol}(B)$.
4.  **Pad & Store Data:**
    *   Pad atoms in $B$ with zeros to full N-voxel space if $\Omega$ was a sub-mask.
    *   Basis matrix $B$ (now $N \times k_{actual}$ or $k_{actual} \times N$) to `/basis/hkwp/global/matrix` (or per-block path).
    *   Coefficients $C = X_{original} \cdot B_{padded}^\top$ (if $B$ is $k \times N$) to `/scans/{run_id}/embedding/coefficients`.
    *   Record `k_actual` in `desc$params`. If basis is orthonormal, `desc$params$eigenvalues = rep(1.0, k_actual)`.
    *   Store `atom_scales` (vector of length `k_actual`) to `/basis/hkwp/atom_scales`.
    *   Store `atom_info_gain` (vector of length `k_actual`) to `/basis/hkwp/atom_info_gain`.
    *   (Optional) Compute and store atom bounding boxes ($k_{actual} \times 6$) to `/basis/hkwp/atom_bounding_boxes`.
    *   (Optional) Serialize HKWP tree structure to JSON, write to `/transforms/NN_spat.hkwp.tree.json` (e.g., `00_spat.hkwp.tree.json`), and store this path in `desc$params$hkwp_tree_structure_path`.
5.  **Update Plan & Stash:** Add dataset definitions to `handle$plan`. Output coefficient path/reference to `handle$stash`.

#### 3.4. `invert_step.spat.hkwp` (Reader-Side)

1.  Load basis matrix `B` (from path in `desc$datasets`).
2.  Load coefficient matrix `C` (from path in `desc$datasets`).
3.  **Progressive Reconstruction Options (from `handle$subset` or direct args):**
    *   `max_scale`: If `handle$subset$max_scale` is provided, load `/basis/hkwp/atom_scales`. Filter columns of `B` and `C` to include only atoms where `atom_scales <= max_scale`.
    *   `info_atoms` or `info_thresh`: If `handle$subset$info_atoms` or `info_thresh` provided, load `/basis/hkwp/atom_info_gain`. Order atoms by `info_gain`. Filter columns of `B` and `C`.
        ```R
        # Pseudocode for filtering inside invert_step
        if (!is.null(handle$subset$info_atoms) || !is.null(handle$subset$info_thresh)) {
            gain_path <- # find path for role 'atom_information_metric' in desc$datasets
            atom_info_gain_values <- handle$h5[[gain_path]]$read()
            if (!is.null(handle$subset$info_atoms)) {
                keep_indices <- order(atom_info_gain_values, decreasing = TRUE)[seq_len(min(handle$subset$info_atoms, length(atom_info_gain_values)))]
            } else { # info_thresh
                keep_indices <- which(atom_info_gain_values >= handle$subset$info_thresh)
            }
            # Assuming B is atoms x voxels, C is time x atoms
            B_matrix <- B_matrix[keep_indices, , drop=FALSE] 
            C_matrix <- C_matrix[, keep_indices, drop=FALSE]
        }
        ```
4.  Handle standard `handle$subset` (ROI slicing on rows of `B` if $B$ is voxels $\times$ atoms, or columns if atoms $\times$ voxels; time slicing on rows of `C`). The atom bounding boxes in `/basis/hkwp/atom_bounding_boxes` can be used to quickly pre-filter atoms for ROI queries before loading full basis atom data.
5.  Reconstruct: $X_{hat} = C \cdot B$ (if $B$ is atoms $\times$ voxels). Or $X_{hat} = B \cdot t(C)$ (if $B$ is voxels $\times$ atoms). Ensure dimensions align.
6.  Place $X_{hat}$ into `handle$stash`.

#### 3.5. Capabilities Flags in Descriptor

```json
{
  "capabilities": {
    "supports_spatial_subsetting": "block-aware", 
    "supports_temporal_subsetting": true,
    "supports_scale_subsetting": true,
    "supports_information_subsetting": true 
  }
}
```

#### 3.6. Coefficient Offset Bookkeeping (If HKWP is used block-wise)

If `spat.hkwp` is adapted for an octree approach or a "global + HKWP_ROI_residual" method (where HKWP provides one or more blocks alongside other basis types):

*   Global HKWP coefficients (if any) start at column offset 0.
*   Each HKWP-derived block (e.g., an octree leaf processed by HKWP, or an HKWP-processed ROI) will have an entry in `/spatial/block_table` (e.g., `block_id = i`).
*   The `coeff_offset_block` array (stored in the descriptor of the transform *producing the final embedding*, which could be `spat.hkwp` itself or a subsequent `embed` transform) maps these `block_id`s to their starting column in the concatenated `/scans/<run>/embedding` HDF5 dataset.
*   **Spec Note:** `coeff_offset_block` is a vector whose length equals `nrow(/spatial/block_table)`. `coeff_offset_block[i]` stores the 0-based starting column index for the coefficients corresponding to the block defined in the $i$-th row of `/spatial/block_table` (i.e., `block_id = i`). The implicit global slice always has an offset of 0.
*   Reader logic for coefficient selection:
    ```R
    # part_id is the 1-based block_id from /spatial/partition for the current voxel
    # k_global is known from the global part of the descriptor
    # k_block_vec holds k_actual for each block_id in /spatial/block_table
    # coeff_offset_block is 1-indexed by block_id
    
    if (part_id == 0L) { # Voxel belongs to the implicit global basis
        col_start_0based <- 0L
        k_for_this_part  <- k_global
    } else { # Voxel belongs to a specific block_id > 0
        col_start_0based <- coeff_offset_block[part_id] 
        k_for_this_part  <- k_block_vec[part_id] 
    }
    # R uses 1-based slicing
    cols_to_select_R <- (col_start_0based + 1L):(col_start_0based + k_for_this_part)
    
    coefficients_for_this_part <- embedding_matrix_full[time_index_t, cols_to_select_R, drop = FALSE]
    ```

### 4. Advanced Capabilities & Usage

#### 4.A. Progressive Reconstruction (by Scale and Information Content)

Enabled by storing `atom_scales` and `atom_info_gain`, and respective reader API arguments (`max_scale`, `info_atoms`, `info_thresh`).

#### 4.B. HKWP Fingerprint Export & Analysis

The serialized tree structure (path stored in `desc$params$hkwp_tree_structure_path`, JSON data at e.g., `/transforms/00_spat.hkwp.tree.json`) allows for rich inter-dataset comparisons.
```R
extract_hkwp_fingerprint <- function(lna_file, hkwp_transform_descriptor_name = "00_spat.hkwp.json") {
  h5 <- hdf5r::H5File$new(lna_file, "r")
  on.exit(h5$close_all())
  main_hkwp_desc <- lna:::read_json_descriptor(h5[["/transforms"]], hkwp_transform_descriptor_name)
  tree_path <- main_hkwp_desc$params$hkwp_tree_structure_path
  if (is.null(tree_path) || !h5$exists(tree_path)) {
    stop("HKWP tree structure path not found or invalid in descriptor.")
  }
  tree_json_str <- lna:::h5_read(h5, tree_path) # Corrected to use lna internal h5_read
  jsonlite::fromJSON(tree_json_str, simplifyVector = FALSE)
}
```

#### 4.C. Basis Re-use via `lap_spec_hash`

The `desc$params$lap_spec_hash` (SHA-1 of first ~32 Laplacian eigenvalues) allows writers to detect if an identical geometry (and thus potentially HKWP basis, if other params match) has been processed, enabling re-use.

#### 4.D. Fast ROI Queries via Atom Bounding Boxes

Storing `/basis/hkwp/atom_bounding_boxes` ($k_{actual} \times 6$ array of `x0,x1,y0,y1,z0,z1` for each atom) allows readers to quickly identify atoms intersecting a query ROI, loading only relevant basis data for partial reconstructions.

### 5. Implementation Feasibility & Testing

*   **R Implementation:** The core algorithm relies on `igraph` (graph construction), `Matrix` (sparse matrices), `RSpectra` (eigen-solvers for Chebyshev operator), `hdf5r`, and `jsonlite`. A prototype is achievable in pure R.
*   **Sanity-Check Test Battery:**
    1.  **Laplacian Trace Check:** Compare `trace(e^{-t\tilde{\Delta}})` from Chebyshev sum vs. sum of $e^{-t\lambda_i}$ from `RSpectra::eigs_sym(tilde_Delta, ...)` eigenvalues. Relative error < 1%.
    2.  **Wavelet Packet Partition of Identity (Orthonormality):** After re-orthogonalization, verify $B B^\top \approx I_{N_\Omega}$ (if $B$ is $N_\Omega \times k_{actual}$) or $B^\top B \approx I_k$ (if $B$ is $k_{actual} \times N_\Omega$). For a tight frame before re-orthogonalization, verify $\sum_{leaf} W_{leaf}^\top W_{leaf} \approx I$ on Ω by applying to random vectors; max deviation < $10^{-9}$. After re-orthogonalization of $B$, verify $B B^T \approx P_\Omega$ (projector onto span of atoms in $\Omega$) and $B^T B = I_k$.
    3.  **Perfect Reconstruction (float64, no subsetting):** $X_{hat} = C B$ (or $B t(C)$) should match original $X$ on Ω with RMSE/σ_X < $10^{-12}$.
    4.  **Energy Compaction (Info Gain):** Cumulative energy of `atom_info_gain` (sorted) should reach >95% well before $k_{actual}/2$ for typical fMRI.
    5.  **Scale Consistency (Progressive Reconstruction):** Reconstructing with atoms up to scale $j$, then adding atoms of scale $j+1$, should yield same result as direct reconstruction up to scale $j+1$.
    6.  **ROI Locality (Bounding Box):** Reconstructing ROI using only atoms whose bounding boxes intersect it should have <1% energy deviation from full ROI reconstruction (for well-localized atoms).
    7.  **Chebyshev Order Sensitivity:** Decreasing `cheb_order` should monotonically increase heat kernel approximation error.
    8.  **Fingerprint Invariance:** HKWP on identical data (fixed seed for any stochastic parts, though ideally none) should yield identical tree structure & `lap_spec_hash`.
    9.  **Progressive Preview (Info Atoms):** RMSE should monotonically decrease as `max_info_atoms` increases.
    10. **Lambda-Map Functionality:** Test with a lambda map having distinct regions; verify denser splitting (more leaf atoms per unit volume) in low-lambda regions.

### 6. (Future Vista) Diffusion-Time Adaptive RCNN (`spat.hk_resnet`)

A potential future transform using stacked, fixed-weight graph convolutions ($e^{-t_0\Delta}$) followed by orthogonalization (QR) to create a compact, hierarchical basis, mimicking a Graph ResNet without explicit training.

### Conclusion

The `spat.hkwp` transform provides a powerful, physics-informed, and data-adaptive spatial compression and analysis tool for LNA. Its implementation, while sophisticated, leverages standard numerical techniques available in R and integrates cleanly with the LNA v1.4 specification. The additional features for progressive reconstruction, fingerprinting, and basis re-use significantly enhance its utility.