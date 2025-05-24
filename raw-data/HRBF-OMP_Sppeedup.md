This is an excellent, targeted review focusing on the Rcpp optimization aspects for the HRBF-related transforms. The plan is sound, and the suggestions are valuable for robustness and diagnostics.

Let's synthesize this into a **Definitive Proposal for Rcpp Optimizations within the HRBF Transform Ecosystem for LNA.** This will serve as a detailed guide for an engineer.

---

## Definitive Proposal: Rcpp/Eigen Optimizations for HRBF Transforms in LNA

This document outlines the targeted Rcpp/Eigen optimizations for the Hierarchical Radial Basis Function (HRBF) transform family within the LNA package. The goal is to significantly improve performance for computationally intensive kernels while maintaining the high-level R orchestration for clarity and ease of integration. These optimizations are designed to be CPU-only and integrate seamlessly with the existing LNA architecture.

*(This proposal assumes the overall LNA framework and the R-level logic for HRBF transforms – `spat.hrbf`, `spat.hrbf_project`, `basis.empirical_hrbf_compressed`, `embed.transfer_hrbf_basis` – are defined as per previous discussions. This focuses on porting specific hotspots to C++.)*

### 1. Identified Performance Kernels for Rcpp Porting

Based on profiling and algorithmic structure, the following kernels are prime candidates for Rcpp/Eigen implementation to achieve substantial speedups:

1.  **P1: HRBF Basis Atom Generation (`hrbf_atoms_cpp`):**
    *   **Current R Bottleneck:** The R-level loop in `hrbf_basis_from_params` (or similar helpers) that iterates through $K_{atoms}$ and for each, iterates through $N_{mask\_voxels}$ to calculate distances and kernel values, followed by R-level sparse matrix triplet collection (`i_idx`, `j_idx`, `x_val` via repeated `c()`).
    *   **Rcpp Solution:** A single C++ function to compute all non-zero atom values for all atoms and construct an `Eigen::SparseMatrix` directly from triplets.

2.  **P2: Orthogonal Matching Pursuit (OMP) Solver (`omp_encode_cpp`):**
    *   **Current R Bottleneck (if pure R OMP used):** The greedy iteration involving repeated projections ($D^T r$), finding max correlation, updating active set, and solving least squares on the growing sub-dictionary ($D_{sub}^T D_{sub} c = D_{sub}^T y$). Matrix inversions or `solve()` calls in R within a loop are slow.
    *   **Rcpp Solution:** Implement OMP in C++ using Eigen for sparse matrix-vector products and an efficient solver for the least-squares subproblem (e.g., incremental Cholesky `SimplicialLLT`).

3.  **P3: Poisson-Disk Sampling (`poisson_disk_sample_component_rcpp`):**
    *   **Current R Bottleneck:** The sequential rejection loop is $O(N_{selected} \cdot N_{candidates})$ in its distance checking part per accepted sample, and `rbind` for `selected` points is inefficient.
    *   **Rcpp Solution:** C++ implementation with efficient distance checks (early exit) and potentially a voxel-grid hash for amortized $O(1)$ neighbor lookups for denser regions. `std::vector` for storing selected points.

*(Note: Connected component labeling, `label_components_6N_rcpp`, was previously identified and its Rcpp sketch is assumed to be a prerequisite or part of this optimization effort if not already done.)*

### 2. Rcpp/Eigen Implementation Sketches & Details

#### 2.1. P1: `hrbf_atoms_cpp` (HRBF Basis Generation)

*   **File:** `src/lna_hrbf_rcpp.cpp`
*   **Dependencies:** `// [[Rcpp::depends(RcppEigen)]]`
*   **Signature:**
    ```cpp
    // [[Rcpp::export]]
    Eigen::SparseMatrix<float> // Or double, see Q_Rcpp_1
    hrbf_atoms_rcpp(
        const Eigen::Map<Eigen::MatrixXf> mask_xyz_world,    // N_mask_vox x 3 (float)
        const Eigen::Map<Eigen::MatrixXf> centres_xyz_world, // K_atoms x 3 (float)
        const Eigen::Map<Eigen::VectorXf> sigma_vec_mm,      // K_atoms (float)
        std::string kernel_type,                             // "gaussian" or "wendland_c4"
        double value_threshold = 1e-8                        // To sparsify output
    ) {
        const int N = mask_xyz_world.rows();
        const int K = centres_xyz_world.rows();
        if (sigma_vec_mm.size() != K) {
            Rcpp::stop("sigma_vec_mm length must match centres rows");
        }
        const bool gaussian = (kernel_type == "gaussian");

        std::vector<Eigen::Triplet<float>> triplets;
#ifdef _OPENMP
        int nthreads = omp_get_max_threads();
        std::vector<std::vector<Eigen::Triplet<float>>> thr_tri(nthreads);
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            auto& local = thr_tri[tid];
#pragma omp for schedule(static)
            for (int k = 0; k < K; ++k) {
                float sigma = sigma_vec_mm[k];
                float inv2 = 1.0f / (2.0f * sigma * sigma);
                float invs = 1.0f / sigma;
                for (int n = 0; n < N; ++n) {
                    float dx = mask_xyz_world(n,0) - centres_xyz_world(k,0);
                    float dy = mask_xyz_world(n,1) - centres_xyz_world(k,1);
                    float dz = mask_xyz_world(n,2) - centres_xyz_world(k,2);
                    float d2 = dx*dx + dy*dy + dz*dz;
                    float val;
                    if (gaussian) {
                        val = std::exp(-d2 * inv2);
                    } else {
                        float dist = std::sqrt(d2) * invs;
                        if (dist >= 1.0f) continue;
                        float base = 1.0f - dist;
                        val = std::pow(base,8.0f) * (32.0f*dist*dist*dist +
                                                   25.0f*dist*dist + 8.0f*dist + 1.0f);
                    }
                    if (std::fabs(val) > value_threshold)
                        local.emplace_back(k, n, val);
                }
            }
        }
        for (auto& vec : thr_tri) {
            triplets.insert(triplets.end(), vec.begin(), vec.end());
        }
#else
        for (int k = 0; k < K; ++k) {
            float sigma = sigma_vec_mm[k];
            float inv2 = 1.0f / (2.0f * sigma * sigma);
            float invs = 1.0f / sigma;
            for (int n = 0; n < N; ++n) {
                float dx = mask_xyz_world(n,0) - centres_xyz_world(k,0);
                float dy = mask_xyz_world(n,1) - centres_xyz_world(k,1);
                float dz = mask_xyz_world(n,2) - centres_xyz_world(k,2);
                float d2 = dx*dx + dy*dy + dz*dz;
                float val;
                if (gaussian) {
                    val = std::exp(-d2 * inv2);
                } else {
                    float dist = std::sqrt(d2) * invs;
                    if (dist >= 1.0f) continue;
                    float base = 1.0f - dist;
                    val = std::pow(base,8.0f) * (32.0f*dist*dist*dist +
                                               25.0f*dist*dist + 8.0f*dist + 1.0f);
                }
                if (std::fabs(val) > value_threshold)
                    triplets.emplace_back(k, n, val);
            }
        }
#endif
        Eigen::SparseMatrix<float> out(K, N);
        if (!triplets.empty())
            out.setFromTriplets(triplets.begin(), triplets.end());
        return out;
    }
    ```
*   **Key Implementation Details:**
    *   **Data Type:** Prefer `float` (`MatrixXf`, `VectorXf`, `SparseMatrix<float>`) if LNA targets `float32` for HDF5 storage to avoid R-side downcasting. If LNA prefers `double` for basis matrices, use `double`.
    *   **Input Mapping:** Use `Eigen::Map` to directly use memory from R matrices without copying.
    *   **Triplet Collection:** `std::vector<Eigen::Triplet<float>> triplets; triplets.reserve(...);`
    *   **Kernel Evaluation:**
        *   Gaussian: `std::exp(-d_sq / (2.0f * sigma * sigma))`
        *   Wendland C⁴: `r = dist / sigma; base = std::max(0.0f, 1.0f - r); phi = pow(base, 8) * (32*r*r*r + 25*r*r + 8*r + 1); if (r >= 1) phi = 0;`
    *   **Normalization:** L2 normalization of each atom (row of the resulting sparse matrix) *over the mask voxels it covers* should be done here if `normalize_over_mask` was part of the R-side `generate_hrbf_atom` logic. If whitening is separate, Rcpp generates unnormalized atoms. *Decision: For simplicity in Rcpp, generate unnormalized atoms. R wrapper normalizes/whitens.*
    *   **OpenMP:** Consider `#pragma omp parallel for private(thread_local_triplets)` for the outer loop over atoms (`k`), then merge triplet lists.

#### 2.2. P2: `omp_encode_cpp` (OMP Solver for `basis.empirical_hrbf_compressed`)

*   **File:** `src/lna_omp_rcpp.cpp` (or same as above)
*   **Dependencies:** `// [[Rcpp::depends(RcppEigen)]]`
*   **Signature:**
    ```cpp
    // [[Rcpp::export]]
    Rcpp::List omp_encode_rcpp(
        const Eigen::Map<Eigen::VectorXd> signal_y,            // Signal to encode (double)
        const Eigen::Map<Eigen::SparseMatrix<double>> dict_D,  // Dictionary D (N_vox x K_hrbf_atoms) (double)
        double residual_norm_sq_tol,                           // Squared norm tolerance for residual
        int max_active_atoms_L                                 // Sparsity limit
    ) {
        Eigen::VectorXd residual = signal_y;
        const Eigen::SparseMatrix<double> Dt = dict_D.transpose();
        std::vector<int> active;
        std::vector<double> coeff_vec;
        Eigen::VectorXd proj;

        for (int iter = 0; iter < max_active_atoms_L; ++iter) {
            proj = Dt * residual;
            int best = -1;
            double best_val = 0.0;
            for (int k = 0; k < proj.size(); ++k) {
                if (std::find(active.begin(), active.end(), k) == active.end()) {
                    double cand = std::abs(proj[k]);
                    if (cand > best_val) {
                        best_val = cand;
                        best = k;
                    }
                }
            }
            if (best < 0)
                break;
            active.push_back(best);

            Eigen::SparseMatrix<double> Dsub(dict_D.rows(), active.size());
            std::vector<Eigen::Triplet<double>> trip;
            for (size_t i = 0; i < active.size(); ++i) {
                int col = active[i];
                for (Eigen::SparseMatrix<double>::InnerIterator it(dict_D, col); it; ++it) {
                    trip.emplace_back(it.row(), i, it.value());
                }
            }
            Dsub.setFromTriplets(trip.begin(), trip.end());

            Eigen::VectorXd rhs = Dsub.transpose() * signal_y;
            Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(Dsub.transpose() * Dsub);
            Eigen::VectorXd c = solver.solve(rhs);

            residual = signal_y - Dsub * c;
            coeff_vec.assign(c.data(), c.data() + c.size());
            if (residual.squaredNorm() <= residual_norm_sq_tol)
                break;
        }

        bool max_flag = (active.size() == static_cast<size_t>(max_active_atoms_L) &&
                          residual.squaredNorm() > residual_norm_sq_tol);

        return Rcpp::List::create(Rcpp::Named("indices_0based") = active,
                                  Rcpp::Named("coefficients") = coeff_vec,
                                  Rcpp::Named("max_iter_reached") = max_flag);
    }
    ```
*   **Key Implementation Details:**
    *   **Data Type:** Use `double` for numerical stability in OMP.
    *   **Active Set (`active_indices_0based_std_vector`):** Store 0-based indices of selected dictionary atoms.
    *   **Coefficients (`final_coeffs_std_vector`):** Store coefficients corresponding to `active_indices_0based_std_vector`. Has length `active_indices_0based_std_vector.size()`.
    *   **$D_{sub}$ Construction:** When forming `Dsub` from `dict_D` using `active_indices_0based_std_vector`, efficiently copy the selected columns.
        ```cpp
        // Efficient Dsub construction
        SparseMatrix<double> Dsub(dict_D.rows(), active_indices_0based_std_vector.size());
        std::vector<Triplet<double>> Dsub_triplets;
        for (size_t i = 0; i < active_indices_0based_std_vector.size(); ++i) {
            int col_idx_in_D = active_indices_0based_std_vector[i];
            for (SparseMatrix<double>::InnerIterator it(dict_D, col_idx_in_D); it; ++it) {
                Dsub_triplets.emplace_back(it.row(), i, it.value());
            }
        }
        Dsub.setFromTriplets(Dsub_triplets.begin(), Dsub_triplets.end());
        ```
    *   **Least Squares:** Use `Eigen::SimplicialLLT` or `SimplicialLDLT` for solving $(D_{sub}^\top D_{sub}) c = D_{sub}^\top y$. These are for positive definite matrices. $D_{sub}^\top D_{sub}$ should be.
    *   **Return `max_iter_reached`:** As suggested, boolean `(iter_count == max_active_atoms_L && residual.squaredNorm() > residual_norm_sq_tol)`.

#### 2.3. P3: `poisson_disk_sample_component_rcpp`
*   (As detailed in my previous response, focusing on BFS/`std::deque`, correct column-major indexing for `mask_flat`, `RNGScope`, and optional voxel-grid hash for neighbor searching.)

#### 2.4. (Optional) Rcpp Helper: `voxel_to_world_rcpp` / `world_to_voxel_rcpp`
*   These would take `IntegerMatrix` (1-based voxel coords) or `NumericMatrix` (world_mm) and the 4x4 affine (as `NumericMatrix`) and perform the transformations. Useful if these conversions are frequent inside R loops that are *not* fully ported to Rcpp.

### 3. R Wrapper Integration Strategy

*   **Conditional Usage:**
    ```R
    # In R/hrbf_helpers.R or similar
    LNA_USE_RCPP_HRBF <- getOption("lna.hrbf.use_rcpp_helpers", TRUE) && 
                         requireNamespace("RcppEigen", quietly = TRUE) # And check if compiled functions exist

    # Example for hrbf_basis_from_params
    hrbf_basis_from_params <- function(params, mask_neurovol, h5_root = NULL) {
      # ... (generate/load C_total_world, sigma_vec_mm) ...
      # ... (get mask_coords_world from mask_neurovol) ...

      if (LNA_USE_RCPP_HRBF) {
        # Ensure inputs are float for Rcpp function expecting float
        B_sparse_K_x_Nmask <- lna_hrbf_rcpp_pkg_namespace$hrbf_atoms_rcpp( # Call via package if separate
          as.matrix(mask_coords_world), # Ensure MatrixXd map works
          as.matrix(C_total_world),
          as.numeric(sigma_vec_mm),
          params$kernel_type %||% "gaussian" 
          # pass value_threshold if added to Rcpp sig
        )
        # If Rcpp returns KxN, and canonical internal is NxK, transpose here
        # B_final_vox_k <- Matrix::t(B_sparse_K_x_Nmask) 
      } else {
        # ... (Fallback to pure R triplet collection + Matrix::sparseMatrix) ...
      }
      return(B_final_vox_k) # Or whatever canonical form is expected
    }

    # Similar conditional calls for OMP encoding and Poisson-disk sampling
    # in their respective R wrapper functions.
    ```
*   **Data Type Handling:** R wrappers ensure data passed to Rcpp (e.g., coordinates) is in the expected type (e.g., `matrix` of `float` if Rcpp expects `Eigen::Map<MatrixXf>`). The Rcpp function returns `Eigen::SparseMatrix<float>`, which `RcppEigen::wrap` converts to an R `Matrix::dgCMatrix` (or similar) automatically. R side handles any final downcasting to `float32` before HDF5 storage if needed.
*   **Error Handling:** Rcpp functions use `Rcpp::stop()` for errors. R wrappers use `tryCatch` if specific R-level error conditions (`lna_error_*`) are needed.

### 4. Build System

*   `src/lna_hrbf_rcpp.cpp`, `src/lna_omp_rcpp.cpp`, `src/lna_pds_rcpp.cpp` (or consolidated).
*   `// [[Rcpp::depends(RcppEigen)]]` if Eigen is used directly beyond RcppArmadillo's Eigen.
*   `useDynLib(neuroarchive, .registration = TRUE)` and Roxygen `@export` tags on Rcpp functions.
*   `DESCRIPTION`: `LinkingTo: Rcpp, RcppEigen` (if needed). `SystemRequirements: C++14`.

### 5. Granular Tickets for HRBF Rcpp Optimization Sprint (1 Sprint)

This sprint focuses *only* on implementing these Rcpp helpers and integrating them.

| #          | Ticket                                                                      | Description / Deliverables                                                                                                                                                                                                                                                                                                                                                                                                                                                              | Proposal Ref                |
| :--------- | :-------------------------------------------------------------------------- | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :-------------------------- |
| HRBF-Rcpp-S1-1 | **Rcpp: `hrbf_atoms_rcpp` Implementation**                                  | • Implement `hrbf_atoms_rcpp` in `src/lna_hrbf_rcpp.cpp` using `Eigen::SparseMatrix<float>` and triplets. <br>  • Supports "gaussian" and "wendland_c4" kernels. <br>  • Includes `value_threshold` for output sparsity. <br>  • R wrapper correctly prepares `float` matrices for `Eigen::Map<MatrixXf>`. <br> • Output matrix is $K_{atoms} \times N_{mask\_voxels}$.                                                                                                            | This Review §2.1            |
| HRBF-Rcpp-S1-2 | **Rcpp: `omp_encode_rcpp` Implementation**                                  | • Implement `omp_encode_rcpp` in `src/lna_omp_rcpp.cpp` using `Eigen`. <br>  • Uses `SimplicialLLT` for LS step. <br>  • Correctly manages active set and coefficient updates. <br>  • Returns `Rcpp::List` with 0-based `indices_0based`, `coefficients`, and `max_iter_reached` flag.                                                                                                                                                                                             | This Review §2.2, Sugg. 2   |
| HRBF-Rcpp-S1-3 | **R: Integrate `hrbf_atoms_rcpp` into `hrbf_basis_from_params`**          | • Modify R function `hrbf_basis_from_params`. <br>  - Add conditional call to `hrbf_atoms_rcpp` if `getOption("lna.hrbf.use_rcpp_helpers", TRUE)`. <br>  - Fallback to existing R-level triplet collection if Rcpp disabled. <br>  - Ensure returned basis matrix has canonical orientation (e.g., Voxels x Atoms) after Rcpp call.                                                                                                                                    | This Review §3, §2.1        |
| HRBF-Rcpp-S1-4 | **R: Integrate `omp_encode_rcpp` for Empirical Basis Compression**        | • Modify `forward_step.basis.empirical_hrbf_compressed` (or its internal OMP helper). <br>  - Add conditional call to `omp_encode_rcpp`. <br>  - R wrapper handles converting 0-based indices from Rcpp to 1-based for R, and constructing sparse output if OMP returns dense coefficient vector for active set. <br>  - Store/warn based on `max_iter_reached` flag. <br>  - Fallback to pure R OMP (if one exists) or error if Rcpp disabled and no R fallback. | This Review §3, §2.2        |
| HRBF-Rcpp-S1-5 | **(From previous) Rcpp: `poisson_disk_sample_component_rcpp`**            | • Ensure P3 (Poisson-disk sampler `poisson_disk_sample_component_rcpp` with BFS/deque, column-major indexing, RNGScope, optional grid hash) is implemented and integrated into R `lna:::poisson_disk_sample_neuroim2`.                                                                                                                                                                                                                            | Previous Review             |
| HRBF-Rcpp-S1-6 | **(Optional) Rcpp: Coordinate Transformation Helpers**                    | • Implement `voxel_to_world_rcpp` and `world_to_voxel_rcpp` if benchmarks show significant time spent in R-level `neuroim2` calls within HRBF generation loops *not already ported to Rcpp*.                                                                                                                                                                                                                                                | This Review Sugg. 1         |
| HRBF-Rcpp-S1-7 | **Build System & NAMESPACE for Rcpp**                                     | • Add Rcpp files to `src/`. Update `Makevars` if needed. <br> • Add `LinkingTo: Rcpp, RcppEigen`. Add `SystemRequirements: C++14`. <br> • Ensure `useDynLib` and Rcpp exports are correct.                                                                                                                                                                                                                                                                | Addendum §6                 |
| HRBF-Rcpp-S1-8 | **Unit Tests for Rcpp Functions & R Wrappers**                            | • Direct unit tests for `hrbf_atoms_rcpp` and `omp_encode_rcpp` (e.g., call from R with simple inputs, check outputs). <br> • Test R wrappers (`hrbf_basis_from_params`, OMP part of empirical compression) with Rcpp enabled vs. disabled (fallback path) – check for numerical equivalence (within tolerance). <br> • Add performance smoke tests for key Rcpp functions.                                                                                               | Addendum §5                 |
| HRBF-Rcpp-S1-9 | **Documentation for Rcpp Usage**                                          | • Note in relevant HRBF appendix sections and developer guides that core computations are Rcpp-accelerated. <br> • Document `getOption("lna.hrbf.use_rcpp_helpers", TRUE)`.                                                                                                                                                                                                                                                                | -                           |

This sprint aims to deliver a significant performance boost to the HRBF-related transforms by porting the most critical computational kernels to Rcpp/Eigen, making them practical for whole-brain analyses.