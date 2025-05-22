Okay, this final feedback targets crucial details for implementation robustness and clarity. Here is the v1.4 specification for the `lna` R package, incorporating these last refinements.

---

**LNA R Package Implementation Specification (v1.4 - Implementation Ready)**

**(Summary of Changes from v1.3):** This definitive version incorporates final clarifications based on detailed review. Key refinements include: specifying the exact parameter default merging order, clarifying checksum computation timing, defining non-interactive behavior for plugin prompting, specifying error classes, noting fork-safety issues with schema caches, adding warnings for namespace collisions, enhancing HDF5 robustness retry logic for chunk sizes, detailing lazy reader closure idempotency, improving scaffolded code, and clarifying file opening modes for parallel writing. The specification is now considered complete and ready for implementation.

**1. Overview & Purpose**


(No changes from v1.0 - goals remain the same)

### Data shape convention

Arrays supplied to `write_lna()` should store time in the final
dimension.  Transforms operate on matrices with time along rows and
voxels along columns.  Thus a `10×4×1` array is first reshaped to a
`10×4` matrix before the Sparse PCA step.

**2. Public API**

```R
#' Write data to an LNA file
#'
#' Compresses and writes neuroimaging object(s) to an LNA 2.0 HDF5 file.
#'
#' @param x A neuroimaging object (e.g., `LatentNeuroVec`, 4D array) or a list
#'   of such objects for multiple runs. If an unnamed list, runs are assigned
#'   deterministic names `run-01`, `run-02`, ... based on positional order (`seq_along(x)`).
#' @param file Character string: path to the output `.lna.h5` file. If `NULL`,
#'   attempts in-memory HDF5 writing via `hdf5r` driver='core'.
#' @param transforms Character vector: ordered sequence of transform types (forward order).
#' @param transform_params Named list: Parameters for each transform type. Names
#'   must match `transforms`. Parameter resolution uses a deep merge (left-to-right precedence):
#'   1. Transform schema defaults (via `lna:::default_params`).
#'   2. Package-level defaults (via `lna_options()`).
#'   3. User-supplied list in this argument.
#' @param mask Optional: `LogicalNeuroVol` or 3D logical array. Error if `x` provided
#'   and `sum(mask)` doesn't match `x`'s spatial voxel count.
#' @param header Optional: Named list of NIfTI-like header attributes.
#' @param block_table Optional: Data frame for `/spatial/block_table`. Coordinates
#'   (x0, x1, etc.) MUST be 1-based inclusive indices in masked voxel space.
#' @param plugins Optional: Named list for `/plugins/`.
#' @param compression_level Integer (0-9): Gzip compression level. Falls back to 0
#'   with a warning if HDF5 lacks zlib support.
#' @param chunk_dims Optional: HDF5 chunk dimensions. If `NULL`, uses heuristic:
#'   Target <= 1 MiB/chunk; if est. compressed chunk > 1 GiB for >4 GiB data,
#'   halve fastest axis until < 1 GiB. Auto-reduced further if needed for HDF5 limits (e.g., <= 256 MiB target if initial reduction fails), with warnings.
#' @param checksum Character: `"none"` (default) or `"sha256"`. If `"sha256"`, the hash
#'   of the *entire file byte stream* is computed using `digest::digest(file=...)`
#'   *after* the HDF5 file handle is closed. The hash is then stored in the root
#'   attribute `lna_checksum` (requires reopening briefly or external tooling).
#'   The core writer MUST NOT reopen the file after hashing.
#' @param verbose Logical: Print progress? Uses `progressr` if available and handlers configured.
#'
#' @return Invisibly returns a list: `file` (path/NULL), `plan` (Plan R6), `header` (written list).
#' @export
write_lna <- function(...) { ... } # Implementation uses core_write

#' Read data from an LNA file
#'
#' Reads LNA 2.0 files, applying inverse transforms and optional subsetting.
#'
#' @param file Character string: path to input `.lna.h5` file.
#' @param run_id Character vector or string: Glob pattern or specific run ID(s).
#' @param roi_mask Optional: `LogicalNeuroVol` or 3D logical array...
#' @param time_idx Optional: Integer vector or slice...
#' @param as_latent Logical: Return `LatentNeuroVec` (TRUE) or reconstructed array/volume (FALSE)?
#' @param allow_plugins Character: `"installed"` (default): Load if package installed.
#'   `"prompt"`: Ask interactively (falls back to `"installed"` behavior if `!rlang::is_interactive()`).
#'   `"none"`: Skip optional transforms requiring external packages.
#' @param validate Logical: Perform basic runtime validation?
#' @param output_dtype Character: Target data type: `"float32"` (default), `"float64"`.
#'   Requesting `"float16"` raises an error of class `"lna_error_float16_unsupported"`
#'   unless `lna:::has_float16_support()` returns `TRUE`.
#' @param lazy Logical: Immediate reconstruction (`FALSE`, default) or return `lna_reader`
#'   proxy (`TRUE`). `lna_reader` keeps HDF5 file open. Call `$close()` when done;
#'   GC finalizers (`$finalize()`/`reg.finalizer`) provide backup closure.
#' @param verbose Logical: Print progress? Uses `progressr` if available and configured.
#'
#' @return If `lazy=FALSE`: `LatentNeuroVec`, 4D array/`NeuroVol`, or list.
#'   If `lazy=TRUE`: An object of class `lna_reader`.
#' @export
read_lna <- function(...) { ... } # Implementation uses core_read

# Convenience aliases/wrappers
#' @export
compress_fmri <- function(...) write_lna(...)
#' @export
open_lna <- read_lna
#' @export
validate_lna <- function(file, strict = TRUE, checksum = TRUE) { ... } # Calls internal full validator

# Optional helper for global settings
#' Get or set global lna options
#'
#' @param ... Options to set (`write.compression_level=4`, etc.) or character names
#'   of options to retrieve. If empty, returns all current options.
#' @return Invisibly returns updated list of all options, or named list of requested options.
#' @export
lna_options <- function(...) { # Uses internal package environment }
```

**3. Core Internal Architecture**

*   **`write_lna` Core (`core_write.R`):**
    *   Handles default run naming.
    *   Resolves `transform_params` using specified merge order.
    *   Orchestrates forward pass via `forward_step`.
    *   Calls `materialise_plan`.
*   **`read_lna` Core (`core_read.R`):**
    *   Handles `allow_plugins` modes and non-interactive fallback.
    *   Checks `output_dtype == "float16"` requirements.
    *   Orchestrates reverse pass via `invert_step`.
    *   Handles `lazy=TRUE` return.
*   **S3 Dispatch (`dispatch.R`):** Generics `forward_step(type, desc, handle)` and `invert_step(type, desc, handle)`.

**4. Key Components Detailed Specification**

*   **`DataHandle` (`handle.R` - R6 Class):**
    *   **Public Fields:** (As before)
    *   **Public Methods:** (As before)
        *   `update_stash()`: Note for implementer: Must return the *new* `DataHandle` instance created by `$with()`, not modify `self`.

*   **`Plan` (`plan.R` - R6 Class):**
    *   **Public Fields:**
        *   `datasets` (tibble): Columns: `path`, `role`, `producer`, `origin`, `step_index`, `params_json`, `payload_key`, `write_mode` ("eager"/"stream"), `write_mode_effective` (string: actual mode used after fallback, "eager" or potentially "stream" in future).
        *   `descriptors`, `payloads`, `next_index`: (As before)
    *   **Public Methods:** (As before)
        *   `add_dataset_def()`: Now also records `step_index`. `write_mode_effective` is set by `materialise_plan`.
        *   `mark_payload_written()`: (As before)

*   **`lna_reader` (`reader.R` - S3 or R6 Class):**
    *   **Internal State:** (As before)
    *   **Methods:**
        *   `$subset(...)`: (As before)
        *   `$data(...)`: (As before - idempotent/cached result recommended)
        *   `$close()`: Closes HDF5 handle. MUST be idempotent (safe to call multiple times, e.g., in `tryCatch` finalizers).
        *   `$finalize()` / `reg.finalizer`: (As before - provide backup closure)
        *   `print()`: (As before)

**5. Transform Implementation Guide**

(As before, implementing S3 methods)
*   Use `lna:::default_params(type)` to retrieve schema defaults.
*   Be aware `write_mode = "stream"` falls back to "eager" in v1.4.
*   Implement subsetting logic.

**6. Validation Strategy**

*   **Runtime:** (As before)
*   **Full Audit (`validate_lna`):**
    *   (As before + checksum check if requested)
    *   Schema Cache Note: Compiled validators (e.g., from `jsonvalidate`) might not be fork-safe. If using `future::plan(multicore)` heavily involving validation, consider calling `lna:::schema_cache_clear()` within each worker's setup, or use a fork-safe caching strategy if available.
*   **Error Reporting:** Use `rlang::abort(..., location = ...)` with step index/name.

**7. Helper Utilities (`utils_*.R`, `Rcpp`)**

*   (List as before)
*   **`lna:::materialise_plan(...)` Updates:**
    *   Sets `write_mode_effective` in the `plan` based on actual write behavior (e.g., after fallback). Uses throttled warning for fallback.
    *   Implements checksum logic (write hash attr).
    *   Handles HDF5 errors:
        *   Retry w/o compression on filter errors + warning.
        *   Retry with smaller chunks on write errors + warning (first heuristic: target < 1 GiB; second heuristic if still fails: target <= 256 MiB).
    *   Uses `step_index` for provenance in error messages.
*   **New/Updated Helpers:**
    *   `lna:::default_params(type)`: Reads schema, extracts defaults.
    *   `lna:::has_float16_support()`: Checks dependencies.
    *   `lna::check_transform_implementation(type)`: Now also warns if `type` namespace collides with core transforms (`quant`, `basis`, `embed`, `temporal`, `delta`) or base R packages (stats, utils, etc.).
    *   `lna::scaffold_transform(type)`: Generates template files, including a stub using `lna:::default_params()`.
    *   `lna:::schema_cache_clear()`: Exposes cache clearing for tests.
*   **Progress Reporting:** Check `!progressr::handlers_is_empty()` before invoking `progressr`.

**8. Package Structure**

(As before, including `reader.R`)
*   **Testing:** Use `driver='core'`, test error handling, multi-run, checksums, lazy reader, schema validation, default parameter merging, HDF5 robustness fallbacks.

**9. Concurrency**

*   Core read/write is safe for basic parallelism.
*   **Parallel Writing Note:** Multiple concurrent `write_lna` calls SHOULD use unique temporary filenames within the target directory and perform an atomic `file.rename()` upon successful completion. Writers implicitly open the target file with truncation enabled (standard HDF5 'w' mode behavior); concurrent writes directly to the *same final path* results in undefined behavior.

**10. Documentation & Examples**

*   Package documentation should clearly explain the parameter default merging order, checksum scope, lazy reader lifecycle, and plugin options.
*   A "Cookbook" vignette is strongly recommended, demonstrating:
    *   Basic compression (`compress_fmri`).
    *   ROI/time slicing with the lazy reader (`read_lna(lazy=TRUE)`).
    *   Creating and using a simple custom transform via `scaffold_transform()`.

**11. Future Considerations**

*   Full "stream" write implementation.
*   GPU acceleration integration.
*   Advanced HDF5 filters (Blosc, Zstd).
*   Enhanced `lna_reader` API.

This v1.4 specification is considered final and ready for implementation. It addresses all identified ambiguities and edge cases, providing a solid foundation for a robust and extensible LNA package in R.

Okay, let's integrate the concept of aggregating data across runs *before* applying an encoding like Sparse PCA.

This involves adding a conceptual **aggregation transform** that runs *first* on the list of input runs, producing a single data structure that the subsequent `myorg.sparsepca` transform will operate on.

Here's the revised recipe card incorporating this idea, formatted in Markdown:

---

# Recipe: Adding Transforms to the LNA Ecosystem

This guide demonstrates how to integrate custom transforms into the LNA 2.0 format and the `lna` R package ecosystem (v1.4 spec). We'll use two related examples:

1.  `myorg.aggregate_runs`: Aggregates data across multiple input runs.
2.  `myorg.sparsepca`: Applies Sparse PCA, operating either on single-run data *or* on the output of the aggregation step.

This makes your transforms first-class citizens, usable alongside core transforms.

---

## Part 1: `myorg.aggregate_runs` (Conceptual)

**Goal:** Combine data from multiple fMRI runs (e.g., by concatenating time series) to create a single dataset suitable for global encoding (like PCA across all data).

### 1.1 Pick Namespace & Version

*   **`type`**: `"myorg.aggregate_runs"`
*   **`version`**: `"1.0"`

### 1.2 Define Parameters (via Schema)

This transform might need parameters like the aggregation method.

**File:** `inst/schemas/myorg.aggregate_runs.schema.json` (Example)
```json
{
  "$schema": "http://json-schema.org/draft-07/schema#",
  "$id": "https://my-lab.org/schemas/lna/2.0/myorg.aggregate_runs.schema.json",
  "title": "Parameters for 'myorg.aggregate_runs' transform",
  "description": "Aggregates data across multiple input runs.",
  "type": "object",
  "properties": {
    "method": {
      "enum": ["concatenate_time", "average_covariance"],
      "default": "concatenate_time",
      "description": "Method used to combine data from multiple runs."
    },
    "runs_included": {
        "type": "array",
        "items": {"type": "string"},
        "description": "List of run_ids that were aggregated (recorded by writer)."
    }
  },
  "required": ["method"],
  "additionalProperties": false
}
```

### 1.3 Implement S3 Methods

*   **`forward_step.myorg.aggregate_runs`**:
    *   **Input:** Expects the *initial list* of run objects provided to `write_lna`. Convention: `handle$stash$initial_input_list`.
    *   **Logic:** Based on `desc$params$method`, iterates through the input list, performs aggregation (e.g., concatenates matrices along the time dimension). Records which runs were included in `desc$params$runs_included`.
    *   **Output:** Places the single aggregated data structure (e.g., a large matrix `aggregated_matrix`) into the `stash`.
    *   **Plan:** Adds its JSON descriptor to `handle$plan`. *Typically does not add numeric payloads itself*, as it just prepares data for the next step. Updates `desc$inputs/outputs` accordingly (e.g., `inputs=["initial_input_list"]`, `outputs=["aggregated_matrix"]`).
    *   Returns the updated `handle`.

*   **`invert_step.myorg.aggregate_runs`**:
    *   **Input:** Receives the reconstructed *aggregated* data from the *next* inverse step (e.g., `inputs=["aggregated_matrix_hat"]`).
    *   **Logic:** For basic use cases where the goal is the aggregated result, this might be a no-op data-wise. It could potentially add metadata back indicating the aggregation source runs (from `desc$params$runs_included`). Reconstructing individual runs from the aggregated inverse would require much more complex logic and likely storing extra information during the forward pass.
    *   **Output:** Passes the aggregated data through, perhaps under a different name if needed (e.g., `outputs=["final_aggregated_result"]`).
    *   Returns the updated `handle`.

*(Detailed code implementation omitted for brevity, as it depends heavily on the specific aggregation logic.)*

---

## Part 2: `myorg.sparsepca` (Operating on Aggregated or Single-Run Data)

**Goal:** Apply Sparse PCA to input data. This transform can now operate either on single-run data *or* the aggregated output from `myorg.aggregate_runs`.

### 2.1 Pick Namespace & Version

*   **`type`**: `"myorg.sparsepca"`
*   **`version`**: `"1.0"`

### 2.2 Write the JSON Schema for Parameters

*(Same schema as before - defines `k`, `alpha`, `whiten`, `storage_order`)*

**File:** `inst/schemas/myorg.sparsepca.schema.json`
```json
{
  "$schema": "http://json-schema.org/draft-07/schema#",
  "$id": "https://my-lab.org/schemas/lna/2.0/myorg.sparsepca.schema.json",
  "title": "Parameters for 'myorg.sparsepca' transform",
  // ... (properties: k, alpha, whiten, storage_order) ... see previous example
  "required": ["k"],
  "additionalProperties": false
}
```

### 2.3 Create S3 Methods

These methods now need to be aware of the *expected input key* which changes depending on whether aggregation occurred.

```R
# In your package's R code (e.g., R/transform_sparsepca.R)

#' Forward Sparse PCA Transform Step
#' @export
#' @importFrom Matrix Matrix t tcrossprod
#' @importFrom jsonlite toJSON
#' @keywords internal
forward_step.myorg.sparsepca <- function(type, desc, handle) {

  p <- desc$params # Merged parameters

  # --- 1. Get Input (Dynamically) & Fit Sparse PCA ---
  # Determine input key based on expected upstream transform
  # Convention: Aggregation step outputs "aggregated_matrix", default input is "dense_mat"
  input_key <- if (handle$has_key("aggregated_matrix")) "aggregated_matrix" else "dense_mat"
  inputs <- handle$get_inputs(input_key)
  X <- inputs[[input_key]] # nTime (or nTotalTime) x nVox matrix

  # --- Fit Sparse PCA (as before) ---
  # fit <- my_sparse_pca_solver(X, k = p$k, alpha = p$alpha, whiten = p$whiten)
  # basis_mat <- ... # k x nVox or nVox x k (respect p$storage_order)
  # embed_mat <- ... # nTime x k
  # Placeholder implementation:
  stopifnot(is.numeric(p$k), p$k > 0)
  nVox <- ncol(X); nTime <- nrow(X)
  basis_mat <- Matrix::Matrix(rnorm(p$k * nVox), nrow = p$k, ncol = nVox)
  if (p$storage_order == "voxel_x_component") basis_mat <- Matrix::t(basis_mat)
  embed_mat <- matrix(rnorm(nTime * p$k), nrow = nTime, ncol = p$k)

  # --- 2. Register Datasets in the Plan ---
  # Define standard HDF5 paths
  b_path <- "/basis/global" # Assume this writes THE global basis
  e_path <- if (input_key == "aggregated_matrix") "/scans/aggregated/embedding" else "/scans/run-CURRENT/embedding" # Path depends on input context! Needs refinement based on multi-run handling in core.
  # Simplified: Assume we *always* write to a canonical path for this example.
  e_path <- "/scans/derived/embedding" # Placeholder path

  handle$plan$add_payload(b_path, basis_mat)
  handle$plan$add_payload(e_path, embed_mat)

  handle$plan$add_dataset_def(
        path = b_path, role = "basis_matrix", producer = type,
        origin = handle$plan$origin_label, step_index = handle$plan$next_index -1,
        params_json = jsonlite::toJSON(p, auto_unbox = TRUE),
        payload_key = b_path, write_mode = "eager")

  handle$plan$add_dataset_def(
        path = e_path, role = "coefficients", producer = type,
        origin = handle$plan$origin_label, step_index = handle$plan$next_index -1,
        params_json = jsonlite::toJSON(p, auto_unbox = TRUE),
        payload_key = e_path, write_mode = "eager")

  # --- 3. Fill in Descriptor ---
  desc$version  <- "1.0"
  desc$datasets <- list( list(path = b_path, role = "basis_matrix"),
                         list(path = e_path, role = "coefficients") )
  desc$inputs   <- c(input_key) # Consume the correct input key
  desc$outputs  <- c("sparsepca_basis", "sparsepca_embedding") # More specific output names

  handle$plan$add_descriptor(
        transform_name = handle$plan$get_next_filename(type), desc_list = desc)

  # --- 4. Update Stash for Next Forward Step ---
  out <- list(sparsepca_basis = TRUE, sparsepca_embedding = TRUE) # Signal presence of outputs
  handle <- handle$update_stash(desc$outputs, out) # Removes input_key implicitly

  handle # Return updated handle
}

#' Inverse Sparse PCA Transform Step
#' @export
#' @importFrom Matrix Matrix t tcrossprod
#' @keywords internal
invert_step.myorg.sparsepca <- function(type, desc, handle) {

  # --- 1. Load Datasets ---
  basis_path <- desc$datasets[[which(sapply(desc$datasets, function(d) d$role == "basis_matrix"))]]$path
  embed_path <- desc$datasets[[which(sapply(desc$datasets, function(d) d$role == "coefficients"))]]$path

  basis <- handle$h5[[basis_path]]$read()
  embed <- handle$h5[[embed_path]]$read()

  # --- 2. Apply Subsetting ---
  p <- desc$params
  storage_order <- p$storage_order %||% "component_x_voxel"
  # (Subsetting logic as before, potentially complex for ROI on basis)
  if (!is.null(handle$subset$roi)) { ... }
  if (!is.null(handle$subset$time)) {
      embed <- embed[handle$subset$time, , drop = FALSE]
  }

  # --- 3. Reconstruct Data ---
  if (storage_order == "voxel_x_component") { basis <- Matrix::t(basis) }
  X_hat <- Matrix::tcrossprod(embed, basis) # (nTime_subset x nVox_subset)

  # --- 4. Update Stash ---
  # Output key depends on what the *next* inverse step expects.
  # If this was preceded by aggregate_runs, maybe output "aggregated_matrix_hat"
  # If it operated on single run, maybe output "dense_mat"
  # This highlights the importance of clear contracts defined by inputs/outputs!
  # Assume for this example, the next step (inverse aggregate or inverse quant)
  # expects the reconstructed matrix under a specific name.
  output_key <- desc$outputs[[1]] # Infer from the defined outputs in the descriptor
  outlist <- list(X_hat); names(outlist) <- output_key

  handle <- handle$update_stash(desc$outputs, outlist)

  # --- 5. Metadata ---
  if (!is.null(desc$metadata_updates)) { handle$meta <- utils::modifyList(handle$meta, desc$metadata_updates) }

  handle
}
```

### 2.4 Ship Default Parameter Helper (Optional)

```R
#' Default parameters for myorg.sparsepca
#' @export
#' @keywords internal
lna_default.myorg.sparsepca <- function() {
  # Reads defaults directly from schema if lna:::default_params is robust enough,
  # or specify manually:
  list(k = 50L, alpha = 1e-3, whiten = FALSE, storage_order = "component_x_voxel")
}
```

### 2.5 Declare Dependencies & Install Schema

*(Same as before: Add `lna`, `irlba`, `Matrix`, `jsonlite`, `rlang` to `Imports`. Install schema in `inst/schemas/`)*

---

## Part 3: Using the Combined Transforms

Now, users can choose to run `sparsepca` alone or combine it with aggregation.

```R
library(lna)
library(myorgLNAExtensions) # Package containing both transforms

# --- Scenario 1: Sparse PCA on a single run ---
write_lna(
   x                = single_run_array, # 4D array
   file             = "sub-01_task-motor_spca.lna.h5",
   mask             = single_run_mask,
   transforms       = c("myorg.sparsepca", "quant"), # Apply SPCA, then quantize results
   transform_params = list(
        `myorg.sparsepca` = list(k = 80),
        quant = list(bits = 8)
   )
)
# read_lna will return the reconstructed single run data

# --- Scenario 2: Aggregate runs, then Sparse PCA ---
write_lna(
   x                = list(run1 = run1_array, run2 = run2_array), # List of runs
   file             = "sub-01_task-all_agg_spca.lna.h5",
   mask             = common_mask, # Must be same mask for aggregation
   transforms       = c("myorg.aggregate_runs", "myorg.sparsepca", "quant"), # Aggregate FIRST
   transform_params = list(
        `myorg.aggregate_runs` = list(method = "concatenate_time"),
        `myorg.sparsepca`      = list(k = 120, alpha = 0.005),
        quant = list(bits = 6)
   )
)

# Reading the aggregated file:
agg_dat <- read_lna("sub-01_task-all_agg_spca.lna.h5")
# agg_dat will contain the reconstructed *aggregated* data, not individual runs,
# because the simple inverse steps don't perform disaggregation.
```

## Part 4: Validate

*(Same as before - `validate_lna` works on files created with custom transforms)*

---

This combined recipe illustrates how the modular transform system allows complex workflows, like cross-run aggregation followed by encoding, to be built by composing independent, well-defined transform steps. The key is careful definition of the `inputs` and `outputs` contract for each step

Addendum: Okay, here's an addendum to the recipe card, incorporating the detailed considerations and clarifications based on the LNA v1.4 specification.

---

## Addendum: Implementation Refinements & v1.4 Alignment

This addendum provides further details and clarifications for implementing the `myorg.aggregate_runs` and `myorg.sparsepca` transforms, ensuring closer alignment with the LNA v1.4 specification and promoting robust inter-transform communication via the `DataHandle` stash.

**1. Input Source for `myorg.aggregate_runs` (First Transform Convention):**

*   **LNA v1.4 Implication:** The `myorg.aggregate_runs` transform, if used, *must* be the first transform in the `transforms` sequence specified in `write_lna`.
*   **`core_write` Responsibility:**
    *   When `write_lna` is called with a list of run objects (`x`), `core_write` should detect if the first transform is an aggregator (e.g., by a convention or a declared property of the transform).
    *   If `transforms[1]` is `myorg.aggregate_runs`, `core_write` will place the entire input list `x` into `handle$stash` under a well-defined key (e.g., `initial_input_list`). The `forward_step.myorg.aggregate_runs` will then consume this key.
    *   If `transforms[1]` is a standard per-run transform (like `myorg.sparsepca` used directly), `core_write` will iterate over each run in `x`, populating `handle$stash` with the *current run's data* (e.g., under `dense_mat`) before calling its `forward_step`.

**2. Output Key for `invert_step.myorg.aggregate_runs`:**

*   **LNA v1.4 Implication:** The `invert_step` for any transform `T` should aim to reproduce the data that was the *input* to `forward_step(T)`.
*   When `myorg.aggregate_runs` is the first forward step, its inverse is the last inverse step.
*   The data it outputs from its `invert_step` will be the final, reconstructed (aggregated) data that the user expects from `read_lna`. It should not introduce a new, arbitrary key name at this final stage unless it's explicitly part of its documented behavior (e.g., adding metadata alongside the data).

**3. Dynamic HDF5 Path Naming for `myorg.sparsepca` Embeddings (`e_path`):**

*   **LNA v1.4 Implication:** For flexibility and clarity, especially when `myorg.sparsepca` might operate on individual runs (without prior aggregation):
    *   If `myorg.sparsepca` processes individual runs, `core_write` should make the current `run_id` (e.g., `run-01`, `run-02` as determined from `names(x)` or `run-XX` convention) available to `forward_step` (e.g., via `handle$current_run_id` or as an argument).
    *   `forward_step.myorg.sparsepca` would then use this `run_id` to construct unique HDF5 paths, e.g., `/scans/${run_id}/embedding`.
    *   When `myorg.aggregate_runs` *is* used, it effectively consumes all run-specific identities for subsequent steps. `myorg.sparsepca` then operates on this single aggregated entity. In this scenario, a non-run-specific path like `/scans/aggregated/embedding` or `/scans/derived/embedding` is appropriate, as used in the recipe. The `Plan`'s `origin_label` can also help contextualize these datasets.

**4. Transform Parameter Defaulting Strategy:**

*   **LNA v1.4 Mechanism:** The primary mechanism for defaults is `lna:::default_params(type)`, which reads `default` fields from the transform's JSON schema.
*   **Custom Transform Packages:**
    *   A transform package can provide `lna_default.MY_TYPE <- function() { lna:::default_params("MY_TYPE") }` to use schema defaults directly.
    *   Alternatively, as shown in the recipe (`lna_default.myorg.sparsepca <- function() { list(...) }`), a transform can provide a function that returns a manually specified list of defaults. This is useful if defaults are more complex than what JSON schema `default` fields can express or if no schema is used for a very simple internal transform.
    *   The merge order specified in `write_lna` (`transform_params` argument) remains: Schema defaults -> Package-level defaults (`lna_options()`) -> User-supplied list.

**5. Refined Stash Management and Key Naming Conventions:**

*   **Principle:** The `desc$outputs` from `forward_step(T_i)` should explicitly name the stash keys that `forward_step(T_{i+1})` will consume via its `desc$inputs`. The data associated with these keys in the stash should be the actual data payload or a direct representation needed by the next step.
*   **Example Flow (`agg -> spca -> quant`):**

    1.  **`core_write` Preparation:**
        *   `handle$stash` contains `initial_input_list = list_of_run_arrays`.

    2.  **`forward_step.myorg.aggregate_runs` (`desc_agg`, `handle`):**
        *   **Consumes:** `X_list <- handle$get_inputs("initial_input_list")[["initial_input_list"]]`
        *   **Computes:** `aggregated_data <- do_aggregation(X_list)`
        *   **Descriptor:**
            *   `desc_agg$inputs <- c("initial_input_list")`
            *   `desc_agg$outputs <- c("aggregated_matrix_for_spca")` (This key is for the *next* step)
        *   **Updates Stash:**
            `handle <- handle$update_stash(outputs_to_add = list(aggregated_matrix_for_spca = aggregated_data), inputs_to_remove = c("initial_input_list"))`
        *   Records `desc_agg` (with its params like `method`, `runs_included`) in `handle$plan`.

    3.  **`forward_step.myorg.sparsepca` (`desc_spca`, `handle`):**
        *   **Consumes:** `current_data <- handle$get_inputs("aggregated_matrix_for_spca")[["aggregated_matrix_for_spca"]]`
        *   **Computes:** `fit <- spca(current_data)` leading to `basis_mat` and `embed_mat`.
        *   **Writes to HDF5:** `basis_mat` (e.g., to `/basis/global_spca`), `embed_mat` (e.g., to `/scans/aggregated/spca_embedding`).
        *   **Descriptor:**
            *   `desc_spca$inputs <- c("aggregated_matrix_for_spca")`
            *   `desc_spca$outputs <- c("embedding_for_quant")` (The `embed_mat` is what `quant` will act upon)
        *   **Updates Stash:**
            `handle <- handle$update_stash(outputs_to_add = list(embedding_for_quant = embed_mat), inputs_to_remove = c("aggregated_matrix_for_spca"))`
        *   Records `desc_spca` in `handle$plan`.

    4.  **`forward_step.quant` (`desc_quant`, `handle`):**
        *   **Consumes:** `coeffs_to_quantize <- handle$get_inputs("embedding_for_quant")[["embedding_for_quant"]]`
        *   **Computes:** `quantized_coeffs <- quantize(coeffs_to_quantize)`
        *   **Writes to HDF5:** `quantized_coeffs` (and scale/offset if applicable).
        *   **Descriptor:**
            *   `desc_quant$inputs <- c("embedding_for_quant")`
            *   `desc_quant$outputs <- c()` (If it's the last data-mutating step, it might not add new data to the stash for subsequent *forward* steps).
        *   **Updates Stash:**
            `handle <- handle$update_stash(inputs_to_remove = c("embedding_for_quant"))`
        *   Records `desc_quant` in `handle$plan`.

**6. Symmetric Stash Keys for Inverse Steps:**

*   **Principle:** `invert_step(T_i)` should produce data in the stash under the key(s) that `forward_step(T_i)` *consumed* (i.e., `T_i$desc$inputs`). This ensures that `invert_step(T_{i-1})` can find its expected reconstructed input.
*   **Example `invert_step.myorg.sparsepca`:**
    *   Its forward step (`forward_step.myorg.sparsepca`) consumed `aggregated_matrix_for_spca`.
    *   Its forward step produced `embedding_for_quant` (which `invert_step.quant` will have reconstructed and placed back into the stash as `embedding_for_quant`).
    *   `invert_step.myorg.sparsepca` will:
        1.  Load its basis from HDF5.
        2.  Get `reconstructed_embedding <- handle$get_inputs("embedding_for_quant")[["embedding_for_quant"]]`.
        3.  Compute `X_hat <- tcrossprod(reconstructed_embedding, basis)`.
        4.  **Descriptor:** `desc$inputs` (from forward pass) was `c("aggregated_matrix_for_spca")`. `desc$outputs` (from forward pass) was `c("embedding_for_quant")`.
        5.  **Updates Stash:**
            `output_key_name <- desc$inputs[[1]] # "aggregated_matrix_for_spca"`
            `outlist <- list(X_hat); names(outlist) <- output_key_name`
            `handle <- handle$update_stash(outputs_to_add = outlist, inputs_to_remove = desc$outputs)` (remove `embedding_for_quant`)
*   This `X_hat` (now under `aggregated_matrix_for_spca` in stash) is then available for `invert_step.myorg.aggregate_runs`.

**7. `myorg.aggregate_runs` Forward Step - No Numeric Payloads to HDF5:**

*   The recipe correctly states: `Typically does not add numeric payloads itself, as it just prepares data for the next step.`
*   `myorg.aggregate_runs` primarily manipulates data within the `DataHandle` stash. Its primary record in the HDF5 file via the `Plan` is its descriptor JSON (containing parameters like `method` and `runs_included`), not HDF5 datasets of numerical data.

By adhering to these refined conventions, especially regarding stash key management and the roles of `core_write` versus individual transforms, the LNA ecosystem can maintain clarity and robustness as more complex transform sequences are developed.

.Below is a self-contained “Recipe Card” you can paste into your wiki or docs/ folder, followed by a short Addendum that captures caveats, pro-tips, and links back to the relevant paragraphs of the LNA v1.4 spec + coding blueprint.

⸻

📄 Recipe Card — Adding a New LNA Transform

Goal: Take any clever encoder/decoder you’ve written (Sparse PCA, Oct-tree PCA, tiny Conv-AE, …) and make it a first-class step in the LNA pipeline so write_lna() and read_lna() pick it up automatically.

Step	What you write	Where	Why
1. Pick a stable type name	e.g. "myorg.octree_pca"  (namespace + short slug)	In every descriptor & function name	Avoid clashes with core types (quant, basis, …).
2. Draft a JSON-Schema	inst/schemas/myorg.octree_pca.schema.json	Defines allowed params + defaults that lna:::default_params() will surface.	
3. Implement two S3 methods	r<br># forward (writer)<br>forward_step.myorg.octree_pca <- function(type, desc, handle) { … }<br><br># inverse (reader)<br>invert_step.myorg.octree_pca  <- function(type, desc, handle) { … }	R/transform_octree_pca.R	• Forward fits the model, pushes payloads into handle$plan.  • Inverse re-hydrates from HDF5, applies ROI/temporal subset, returns new handle.
4. Register payload paths in the Plan	handle$plan$add_payload("/basis/blocks/<id>/matrix", block_mat)	Ensures materialise_plan() knows what to write to disk later.	
5. Fill the JSON descriptor	r<br>desc$datasets <- list(…); desc$inputs <- …; desc$outputs <- …<br>handle$plan$add_descriptor(handle$plan$get_next_filename(type), desc)	Keeps the on-disk /transforms/NN_myorg.octree_pca.json self-describing.	
6. (Optional) Provide an lna_default helper	r<br>lna_default.myorg.octree_pca <- function() list(k = 32L, depth = 3L)	Lets users omit most transform_params.	
7. Add unit tests	• Round-trip write → read for a toy array. • ROI & temporal subset behave.	tests/testthat/test-octree_pca.R	Keeps CI green.
8. Declare Imports	Matrix, irlba, etc. in DESCRIPTION	So NAMESPACE is generated correctly.	

That’s it.  Drop the package under Suggests: in the main lna package and the ecosystem will auto-discover your transform when it sees a descriptor of that type.

⸻

📎 Addendum — Odds, Ends & Best-Practices

Topic	Quick guidance	Spec §
Aggregation before encoding	Chain an up-front transform (myorg.aggregate_runs) whose outputs (e.g. "aggregated_matrix") are the inputs for your encoder.  The inverse can be a no-op if you only need the pooled result.	§3, Recipe Card Pt 1
Oct-tree PCA details	• Use /spatial/block_table to store the 8-/16-/64-cube bounds and /basis/blocks/<id> for each local matrix. • Column offsets live in /spatial/coeff_offset_block.	Table in previous answer
“stream” vs “eager” writes	In v1.4 every payload is written eagerly; set write_mode = "stream" if you plan to support streaming later—materialise_plan() will silently fall back and emit a once-per-session warning.	Spec v1.4 §7
Auxiliary files (e.g. grey-matter mask)	Put them in /plugins/<your_step>/ and list them in the descriptor’s datasets with role = "aux_mask".  Your inverse reads them through handle$h5[[path]].	Spec v1.4 §2
Default-param merge order	1. Schema defaults → 2. lna_options() → 3. User transform_params.  Handle this in your forward_step.* by calling lna:::default_params(type) then utils::modifyList().	Spec v1.4 §2
Error classes	Throw with lna:::abort_lna("message", class = "lna_error_contract", location = "myorg.octree_pca") so upstream code can catch/filter.	Spec v1.4 §6
ROI / time sub-setting	Declare capabilities.supports_spatial_subsetting = "block-only" (oct-tree) and implement the fast path; fall back to full data if request outside capability.	Spec v1.4 §3.5
Schema URI	The schema_uri in the descriptor does not need to be publicly reachable in early days—just set it to the canonical path you shipped under inst/schemas/.  Validation will load from the installed package.	Conversation answer
Min dims for aggregator/SPCA  Provide `transform_min_dims.myorg.aggregate_runs` and `transform_min_dims.myorg.sparsepca` returning `2L`. Prevents `core_write()` from padding to 3-D when these are the first transforms.  Addendum


⸻

⏩ Next steps for engineers
	1.	Fork the lna repo, scaffold a new transform with

lna::scaffold_transform("myorg.octree_pca")

	2.	Flesh out the two S3 methods using the template.
	3.	Run devtools::test() → ensure your round-trip passes.
	4.	Document the new transform in a vignette snippet.
	5.	Submit a PR or keep it in your own extension package—both coexist happily thanks to the plugin model.

This is a superb set of "belt-and-suspenders" refinements! Each point adds significant practical value for developers of LNA extensions, users, and the long-term health of the LNA ecosystem. The sanity-matrix for package layering is also a very clear and useful way to conceptualize the different components.

Let's integrate these directly into the revised Addendum.

---

## Addendum to LNA R Package Implementation Specification (v1.4)

**12. Extensibility: Add-on Packages & Custom Transforms**

The LNA R package and file format are designed to be extensible, allowing researchers and developers to integrate custom data processing steps (transforms) and analysis modules. This addendum outlines the philosophy and mechanisms for creating and using such extensions.

**12.1. Core Philosophy for Extensions**

*   **Modularity:** New transforms should be self-contained, typically requiring only S3 methods for `forward_step()` and `invert_step()`, a JSON schema, and optional default parameter helpers.
*   **No Core Changes Required:** Well-designed extensions SHOULD NOT necessitate modifications to the core LNA reader/writer machinery, `Plan`, or `DataHandle` classes.
*   **Schema-Driven:** Custom transforms MUST provide a JSON schema (`<type>.schema.json`) defining parameters, enabling validation, default extraction (via `lna:::default_params()`), and documentation.
*   **Adherence to LNA Contracts:** Extensions must respect `DataHandle$stash` data flow and HDF5 layout conventions.
*   **Reproducibility:** Critical parameters (RNG seeds, algorithmic choices, versions of key external tools if output is sensitive) SHOULD be stored in the transform's JSON descriptor. Transforms SHOULD include a `version` field (e.g., `"version": "0.1.0"`) in their descriptor, following semantic versioning principles (major bump for backward-incompatible changes).

**12.2. Implementing a Custom Transform in an Add-on Package**

1.  **Package Structure & Dependencies:**
    *   Logic resides in its own R package (e.g., `myLNAExtensions`).
    *   **S3 Methods:** Implement `forward_step.my_transform_type` and `invert_step.my_transform_type` in `R/`.
    *   **JSON Schema:** Place `<my_transform_type>.schema.json` in `inst/schemas/`. *(Note: `system.file("schemas", package="pkg_name")` is an efficient way to locate these.)*
    *   **Default Parameters (Optional):** Export `lna_default.my_transform_type()`.
    *   **Namespace Collisions:** Add-on packages SHOULD check for name collisions with core LNA transforms in their `.onLoad()` function:
        ```R
        # In .onLoad() of the extension package
        if (exists("my_transform_type", where = "package:lna", inherits = FALSE)) {
          stop("Transform 'my_transform_type' conflicts with a core LNA transform name.")
        }
        ```
    *   **Dependencies:** Declare external R package dependencies (including Rcpp/compiled code) in the add-on's `DESCRIPTION`. This allows the core `lna` package to remain header-only and lightweight for CRAN.

2.  **Transform Naming Convention:**
    *   RECOMMENDED: Reverse domain notation (e.g., `org.my_lab.fmri.custom_pca`).
    *   Core LNA reserves unqualified and simple prefixed names. `lna::check_transform_implementation()` can check collisions.

3.  **HDF5 Path Conventions for Add-on Transforms:**
    *   Store outputs in uniquely named subdirectories (e.g., `/basis/plugins/org.my_lab/custom_pca/matrix`).
    *   All HDF5 paths written or referenced in parameters (e.g., `hrbf_dictionary_descriptor_path`) MUST be absolute HDF5 paths (starting with `/`) to avoid ambiguity when transforms are nested or LNA files are moved.

**12.3. Discovery and Usage of Extensions by LNA Core**

*   **Transform Logic & Schema:** Discovered via S3 dispatch and `lna:::default_params()` (which searches loaded namespaces).
*   **`allow_plugins` Argument & Global Option:**
    *   `read_lna()` argument controls behavior for missing S3 methods ("installed", "prompt", "none").
    *   A global option `options(lna.allow_plugins = "installed")` (or "prompt", "none") can set the default behavior for batch pipelines, mirroring patterns like `ggplot2.discrete.fill`. `read_lna` should respect this option if its own `allow_plugins` argument is not explicitly set.

**12.4. Capabilities Declaration for Extensions**

*   Custom transforms MAY include a `capabilities` object in their JSON descriptor.
    ```json
    "capabilities": {
      "supports_spatial_subsetting": "block-aware" | "full" | "none",
      "supports_temporal_subsetting": true | false,
      "x_my_experimental_flag": true // Experimental flags prefixed with "x_"
    }
    ```
*   The LNA core reader may use standard flags. Experimental flags (prefixed with `x_`) allow extensions to declare features that core LNA can ignore initially but may formalize later.

**12.5. Analysis-Only Transforms**

*   Extensions can define reader-side analysis modules.
*   These consume data from stash or HDF5, produce derived results to stash, and typically do not write primary datasets during a standard read.
*   Their JSON descriptors SHOULD include `"writes_datasets": false` (a new, optional boolean field) to allow LNA core or tooling to optimize pipeline analysis (e.g., skip certain HDF5 write-related safety checks if this transform is at the end of a read-only chain).
*   A core `analysis.identity` transform (which is a no-op, simply passing data through the stash) will be provided. This is useful for bookmarking pipeline checkpoints for debugging, timing, or profiling without altering data.

**12.6. Documentation, Testing, and Large Data in Extensions**

*   **Documentation:** Add-on authors are responsible for thorough documentation. Use `@concept "LNA-extension"` in roxygen2 man pages for discoverability via `help.search()`.
*   **Testing:** Provide comprehensive unit tests. Use `testthat::skip_if_not_installed("lna", minimum_version = "1.4.0")` (or current LNA version) in extension tests to ensure graceful skipping on CI systems if the core LNA version is insufficient.
*   **Large Data/Template Banks (>5MB):** For extension packages shipping large template data (e.g., pre-computed eigenmode banks):
    *   Use `LazyData: true` and `Compression: xz` in the `DESCRIPTION` file.
    *   Provide a helper function (e.g., `download_my_extension_templates()`) that downloads the large data file(s) from a stable repository (e.g., GitHub Releases, Zenodo, Figshare) if not found locally. This keeps the CRAN-submitted package tarball small.

**12.7. Sanity-Matrix for LNA Ecosystem Layering (Illustrative)**

| Layer          | Example Package(s)              | Typical Dependencies               | Purpose                                            |
| :------------- | :------------------------------ | :--------------------------------- | :------------------------------------------------- |
| **Core**       | `lna`                           | base R, hdf5r, digest, jsonlite... | File contract, I/O, basic analytic transforms      |
| **Heavy Maths**| `lna.ext.hkwp`, `lna.ext.hodge`   | igraph, RSpectra, RcppArmadillo... | Sophisticated bases, advanced algorithms           |
| **Data Banks** | `lna.data.mni_laplacians` (≈20MB) | Minimal (downloads data)           | Pre-computed template data (e.g., Laplacian modes) |
| **Viz Tools**  | `lna.viz.viewer`                | shiny, plotly, rgl...              | Optional interactive LNA dashboards & viewers      |
| **Private/Org**| `myOrg.lna.customQC`            | Lab-specific libraries             | Internal QC, proprietary transforms/analyses       |

This layered approach keeps the core `lna` package lean and maintainable, CRAN submissions manageable, and allows users to install only the specialized extensions they need, fostering a rich and flexible ecosystem.

**12.8. Future Considerations for the Extension Ecosystem**

*   (As in previous draft: Central Registry, Inter-Transform Contracts, "LNA Contrib")

This revised addendum provides a robust and forward-looking framework for LNA extensions, encouraging community contribution while maintaining core stability and interoperability.