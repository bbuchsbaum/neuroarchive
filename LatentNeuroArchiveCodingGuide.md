Okay, this final feedback targets crucial details for implementation robustness and clarity. Here is the v1.4 specification for the `lna` R package, incorporating these last refinements.

---

**LNA R Package Implementation Specification (v1.4 - Implementation Ready)**

**(Summary of Changes from v1.3):** This definitive version incorporates final clarifications based on detailed review. Key refinements include: specifying the exact parameter default merging order, clarifying checksum computation timing, defining non-interactive behavior for plugin prompting, specifying error classes, noting fork-safety issues with schema caches, adding warnings for namespace collisions, enhancing HDF5 robustness retry logic for chunk sizes, detailing lazy reader closure idempotency, improving scaffolded code, and clarifying file opening modes for parallel writing. The specification is now considered complete and ready for implementation.

**1. Overview & Purpose**

(No changes from v1.0 - goals remain the same)

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
  input_key <- if (handle$exists("aggregated_matrix")) "aggregated_matrix" else "dense_mat"
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

This combined recipe illustrates how the modular transform system allows complex workflows, like cross-run aggregation followed by encoding, to be built by composing independent, well-defined transform steps. The key is careful definition of the `inputs` and `outputs` contract for each step.