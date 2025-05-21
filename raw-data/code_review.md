Okay, this is a substantial codebase, and a static audit like this is a great way to catch potential issues. I'll go through the provided `repomix-output.xml` structure and highlight areas of concern, potential bugs, reasons tests might fail, and then generate a ticketed list for fixes.

**Overall Impression:**
The codebase shows a good adherence to the LNA v1.4 specification, with thoughtful implementation of core classes (`Plan`, `DataHandle`, `lna_reader`), HDF5 utilities, and transform dispatch. Error handling (`abort_lna`), parameter resolution (`resolve_transform_params`), and schema usage are generally well-structured. The addition of helper functions like `run_transform_step` and `handle_missing_methods` improves robustness. The test suite appears to be growing to cover many components.

## Code Audit: Analysis & Potential Issues

Here's a breakdown by common themes and specific files:

**1. HDF5 Resource Management & Paths:**

*   **`R/utils_json.R` (`write_json_descriptor`):** The `on.exit` block correctly closes HDF5 objects. This is good.
*   **`R/materialise.R` (`materialise_plan`):**
    *   Group creation checks (`if (h5$exists(...) ... if (!inherits(obj, "H5Group")) ... obj$close()`) are good for robustness against malformed existing files.
    *   **Potential Issue (Minor):** In `write_payload`'s `attempt` function, `h5_write_dataset` returns a dataset object `ds`. If `h5_write_dataset` itself ensures the `ds` it returns is closed (or if its returned value isn't actually an open dataset object but just a signal of success/failure after closing), then it's fine. The current `h5_write_dataset` seems to return `dset` which it closes internally. However, if `h5_write_dataset` *were* to return an open `H5D` object, `attempt` should close it. The current `h5_write_dataset` *does* close the dataset it creates (`invisible(dset)` is after any explicit close in its own scope, but the last line `invisible(dset)` will return the `dset` object which `hdf5r` might auto-close if not assigned, or keep open if assigned. Best to be explicit: `dset$close(); invisible(TRUE)`). The current `materialise_plan` receives `NULL` from `attempt` due to the structure `ds <- ...; if(inherits(ds, "H5D")) ds$close(); NULL`. This seems okay.
    *   **Potential Bug/Refinement:** The `dtype_size` guessing logic inside `materialise_plan -> write_payload` for chunking heuristics:
        ```R
        dtype <- guess_h5_type(data)
        dtype_size <- dtype$get_size(variable_as_inf = FALSE)
        if (!is.finite(dtype_size) || dtype_size <= 0) {
          dtype_size <- 1L # Default to 1 if invalid
        }
        # This part overrides the above, which might be fine but looks a bit redundant:
        if (is.integer(data)) {
          dtype_size <- 4L
        } else if (is.double(data)) {
          dtype_size <- 8L
        }
        if (inherits(dtype, "H5T")) dtype$close()
        ```
        The explicit `is.integer/is.double` check is fine, but if `guess_h5_type` is robust, the initial `dtype$get_size` should be sufficient. Ensure `guess_h5_type` correctly maps R types to HDF5 types whose sizes are then used. The current `guess_h5_type` returns an `H5T` object, so `dtype$get_size()` *should* work. The `if (is.integer...)` might be a safeguard.
    *   Plugin name validation in `materialise_plan` (`if (grepl("/", nm))`) is good to prevent path manipulation.
*   **`R/utils_hdf5.R`:**
    *   `h5_write_dataset`: Creates parent groups. Uses `guess_chunk_dims`. Handles compression. Closes the dataset it creates via `invisible(dset)` if `dset` is the last expression and not assigned, or if `dset$close()` is called before `invisible(dset)`. Currently, it seems the `dset` object is returned and relies on R's GC or `hdf5r`'s auto-closing. It's safer for `h5_write_dataset` to explicitly close `dset` before returning (e.g., `dset$close(); invisible(TRUE)` or similar signal of success).
    *   `h5_read`, `h5_read_subset`: Correctly use `tryCatch` and `finally` to close dataset handles.
    *   `sanitize_run_id`: Regex `^run-[0-9]{2}$` is good.
*   **`R/reader.R` (`lna_reader$initialize`):** The `on.exit(close_h5_safely(h5))` *before* `self$h5 <- h5` is crucial. If `discover_run_ids` or `resolve_run_ids` fails *after* `open_h5` but *before* `self$h5` is assigned, the file would be left open without this. This is good.
*   **`R/validate.R` (`validate_lna`):** The nested `on.exit(if (inherits(dset, "H5D")) dset$close(), add = TRUE)` and `on.exit(if (inherits(dt, "H5T")) dt$close(), add = TRUE)` inside the loop for dataset validation are good for resource cleanup.

**2. Parameter Handling & Defaults:**

*   **`R/utils_defaults.R`:**
    *   `resolve_transform_params`: Merge order (schema -> package opts -> user) seems correct. `utils::modifyList(..., keep.null = TRUE)` is appropriate for deep merging while preserving explicit NULLs.
    *   `.extract_schema_defaults`: Handles `properties` and `items` (for arrays). This is good for complex schemas.
    *   `default_params` / `required_params`: Schema finding via `system.file` across `loadedNamespaces()` is a reasonable approach.
*   **`R/options.R`:** Default options list (`write.compression_level`, `write.chunk_target_mib`, `quant`, `delta`) is a good start.
*   **Transforms:** Most transforms correctly use `p <- desc$params %||% list()` and then `%||%` for individual parameters, which is robust.

**3. Control Flow & Logic:**

*   **`R/core_write.R`:**
    *   Input validation (3D check, mask, header, plugins) added later improves robustness.
    *   `run_transform_step` for consistent error provenance is excellent.
*   **`R/core_read.R`:**
    *   `allow_plugins` logic (including interactive fallback and `handle_missing_methods`) seems robust.
    *   Subset parameter validation (e.g., `roi_mask` logical, `time_idx` numeric) is important.
    *   `runtime_validate_step` inside the loop (added later) is a good proactive check.
*   **`R/reader.R` (`lna_reader$data`):**
    *   Check for closed reader (`if (is.null(self$h5) || !self$h5$is_valid())`) is crucial.
    *   Cache logic using `identical(params, self$cache_params)` is correct.
    *   Replicates much of `core_read`'s inverse loop logic, which is necessary.
*   **`R/transform_quant.R`:**
    *   Parameter validation at the start of `forward_step.quant` (bits, method, center, scope) is good.
    *   Fallback for `scale_scope="voxel"` on non-4D data is a sensible default.
    *   Handling of `scale == 0` in `.quantize_global` and `.quantize_voxel` (setting `scale=1` and `q=0L`) prevents division by zero and ensures a defined output.
*   **`R/transform_temporal.R`:**
    *   `temporal_basis` S3 generic with methods for "dct", "bspline", "dpss", "polynomial", "wavelet" is a flexible design.
    *   Wavelet basis requiring power-of-two length is correctly checked.
*   **`R/discover.R` (`discover_transforms`):**
    *   Regex `^(\\d+)_([^_]+)\\.json$` seems correct for `NN_type.json`.
    *   Error handling for invalid names, non-numeric index, non-contiguous sequence, and sequence not starting at 0 is comprehensive and uses specific error classes.
*   **`R/plan.R` (`Plan$get_next_filename`):** Validation for invalid characters in `type` (`grepl("\.\.", type)` etc.) is a good security/robustness measure.

**4. Error Handling & Validation:**

*   **`R/utils_error.R` (`abort_lna`):** Consistent use of `rlang::abort` with custom subclasses and `location`/`parent` is best practice.
*   **`R/utils_transform.R` (`handle_missing_methods`, `run_transform_step`):** These centralize logic for dealing with missing S3 methods and add error provenance, improving robustness and debuggability.
*   **`R/validate.R`:** `validate_lna` is becoming quite comprehensive. `runtime_validate_step` is a good addition for `core_read`.
*   **Test: `tests/testthat/test-core_read.R`:**
    *   `core_read validate=TRUE checks dataset existence`: This tests `runtime_validate_step`. If `desc$datasets[[1]]$path` is `/scans/run-01/missing`, `assert_h5_path` (called by `runtime_validate_step`) should throw `lna_error_missing_path`.
    *   `core_read validate=TRUE checks required params`: This also tests `runtime_validate_step`. For `embed`, `basis_path` is required. If `desc$params` is empty, `required_params("embed")` should identify `basis_path` as missing.
*   **Test: `tests/testthat/test-discover.R`:** Errors on invalid names use `stop()`, while spec-defined sequence errors use `abort_lna(..., .subclass = "lna_error_sequence")`. The tests for invalid names check for `stop()` (e.g., `expect_error(..., "Invalid object name")`). This is slightly inconsistent but might be intentional (internal validation vs. spec-defined file format error). The latest version of `discover.R` uses `abort_lna` for invalid names too (`lna_error_descriptor`), so tests should be updated.

**5. R6 Class Specifics:**

*   **`R/handle.R` (`DataHandle`):** The `$with` method uses `self$clone(deep = TRUE)`. The logic for finding `allowed_fields` using `self$.__enclos_env__$private$.class` or `self$.__enclos_env__$self` is a bit unusual but might be a workaround for some R6 internal behavior or specific environment where `class(self)[1]` might not behave as expected. Standard R6 usually allows direct access to `self$public_fields` or similar introspection if needed, but this works.
*   **`R/reader.R` (`lna_reader`):** `finalize` method calling `$close()` is good for GC-based cleanup. `$close()` also clearing caches is good.

**6. Specific Potential Bugs / Test Issues:**

*   **`R/discover.R` & `tests/testthat/test-discover.R`:**
    *   The `discover_transforms` function now uses `abort_lna(..., .subclass = "lna_error_descriptor")` for invalid names. The tests (`test-discover.R`) for invalid names (e.g., "invalid_name.txt", "aa_mask.json") expect a simple `stop()` message like `"Invalid object name..."`. These tests will fail because they'll get an `lna_error_descriptor` condition object instead of matching the string.
    *   **Fix:** Update `test-discover.R` to use `expect_error(..., class = "lna_error_descriptor")`.

*   **`R/utils_json.R` & `tests/testthat/test-utils_json.R`:**
    *   `read_json_descriptor`: If `assert_h5_path` fails, it throws `lna_error_missing_path`. The test `read_json_descriptor error handling works -> Read non-existent descriptor` expects a specific string `"JSON descriptor dataset 'missing_desc' not found"`. This test will likely fail as `assert_h5_path` will error out first with its own message and class.
    *   **Fix:** The `read_json_descriptor` error test for missing descriptor should expect `class = "lna_error_missing_path"`.
    *   `read_json_descriptor`: For invalid JSON and non-string dataset, it now uses `abort_lna` with specific subclasses (`lna_error_json_parse`, `lna_error_invalid_descriptor`). The tests need to be updated to expect these classes.

*   **`R/transform_quant.R` (`forward_step.quant`):**
    *   The parameter validation for `scope`:
        ```R
        if (!(is.character(scope) && length(scope) == 1 &&
              scope %in% c("global", "voxel"))) {
          abort_lna(
            sprintf("Invalid scale_scope '%s'", scope),
            .subclass = "lna_error_validation",
            location = "forward_step.quant:scale_scope"
          )
        # This block is now dead code due to the abort_lna above it.
        if (!scope %in% c("global", "voxel")) {
          warning(sprintf("unknown scale_scope '%s'; falling back to 'global'", scope))
          scope <- "global"
        }
        ```
    *   **Fix:** The second `if (!scope %in% ...)` block is unreachable because the `abort_lna` will have already exited if `scope` is invalid. If the intention was to allow other `scope` values with a warning and fallback, the `abort_lna` is too strict. If only "global" and "voxel" are allowed, the warning block is correctly dead. Assuming strict validation, the warning block can be removed. The spec doesn't mention a fallback for `scale_scope`.

*   **`R/transform_temporal.R` (`forward_step.temporal`):**
    *   `sanitize_run_id(run_id)` is called. If `handle$current_run_id` is NULL (e.g., for a global temporal basis not tied to a specific run scan), `run_id` becomes `"run-01"`. If this default "run-01" is then sanitized and `sanitize_run_id` *strictly* requires the "run-XX" format, this could be an issue if `handle$current_run_id` was legitimately `NULL` for a global resource. However, `sanitize_run_id` seems to be for path construction, and `/temporal/...` paths don't use `run_id`. `/scans/...` paths *do* use `run_id`. This is likely fine as `sanitize_run_id` is defensive.
    *   The test `test-transform_temporal.R` -> `temporal transform rejects unsupported kind` directly calls `core_write`. If an unsupported `kind` (like "dpss" before it was implemented) is passed, `temporal_basis.default` will be called, which correctly uses `abort_lna`. The test should catch this.

*   **`tests/testthat/test-transform_embed_inverse.R`:**
    *   The first test `invert_step.embed reconstructs dense data` has a duplicated block of code after the loop.
        ```R
        # ... inside for loop ...
        h5$close_all()
      } # end for loop
      # This block is outside the loop and h5 is closed
      handle <- DataHandle$new(initial_stash = list(coef = coef_mat), h5 = h5) # h5 is closed here
      h <- invert_step.basis("basis", desc, handle) # This will fail
      # ...
    ```
    *   The lines after the `for` loop will try to use a closed `h5` handle and call `invert_step.basis` instead of `invert_step.embed`. This test block needs to be fixed or removed. It seems like a copy-paste error from `test-transform_basis_inverse.R`.
    *   **Fix:** Remove or correct the duplicated block after the loop.
    *   The test `invert_step.embed errors when datasets are missing` has a `}` closing the `with_mocked_bindings` too early, and the second `expect_error` block is outside. This looks like a syntax error in the test.
    *   **Fix:** Correct the curly brace placement in the "datasets missing" test.

*   **`tests/testthat/test-core_write.R` (`mask is validated and stored`):**
    *   `expect_true(all(res$handle$meta$mask == msk))` compares the mask array. For large masks, this could be slow or print a lot on failure. `expect_identical` or checking properties (dims, sum) might be better if exact value comparison isn't strictly needed for this particular check point (though `core_write` *does* store it in meta). This is minor.

*   **`R/reader.R` (`lna_reader$data`):**
    *   The logic `if (nrow(transforms) > 0) { missing_methods <- ...; skip_types <- handle_missing_methods(...) }` is duplicated from `core_read.R`. This isn't a bug, but could be refactored into a shared helper if desired for DRYness, though keeping it separate might be fine given `lna_reader`'s specific caching and state.

*   **`R/utils_hdf5.R` (`guess_chunk_dims`):**
    *   Uses `lna_options("write.chunk_target_mib")[[1]]`. If this option isn't set or is NULL, `target_mib` will be NULL. The code then tries `as.numeric(target_mib) * 1024^2`. `as.numeric(NULL)` is `numeric(0)`. `numeric(0) * 1024^2` is `numeric(0)`. This might lead to `hdf5r::guess_chunks` using its own default if `target` is not a positive scalar.
    *   **Fix:** The `lna_options("write.chunk_target_mib")[[1]]` should have a fallback if NULL, e.g., `target_mib <- lna_options("write.chunk_target_mib")[[1]] %||% 1`. (The `R/options.R` has been updated to provide a default for this, so this is likely fine now).

**7. Documentation & Comments:**
*   Generally good use of roxygen2.
*   Some internal functions could benefit from slightly more detailed comments explaining non-obvious logic, but most are clear from context and the spec.

## Ticketed List of Fixes & Refinements

This list prioritizes correctness and testability.

**Epic FX-CORE: Core Logic & Utility Fixes**

| #       | Ticket                                             | Description / Deliverables                                                                                                                                                                                                                                                                                               | Files Affected                          | Severity |
| :------ | :------------------------------------------------- | :----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :-------------------------------------- | :------- |
| FX-C-1  | **Test `discover_transforms` for Error Classes**   | Update tests in `test-discover.R` that check for errors on invalid transform names. They currently expect `stop()` messages. They should expect `abort_lna(..., class = "lna_error_descriptor")` as per the updated `discover_transforms` implementation.                                                                | `tests/testthat/test-discover.R`        | Medium   |
| FX-C-2  | **Test `read_json_descriptor` for Error Classes**  | Update error checks in `test-utils_json.R` for `read_json_descriptor`: <br> 1. For non-existent descriptor, expect `class = "lna_error_missing_path"`. <br> 2. For invalid JSON, expect `class = "lna_error_json_parse"`. <br> 3. For non-string dataset, expect `class = "lna_error_invalid_descriptor"`.                       | `tests/testthat/test-utils_json.R`      | Medium   |
| FX-C-3  | **Refine `h5_write_dataset` Return/Close**         | Ensure `h5_write_dataset` in `R/utils_hdf5.R` explicitly closes the `H5D` object it creates before returning, or change its return to a simple success indicator (e.g., `invisible(TRUE)`). This makes its resource management more explicit and less reliant on GC or `hdf5r` auto-closing.                         | `R/utils_hdf5.R`                        | Low      |
| FX-C-4  | **Clarify `dtype_size` in `materialise_plan`**     | Review `dtype_size` calculation within `materialise_plan` -> `write_payload`. Ensure it robustly uses `guess_h5_type(data)$get_size()` and simplify the subsequent `is.integer`/`is.double` overrides if they are truly redundant and covered by `guess_h5_type`. Ensure `dtype$close()` is called.           | `R/materialise.R`                       | Low      |
| FX-C-5  | **Correct `quant` param validation logic**         | In `R/transform_quant.R` (`forward_step.quant`), remove the dead/unreachable `if (!scope %in% c("global", "voxel"))` warning block, as the preceding `abort_lna` for invalid `scale_scope` makes it unreachable. Confirm strict validation (erroring) is the intended behavior for `scale_scope`.                 | `R/transform_quant.R`                   | Low      |
| FX-C-6  | **Fix `test-transform_embed_inverse.R` issues**    | 1. Remove/correct the duplicated code block after the `for` loop in the "reconstructs dense data" test, which attempts to use a closed `h5` handle. <br> 2. Correct curly brace placement in the "errors when datasets are missing" test. Ensure the second `expect_error` is correctly scoped within `with_mocked_bindings`. | `tests/testthat/test-transform_embed_inverse.R` | High (Test Failure) |

**Epic FX-TEST: Test Suite Enhancements & Corrections**

| #       | Ticket                                         | Description / Deliverables                                                                                                                                                                                                                                                                     | Files Affected                                                              | Severity |
| :------ | :--------------------------------------------- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :-------------------------------------------------------------------------- | :------- |
| FX-T-1  | **Align `discover_transforms` error tests**    | (Covered by FX-C-1) Update tests for `discover_transforms` to expect `lna_error_descriptor` for invalid names, not generic `stop()` messages.                                                                                                                                                        | `tests/testthat/test-discover.R`                                            | Medium   |
| FX-T-2  | **Align `read_json_descriptor` error tests**   | (Covered by FX-C-2) Update tests for `read_json_descriptor` to expect specific `lna_error_*` classes for missing files, parse errors, and invalid descriptor types.                                                                                                                                     | `tests/testthat/test-utils_json.R`                                          | Medium   |
| FX-T-3  | **Robustness test for `lna_reader` initialize**| Add a test to specifically verify that if `lna_reader$initialize` fails (e.g., due to `discover_run_ids` erroring because no runs found), the HDF5 file handle opened by `open_h5` is indeed closed by the `on.exit` handler.                                                                      | `tests/testthat/test-reader.R`                                              | Medium   |
| FX-T-4  | **Review `test-transform_embed_inverse.R` syntax** | (Covered by FX-C-6) The test `invert_step.embed errors when datasets are missing` has a syntax error with curly braces that needs fixing.                                                                                                                                                         | `tests/testthat/test-transform_embed_inverse.R`                             | High (Test Failure) |

This audit and ticket list should help in identifying and rectifying potential issues before extensive testing. The codebase is in good shape, and these are mostly refinements or test alignment issues.