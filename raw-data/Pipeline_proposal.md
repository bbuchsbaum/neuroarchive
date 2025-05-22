This is an excellent set of "polishing passes" and a superb summary of the DSL's behavior and benefits. The micro-refinements are all valuable and address practical usability and implementation details. The sanity-check matrix and "What do I actually get back?" sections are perfect for user documentation.

Let's integrate these final refinements into the **Definitive Appendix Proposal for the LNA Tidy Pipeline DSL.**

---

## Appendix: LNA Tidy Pipeline DSL - An Ergonomic Façade

This appendix details a high-level, pipeable Domain Specific Language (DSL) designed to sit atop the core `neuroarchive` (LNA) package. This "tidy-pipeline" façade aims to provide a more readable, discoverable, and composable way to define LNA transform sequences, aligning with common R idioms (e.g., `magrittr` pipes), without sacrificing any of the underlying LNA engine's power or flexibility.

### 1. Motivation: Why a Pipeable DSL?

While the core `write_lna()` function is powerful, specifying complex transform chains via character vectors and deeply nested `transform_params` lists can be verbose. A pipeable DSL offers:

*   **Readability:** Transform chains are expressed linearly, left-to-right:
    ```R
    data |>
      as_pipeline() |>
      hrbf(levels = 3, sigma0 = 6.0) |>
      pca(k = 120) |>
      quant(bits = 8) |>
      lna_write("out.lna.h5")
    ```
*   **Discoverability:** Each transform "verb" (e.g., `hrbf()`, `pca()`) has its own help page (`?hrbf`).
*   **Prototyping & Introspection:** Build pipelines incrementally and inspect/validate them before execution.
*   **Composability:** Use standard R control flow and functions to build/modify pipelines.

### 2. Core Design Principles

*   **Thin Syntactic Sugar:** A lightweight layer; the core LNA engine does all heavy lifting.
*   **Deferred Execution:** DSL verbs append specifications to a `pipeline` object; no computation occurs until `lna_write()` is called.
*   **Final Execution:** `lna_write()` translates the pipeline specification into arguments for the core `write_lna()`.
*   **Input Agnostic:** Handles arrays, matrices, `NeuroVec` objects, lists for multi-run data.
*   **Extensibility:** External packages can register their LNA transforms and DSL verbs.

### 3. Key Components

#### 3.1. The `lna_pipeline` R6 Object

The central object storing the pipeline definition.

*   **Fields:**
    *   `input`: The initial data object (or list for multi-run).
    *   `input_summary`: A compact string describing the input data's shape (e.g., "1 run × (240 TR × 65k vox)").
    *   `runs`: (Optional) Character vector of run identifiers.
    *   `steps`: A list of transform specifications: `list(type = "transform_type", params = list(...), descriptor_basename = "NN_transform_type")`.
    *   `engine_opts`: (Optional) List of hints for `core_write`, e.g., `list(chunk_mb_suggestion = 256)`.
*   **Key Methods:**
    *   `initialize()`: Creates an empty pipeline.
    *   `set_input(x, run_ids = NULL, chunk_mb = NULL)`: Sets input, derives `input_summary`, handles `runs`, stores `chunk_mb` in `engine_opts`.
    *   `add_step(step_spec)`: Appends a transform step.
    *   `get_last_step_spec()`: Returns the last added step's specification.
    *   **Introspection & Modification:**
        *   `print()`: User-friendly display including `input_summary` and steps with resolved parameters (user-supplied parameters highlighted, defaults styled subtly via `pillar::style_subtle()`).
        *   `steps()`: Returns the raw list of step specifications.
        *   `get_step(index_or_type)`: Retrieves a step by 1-based index or the *last* step matching a `type` string.
        *   `modify_step(index_or_type, new_params_list)`: Modifies parameters for a step, re-running parameter merging.
        *   `insert_step(step_spec, after_index_or_type)`: Inserts a new step.
        *   `remove_step(index_or_type)`: Removes a step.
        *   `validate_params(strict = FALSE)`: Dry-runs schema validation for all step parameters. If `strict=FALSE` (default), collects all failures and warns. If `strict=TRUE`, aborts on the first violation.
        *   `diagram(engine = c("grViz", "ascii", "dot"))`: Generates a visual pipeline graph. For ASCII, parameter text is clipped (e.g., ±30 chars with ellipsis) for console readability.
            *   A low-level helper `infer_embed_step(prev_step_spec)` can be exposed to aid authors of context-aware verbs, returning `list(suggested_step_type, auto_filled_params)`.

#### 3.2. Initiating a Pipeline: `as_pipeline()`

```R
as_pipeline <- function(x, run_ids = NULL, chunk_mb_suggestion = NULL) {
  pipe <- lna_pipeline$new()
  # Validate x (is.list, is.array, inherits("NeuroVec"), etc.)
  # ...
  pipe$set_input(x, run_ids, chunk_mb_suggestion) # sets pipe$input, pipe$runs, pipe$input_summary, pipe$engine_opts
  pipe
}
```
*   **Multi-run Input:** If `x` is a list, `run_ids` defaults to `names(x)` (or `run-01`, `run-02`... if unnamed). `pipe$input` stores the list, `pipe$runs` stores the IDs. Core `write_lna()` handles the multi-run logic.

#### 3.3. DSL Verbs (e.g., `hrbf()`, `pca()`)

Each verb is an R function.
*   **Signature Pattern:** `verb_name(data_or_pipe, param1 = default1, ...)`
    *   If `data_or_pipe` is not an `lna_pipeline`, it's treated as input data, and `as_pipeline(data_or_pipe)` is called internally.
*   **Functionality:**
    1.  Ensures it has an `lna_pipeline` object.
    2.  **Parameter Merging (within verb):**
        a.  `pars <- lna:::default_params("lna_transform_type_string")` (memoised via `memoise::memoise()` for performance).
        b.  `pars <- utils::modifyList(pars, lna_options("lna_transform_type_string"))`.
        c.  `pars <- utils::modifyList(pars, list(...))` (user-supplied args).
    3.  Constructs `step_spec = list(type = "lna_transform_type_string", params = pars)`.
    4.  Calls `pipeline_obj$add_step(step_spec)`.
    5.  Returns the `pipeline_obj`.

#### 3.4. Terminal Verb: `lna_write()`

```R
lna_write <- function(pipeline_obj, file, ..., 
                      .verbose = TRUE, .checksum = "sha256") {
  # ... (validation of pipeline_obj) ...
  transform_types <- vapply(pipeline_obj$steps, function(s) s$type, character(1))
  transform_params_list <- lapply(pipeline_obj$steps, function(s) s$params)
  names(transform_params_list) <- transform_types

  # Pass engine_opts if present
  core_write_args_extra <- pipeline_obj$engine_opts %||% list()
  # ... (merge with ..., ensuring .verbose, .checksum take precedence if also in ...)

  tryCatch({
    # Core call
    write_lna(
      x = pipeline_obj$input, file = file,
      run_ids = pipeline_obj$runs,
      transforms = transform_types,
      transform_params = transform_params_list,
      # Pass other args like header, mask, plugins from ...
      # Pass .verbose, .checksum, and content of core_write_args_extra
      ... 
    )
  }, lna_error = function(e) {
    step_idx_core <- attr(e, "step_index", exact = TRUE) # 0-based from core
    ttype_core <- attr(e, "transform_type", exact = TRUE)
    
    bullet <- "Pipeline execution failed."
    if (!is.null(step_idx_core) && !is.null(ttype_core)) {
      # Attempt to map core 0-based index to 1-based DSL step
      # Find first step in pipeline$steps that has matching type and whose position corresponds to step_idx_core
      # This requires careful handling if types are repeated. A simpler approach:
      # dsl_step_num <- step_idx_core + 1 (if types were unique and directly mapped)
      # For robustness, might need to search pipeline_obj$steps
      dsl_step_info <- pipeline_obj$steps[[step_idx_core + 1]] # Assuming simple mapping for now
      if (!is.null(dsl_step_info) && dsl_step_info$type == ttype_core) {
         bullet <- sprintf("Pipeline failure in step %d (type='%s')",
                              step_idx_core + 1, ttype_core)
      } else { # Fallback if direct mapping fails (e.g. error before transform loop)
         bullet <- sprintf("Pipeline failure near transform '%s' (internal index %d)",
                              ttype_core %||% "unknown", step_idx_core %||% "unknown")
      }
    }
    rlang::abort(message = c(bullet, "i" = conditionMessage(e)), parent = e, .subclass = class(e))
  }, error = function(e) { # Generic errors
    rlang::abort(message = c("Pipeline execution failed.", "i" = conditionMessage(e)), parent = e)
  })
}
```

#### 3.5. Core Verbs (Shipped with `lna` package)

Includes `as_pipeline()`, `hrbf()`, `pca()`, `embed()` (context-aware), `quant()`, `delta()`, `temporal()`, and `lna_write()`.

#### 3.6. Extensibility: Registering Verbs & Templates

*   **`register_lna_verb(verb_name_symbol_or_string, lna_transform_type_string, force = FALSE)`:**
    *   Default slugging: `my.org.filter` -> `my_org_filter`. Warns on collision unless `force = TRUE` (silently replaces, useful for package reloading during development).
*   **`register_lna_template(template_name_string, template_function)`:** As before.
*   **`apply_template(pipeline_obj, template_name_string, ...)`:** As before. Supports parameter overrides via `transform_type.param_name = value` or `list(transform_type = list(param_name = value))`. Internally converts dotted names.

### 4. What is Returned?

*   **Pipeline Construction Phase (e.g., `data |> as_pipeline() |> hrbf()`):** Returns the modified `lna_pipeline` R6 object (a lightweight plan).
*   **Materialization Phase (`... |> lna_write("file.lna.h5")`):**
    *   The primary output is the `.lna.h5` archive on disk.
    *   `lna_write()` invisibly returns the result from the core `write_lna()` (a list containing `file` path, `plan` object, `header`).
    *   Intermediate large data objects are streamed to disk and not kept in R's workspace.

### 5. Command-Line Interface (CLI) Parity (Future Consideration)

A small wrapper script could enable execution of serialized `lna_pipeline` objects (e.g., saved as `.Rds` files) from the command line: `Rscript -e "lna::run_pipeline_from_file('pipeline.Rds', 'output.lna.h5')"`, facilitating integration with workflow managers like Snakemake or CWL.

### Conclusion

The LNA Tidy Pipeline DSL offers a significant enhancement in usability and readability for defining LNA transform sequences. It maintains full compatibility with the powerful core LNA engine, ensuring robustness and reproducibility, while providing an intuitive, R-idiomatic interface for users. Its extensibility allows the LNA ecosystem to grow with community-contributed transforms and standardized processing templates.