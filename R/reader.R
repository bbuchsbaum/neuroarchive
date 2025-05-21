#' lna_reader Class for Lazy Reading
#'
#' @description Provides deferred data loading from LNA files. The
#'   reader keeps an HDF5 file handle open and runs the inverse
#'   transform pipeline on demand via `$data()`. Subsetting parameters
#'   can be stored with `$subset()`.
#'
#' @details
#' Create an instance via `read_lna(file, lazy = TRUE)` or directly
#' using `lna_reader$new()`.  Call `$subset()` to store ROI or time
#' indices and `$data()` to materialise the data.  Always call
#' `$close()` when finished.
#'
#' @examples
#' r <- read_lna("example.lna.h5", lazy = TRUE)
#' r$subset(time_idx = 1:10)
#' dat <- r$data()
#' r$close()
#'
#' @keywords internal
lna_reader <- R6::R6Class("lna_reader",
  public = list(
    #' @field file Path to the underlying LNA file
    file = NULL,
    #' @field h5 Open H5File handle
    h5 = NULL,
    #' @field core_args List of arguments forwarded to `core_read`
    core_args = NULL,
    #' @field run_ids Selected run identifiers
    run_ids = NULL,
    #' @field allow_plugins Stored allow_plugins behaviour from `read_lna`
    allow_plugins = "installed",
    #' @field current_run_id Run identifier currently used
    current_run_id = NULL,
    #' @field subset_params Stored subsetting parameters
    subset_params = NULL,
    #' @field data_cache Cached DataHandle from `$data()`
    data_cache = NULL,
    #' @field cache_params Parameters used for `data_cache`
    cache_params = NULL,

    #' @description
    #' Create a new `lna_reader`
    #' @param file Path to an LNA file
    #' @param core_read_args Named list of arguments for `core_read`
    initialize = function(file, core_read_args) {
      stopifnot(is.character(file), length(file) == 1)
      self$file <- file
      self$core_args <- core_read_args
      self$allow_plugins <- core_read_args$allow_plugins %||% "installed"
      subset_params <- list()
      if (!is.null(core_read_args$roi_mask)) {
        roi <- core_read_args$roi_mask
        if (inherits(roi, "LogicalNeuroVol")) roi <- as.array(roi)
        subset_params$roi_mask <- roi
      }
      if (!is.null(core_read_args$time_idx)) {
        subset_params$time_idx <- as.integer(core_read_args$time_idx)
      }
      self$subset_params <- subset_params

      h5 <- open_h5(file, mode = "r")
        on.exit(neuroarchive:::close_h5_safely(h5))

      runs_avail <- discover_run_ids(h5)
      runs <- resolve_run_ids(core_read_args$run_id, runs_avail)
      if (length(runs) == 0) {
        abort_lna("run_id did not match any runs", .subclass = "lna_error_run_id")
      }
      if (length(runs) > 1) {
        warning("Multiple runs matched; using first match for lazy reader")
        runs <- runs[1]
      }

      self$run_ids <- runs
      self$current_run_id <- runs[1]
      self$h5 <- h5
      on.exit(NULL, add = FALSE)
    },

    #' @description
    #' Close the HDF5 handle. Safe to call multiple times.
    close = function() {
      if (!is.null(self$h5)) {
          neuroarchive:::close_h5_safely(self$h5)
        self$h5 <- NULL
      }
      self$data_cache <- NULL
      self$cache_params <- NULL
      invisible(NULL)
    },

    #' @description
    #' Print summary of the reader
    print = function(...) {
      status <- if (!is.null(self$h5) && self$h5$is_valid()) "open" else "closed"
      cat("<lna_reader>", self$file, "[", status, "] runs:", paste(self$run_ids, collapse = ","), "\n")
      invisible(self)
    },

    #' @description
    #' Store subsetting parameters for later `$data()` calls.
    #' Only `roi_mask` and `time_idx` are accepted.
    #' @param ... Named parameters such as `roi_mask`, `time_idx`
    subset = function(...) {
      allowed <- c("roi_mask", "time_idx")
      args <- list(...)
      if (length(args) > 0) {
        if (is.null(names(args)) || any(names(args) == "")) {
          abort_lna(
            "subset parameters must be named",
            .subclass = "lna_error_validation",
            location = "lna_reader:subset"
          )
        }
        unknown <- setdiff(names(args), allowed)
        if (length(unknown) > 0) {
          abort_lna(
            paste0("Unknown subset parameter(s): ",
                   paste(unknown, collapse = ", ")),
            .subclass = "lna_error_validation",
            location = "lna_reader:subset"
          )
        }
        self$subset_params <- utils::modifyList(
          self$subset_params, args, keep.null = TRUE
        )
      }
      invisible(self)
    },

    #' @description
    #' Load and reconstruct data applying current subsetting.
    #' @param ... Optional subsetting parameters overriding stored ones
    #' @return A `DataHandle` object representing the loaded data
      data = function(...) {
        args <- list(...)
        if (is.null(self$h5) || !self$h5$is_valid()) {
          abort_lna(
            "lna_reader is closed",
            .subclass = "lna_error_closed_reader",
            location = sprintf("lna_reader:data:%s", self$file)
          )
        }

        params <- self$subset_params
        if (length(args) > 0) {
          params <- utils::modifyList(params, args, keep.null = TRUE)
        }
      if (!is.null(self$data_cache) && identical(params, self$cache_params)) {
        return(self$data_cache)
      }

      h5 <- self$h5
      handle <- DataHandle$new(
        h5 = h5,
        subset = params,
        run_ids = self$run_ids,
        current_run_id = self$current_run_id
      )
      tf_group <- h5[["transforms"]]
      transforms <- discover_transforms(tf_group)
      allow_plugins <- self$allow_plugins
      if (identical(allow_plugins, "prompt") && !rlang::is_interactive()) {
        allow_plugins <- "installed"
      }
      if (nrow(transforms) > 0) {
        missing_methods <- transforms$type[
          vapply(
            transforms$type,
            function(t) is.null(getS3method("invert_step", t, optional = TRUE)),
            logical(1)
          )
        ]
        skip_types <- handle_missing_methods(
          missing_methods,
          allow_plugins,
          location = sprintf("lna_reader:data:%s", self$file)
        )
        if (length(skip_types) > 0) {
          transforms <- transforms[!transforms$type %in% skip_types, , drop = FALSE]
        }
      }
      if (nrow(transforms) > 0) {
        for (i in rev(seq_len(nrow(transforms)))) {
          name <- transforms$name[[i]]
          type <- transforms$type[[i]]
          step_idx <- transforms$index[[i]]
          desc <- read_json_descriptor(tf_group, name)
          handle <- run_transform_step("invert", type, desc, handle, step_idx)
        }
      }

      output_dtype <- self$core_args$output_dtype
      if (identical(output_dtype, "float16") && !has_float16_support()) {
        abort_lna(
          "float16 output not supported",
          .subclass = "lna_error_float16_unsupported",
          location = sprintf("lna_reader:data:%s", self$file)
        )
      }
      handle$meta$output_dtype <- output_dtype
      handle$meta$allow_plugins <- allow_plugins

      self$data_cache <- handle
      self$cache_params <- params
      handle
    }
  ),
  
  private = list(
    #' @description
    #' Finalizer called by GC
    finalize = function() {
      self$close()
    }
  )
)
