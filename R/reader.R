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
#' tmp <- tempfile(fileext = ".h5")
#' arr <- array(rnorm(4 * 4 * 1 * 10), dim = c(4, 4, 1, 10))
#' write_lna(arr, tmp, transforms = "quant")
#' r <- read_lna(tmp, lazy = TRUE)
#' r$subset(time_idx = 1:10)
#' dat <- r$data()
#' r$close()
#' unlink(tmp)
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
    #' @field qpp_tile_cache Reader-owned bounded cache of decoded QPP tiles.
    qpp_tile_cache = NULL,

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
      self$qpp_tile_cache <- .qpp_new_tile_cache()

      h5 <- open_h5(file, mode = "r")
      on.exit(close_h5_safely(h5))

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
      self$data_cache <- NULL
      self$cache_params <- NULL
      .qpp_tile_cache_clear(self$qpp_tile_cache)
      self$qpp_tile_cache <- NULL
      if (!is.null(self$h5)) {
        # A wildcard lazy read owns one H5File identifier per run. `close_all()`
        # closes every identifier for the same physical file, so use the
        # reader-owned identifier close here to keep sibling readers alive.
        if (isTRUE(self$h5$is_valid)) {
          tryCatch(
            self$h5$close(),
            error = function(e) {
              warning(paste("Error closing lazy LNA handle:", conditionMessage(e)))
            }
          )
        }
        self$h5 <- NULL
      }
      invisible(NULL)
    },

    #' @description
    #' Print summary of the reader
    print = function(...) {
      status <- if (!is.null(self$h5) && self$h5$is_valid) "open" else "closed"
      cat("<lna_reader>", self$file, "[", status, "] runs:", paste(self$run_ids, collapse = ","), "\n")
      invisible(self)
    },

    #' @description
    #' Store subsetting parameters for later `$data()` calls.
    #' `parcel_idx` is also accepted for parcel-space consumers.
    #' @param ... Named parameters such as `roi_mask`, `time_idx`, `parcel_idx`
    subset = function(...) {
      allowed <- c("roi_mask", "time_idx", "parcel_idx")
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
        if (is.null(self$h5) || !self$h5$is_valid) {
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
      run_header <- .read_lna_header(h5, run_id = self$current_run_id)
      handle <- DataHandle$new(
        h5 = h5,
        initial_meta = list(
          qpp_tile_cache = self$qpp_tile_cache,
          header = run_header
        ),
        subset = params,
        run_ids = self$run_ids,
        current_run_id = self$current_run_id
      )
      tf_group <- h5[["transforms"]]
      on.exit(safe_h5_close(tf_group), add = TRUE)
      transforms <- discover_transforms(tf_group)
      transforms <- select_transforms_for_run(
        transforms,
        tf_group,
        self$current_run_id
      )
      descriptor_map <- lapply(transforms$name, function(name) {
        read_json_descriptor(tf_group, name)
      })
      names(descriptor_map) <- transforms$name
      meta <- handle$meta
      meta$transform_descriptors <- descriptor_map
      handle <- handle$with(meta = meta)
      allow_plugins <- normalize_allow_plugins(self$allow_plugins)
      if (nrow(transforms) > 0) {
        missing_methods <- transforms$type[
          vapply(
            transforms$type,
            function(t) is.null(.resolve_transform_method("invert_step", t)),
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
    },

    #' @description
    #' Extract time series for specific voxels.
    #'
    #' Runs the inverse pipeline with an ROI mask selecting only the
    #' requested voxels. Returns a time x voxels matrix.
    #'
    #' @param voxel_indices Integer vector of linear voxel indices (in
    #'   mask order), or a logical 3D array used as an ROI mask.
    #' @param time_idx Optional integer vector of time indices.
    #' @return Numeric matrix (time x selected voxels).
    series = function(voxel_indices, time_idx = NULL) {
      if (is.null(self$h5) || !self$h5$is_valid) {
        abort_lna("lna_reader is closed",
                  .subclass = "lna_error_closed_reader",
                  location = "lna_reader:series")
      }

      params <- list()
      if (is.logical(voxel_indices) || (is.array(voxel_indices) && is.logical(voxel_indices))) {
        params$roi_mask <- voxel_indices
      } else {
        # Convert linear indices to a logical mask
        # Need the mask dimensions from the file
        handle_for_mask <- self$data(roi_mask = NULL, time_idx = 1L)
        mask_arr <- handle_for_mask$mask_info$mask
        if (is.null(mask_arr)) {
          abort_lna("No mask found; cannot convert voxel indices to ROI",
                    .subclass = "lna_error_missing_data",
                    location = "lna_reader:series")
        }
        roi <- array(FALSE, dim = dim(mask_arr))
        # voxel_indices are in masked-voxel order
        mask_idx <- which(as.logical(mask_arr))
        if (max(voxel_indices) > length(mask_idx)) {
          abort_lna("voxel_indices exceed number of masked voxels",
                    .subclass = "lna_error_validation",
                    location = "lna_reader:series")
        }
        roi[mask_idx[voxel_indices]] <- TRUE
        params$roi_mask <- roi
      }
      if (!is.null(time_idx)) params$time_idx <- as.integer(time_idx)

      handle <- self$data(roi_mask = params$roi_mask, time_idx = params$time_idx)
      mat <- as_dense_mat(handle$stash$input)
      # Normalize to time x voxels. The stash may contain voxels x time
      # when ROI-subsetted. Use the requested voxel count to orient.
      if (!is.null(params$roi_mask)) {
        n_roi <- sum(as.logical(params$roi_mask))
        if (nrow(mat) == n_roi && ncol(mat) != n_roi) {
          mat <- t(mat)
        }
      }
      mat
    },

    #' @description
    #' Coerce the reader's data to a neuroim2 NeuroVec.
    #'
    #' Materializes data and returns a spatially-aware neuroim2 object.
    #' Returns `SparseNeuroVec` when an ROI mask is active, otherwise
    #' `DenseNeuroVec`.
    #'
    #' @param assume_identity Logical. Explicitly allow identity geometry when
    #'   the LNA header has no dimensions or affine. Defaults to `FALSE`.
    #' @param ... Subsetting parameters passed to `$data()`.
    #' @return A `neuroim2::SparseNeuroVec` or `neuroim2::DenseNeuroVec`.
    as_neurovec = function(assume_identity = FALSE, ...) {
      handle <- self$data(...)
      .as_neurovec_with_h5(
        handle,
        self$h5,
        assume_identity = assume_identity
      )
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
