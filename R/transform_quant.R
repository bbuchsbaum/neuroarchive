#' Quantization Transform - Forward Step
#'
#' Implements the writer-side quant step. Parameters follow the quant
#' schema (`bits`, `method`, `center`, `scale_scope`, `allow_clip`).
#' Excessive clipping triggers an error unless `allow_clip` is `TRUE`.
#'
#' @param type,desc,handle Internal arguments used by the transform
#'   dispatcher.
#' @return None. Updates the write plan and `handle$meta`.
#' @keywords internal
forward_step.quant <- function(type, desc, handle) {
  opts <- desc$params %||% list()
  bits <- opts$bits %||% 8
  method <- opts$method %||% "range"
  center <- opts$center %||% TRUE
  scope <- opts$scale_scope %||% "global"
  snr_sample_frac <- opts$snr_sample_frac
  if (is.null(snr_sample_frac)) {
    snr_sample_frac <- lna_options("quant")[[1]]$snr_sample_frac
  }
  if (is.null(snr_sample_frac)) snr_sample_frac <- 0.01


  if (!(is.numeric(bits) && length(bits) == 1 &&
        bits == as.integer(bits) && bits >= 1 && bits <= 16)) {
    abort_lna(
      "bits must be an integer between 1 and 16",
      .subclass = "lna_error_validation",
      location = "forward_step.quant:bits"
    )
  }
  if (!(is.character(method) && length(method) == 1 &&
        method %in% c("range", "sd"))) {
    abort_lna(
      sprintf("Invalid method '%s'", method),
      .subclass = "lna_error_validation",
      location = "forward_step.quant:method"
    )
  }
  if (!(is.logical(center) && length(center) == 1)) {
    abort_lna(
      "center must be a single logical",
      .subclass = "lna_error_validation",
      location = "forward_step.quant:center"
    )
  }
  if (!(is.character(scope) && length(scope) == 1 &&
        scope %in% c("global", "voxel"))) {
    abort_lna(
      sprintf("Invalid scale_scope '%s'", scope),
      .subclass = "lna_error_validation",
      location = "forward_step.quant:scale_scope"
    )
  }
  if (!(is.numeric(snr_sample_frac) && length(snr_sample_frac) == 1 &&
        snr_sample_frac > 0 && snr_sample_frac <= 1)) {
    abort_lna(
      "snr_sample_frac must be a number in (0,1]",
      .subclass = "lna_error_validation",
      location = "forward_step.quant:snr_sample_frac"
    )
  }

  input_key <- if (!is.null(desc$inputs)) desc$inputs[[1]] else "input"
  input_data <- handle$get_inputs(input_key)[[1]]

  if (any(!is.finite(input_data))) {
    abort_lna(
      "quant cannot handle non-finite values – run an imputation/filtering transform first.",
      .subclass = "lna_error_validation",
      location = "forward_step.quant"
    )
  }

  if (identical(scope, "voxel")) {
    if (length(dim(input_data)) < 4) {
      warning("scale_scope='voxel' requires 4D data; falling back to global")
      scope <- "global"
    }
  }

  if (identical(scope, "voxel")) {
    params <- .quantize_voxel(input_data, bits, method, center)
  } else {
    params <- .quantize_global(as.numeric(input_data), bits, method, center)
    params$q <- array(params$q, dim = dim(input_data))
  }

  quantized_vals <- params$q
  scale <- params$scale
  offset <- params$offset
  if (identical(scope, "global")) {
    clip_warn_pct <- lna_options("quant.clip_warn_pct")[[1]]
    if (is.null(clip_warn_pct)) clip_warn_pct <- 0.5
    clip_abort_pct <- lna_options("quant.clip_abort_pct")[[1]]
    if (is.null(clip_abort_pct)) clip_abort_pct <- 5.0
    allow_clip <- isTRUE(opts$allow_clip)
    clip_pct <- params$clip_pct %||% 0
    if (clip_pct > clip_abort_pct && !allow_clip) {
      abort_lna(
        sprintf(
          "Clipping %.2f%% exceeds abort threshold of %.1f%%",
          clip_pct, clip_abort_pct
        ),
        .subclass = "lna_error_validation",
        location = "forward_step.quant:clipping"
      )
    } else if (clip_pct > clip_warn_pct) {
      warning(
        sprintf(
          "Clipping %.2f%% exceeds warning threshold %.1f%%",
          clip_pct, clip_warn_pct
        ),
        call. = FALSE
      )
    }
    handle$meta$quant_stats <- list(
      n_clipped_total = params$n_clipped_total %||% 0L,
      clip_pct = params$clip_pct %||% 0,
      scale_val = as.numeric(scale),
      offset_val = as.numeric(offset)
    )
    quantized_vals[quantized_vals < 0] <- 0L
    quantized_vals[quantized_vals > (2^bits - 1)] <- (2^bits - 1L)
  }

  storage.mode(quantized_vals) <- "integer"
  storage_type_str <- if (bits <= 8) "uint8" else "uint16"

  run_id <- handle$current_run_id %||% "run-01"
  run_id <- sanitize_run_id(run_id)
  data_path <- paste0("/scans/", run_id, "/quantized")
  scale_path <- paste0("/scans/", run_id, "/quant_scale")
  offset_path <- paste0("/scans/", run_id, "/quant_offset")

  blockwise <- identical(scope, "voxel") && !is.null(handle$h5) && handle$h5$is_valid
  if (blockwise) {
    dims <- dim(input_data)
    bs <- auto_block_size(dims[1:3],
                          element_size_bytes = if (bits <= 8) 1L else 2L)
    data_chunk <- c(bs$slab_dims, dims[4])
    param_chunk <- bs$slab_dims

    root <- handle$h5[["/"]]
    h5_create_empty_dataset(root, data_path, dims,
                            dtype = storage_type_str,
                            chunk_dims = data_chunk)
    dset_q <- root[[data_path]]
    h5_attr_write(dset_q, "quant_bits", as.integer(bits))
    h5_create_empty_dataset(root, scale_path, dims[1:3],
                            dtype = "float32",
                            chunk_dims = param_chunk)
    h5_create_empty_dataset(root, offset_path, dims[1:3],
                            dtype = "float32",
                            chunk_dims = param_chunk)
    dset_scale <- root[[scale_path]]
    dset_offset <- root[[offset_path]]

    slab <- bs$slab_dims
    n_clipped_total <- 0L
    sc_min <- Inf; sc_max <- -Inf; sc_sum <- 0; sc_sum_sq <- 0
    off_min <- Inf; off_max <- -Inf; off_sum <- 0; off_sum_sq <- 0
    for (z in seq(1, dims[3], by = slab[3])) {
      zi <- z:min(z + slab[3] - 1, dims[3])
      for (y in seq(1, dims[2], by = slab[2])) {
        yi <- y:min(y + slab[2] - 1, dims[2])
        for (x0 in seq(1, dims[1], by = slab[1])) {
          xi <- x0:min(x0 + slab[1] - 1, dims[1])
          block <- input_data[xi, yi, zi, , drop = FALSE]
          res <- .quantize_voxel_block(block, bits, method, center)
          idx <- list(xi, yi, zi, seq_len(dims[4]))
          dset_q$write(args = idx, res$q)
          dset_scale$write(args = idx[1:3], res$scale)
          dset_offset$write(args = idx[1:3], res$offset)
          n_clipped_total <- n_clipped_total + res$n_clipped_total
          sc <- as.numeric(res$scale)
          of <- as.numeric(res$offset)
          sc_min <- min(sc_min, min(sc))
          sc_max <- max(sc_max, max(sc))
          sc_sum <- sc_sum + sum(sc)
          sc_sum_sq <- sc_sum_sq + sum(sc^2)
          off_min <- min(off_min, min(of))
          off_max <- max(off_max, max(of))
          off_sum <- off_sum + sum(of)
          off_sum_sq <- off_sum_sq + sum(of^2)
        }
      }
    }
    dset_q$close(); dset_scale$close(); dset_offset$close()
    clip_pct <- if (length(input_data) > 0) 100 * n_clipped_total / length(input_data) else 0
    vox_total <- prod(dims[1:3])
    sc_mean <- sc_sum / vox_total
    sc_sd <- sqrt(sc_sum_sq / vox_total - sc_mean^2)
    off_mean <- off_sum / vox_total
    off_sd <- sqrt(off_sum_sq / vox_total - off_mean^2)
    handle$meta$quant_stats <- list(
      n_clipped_total = as.integer(n_clipped_total),
      clip_pct = clip_pct,
      scale_min = sc_min,
      scale_max = sc_max,
      scale_mean = sc_mean,
      scale_sd = sc_sd,
      offset_min = off_min,
      offset_max = off_max,
      offset_mean = off_mean,
      offset_sd = off_sd
    )
    quantized_vals <- NULL
    scale <- NULL
    offset <- NULL
  }

  input_range <- range(as.numeric(input_data))
  qs <- handle$meta$quant_stats
  if (!identical(scope, "voxel")) {
    var_x <- stats::var(as.numeric(input_data))
    snr_db <- 10 * log10(var_x / ((qs$scale_val)^2 / 12))
    hist_info <- hist(as.numeric(quantized_vals), breaks = 64, plot = FALSE)
    quant_report <- list(
      report_version = "1.0",
      clipped_samples_count = qs$n_clipped_total,
      clipped_samples_percentage = qs$clip_pct,
      input_data_range = input_range,
      effective_step_size = qs$scale_val,
      effective_offset = qs$offset_val,
      estimated_snr_db = snr_db,
      histogram_quantized_values = list(
        breaks = hist_info$breaks,
        counts = unname(hist_info$counts)
      )
    )
  } else {
    dims <- dim(input_data)
    vox <- prod(dims[1:3])
    # Sample only a fraction of voxels to estimate SNR quickly when the
    # volume is large. This keeps memory usage and runtime manageable.
    sample_n <- max(1L, floor(vox * snr_sample_frac))
    mat <- matrix(as.numeric(input_data), vox, dims[4])
    idx <- sample(vox, sample_n)
    var_x <- stats::var(as.numeric(mat[idx, , drop = FALSE]))
    snr_db <- 10 * log10(var_x / ((qs$scale_mean)^2 / 12))
    quant_report <- list(
      report_version = "1.0",
      clipped_samples_count = qs$n_clipped_total,
      clipped_samples_percentage = qs$clip_pct,
      input_data_range = input_range,
      scale_stats = list(
        min = qs$scale_min,
        max = qs$scale_max,
        mean = qs$scale_mean,
        sd = qs$scale_sd
      ),
      offset_stats = list(
        min = qs$offset_min,
        max = qs$offset_max,
        mean = qs$offset_mean,
        sd = qs$offset_sd
      ),
      estimated_snr_db = snr_db
    )
  }
  handle$meta$quant_report <- quant_report

  plan <- handle$plan
  step_index <- plan$next_index
  fname <- plan$get_next_filename("quant")
  base_name <- tools::file_path_sans_ext(fname)

  report_path <- paste0("/transforms/", base_name, "_report.json")
  opts$report_path <- report_path
  json_report_str <- jsonlite::toJSON(quant_report, auto_unbox = TRUE, pretty = TRUE)
  gzipped_report <- memCompress(charToRaw(json_report_str), type = "gzip")

  payload_key_report <- report_path
  if (!is.null(handle$h5) && handle$h5$is_valid) {
    root <- handle$h5[["/"]]
    tryCatch({
      h5_write_dataset(root, report_path, gzipped_report, dtype = "uint8")
      dset_rep <- root[[report_path]]
      h5_attr_write(dset_rep, "compression", "gzip")
      dset_rep$close()
    }, error = function(e) {
      stop(paste0("Error writing quantization report: ", conditionMessage(e)), call. = FALSE)
    })
    payload_key_report <- ""
  } else {
    plan$add_payload(report_path, gzipped_report)
  }
  params_json <- as.character(jsonlite::toJSON(opts, auto_unbox = TRUE))


  # Use the existing desc object which has inputs/outputs set by run_transform_loop
  desc$type <- "quant"
  desc$params <- opts
  desc$version <- "1.0"
  
  plan$add_descriptor(fname, desc)
  payload_key_data <- data_path
  payload_key_scale <- scale_path
  payload_key_offset <- offset_path
  if (blockwise) {
    payload_key_data <- ""
    payload_key_scale <- ""
    payload_key_offset <- ""
  } else {
    plan$add_payload(data_path, quantized_vals)
    plan$add_payload(scale_path, scale)
    plan$add_payload(offset_path, offset)
  }
  plan$add_dataset_def(data_path, "quantized", as.character(type), run_id,
                       as.integer(step_index), params_json,
                       payload_key_data, "eager", dtype = storage_type_str)
  plan$add_dataset_def(scale_path, "scale", as.character(type), run_id,
                       as.integer(step_index), params_json,
                       payload_key_scale, "eager", dtype = "float32")
  plan$add_dataset_def(offset_path, "offset", as.character(type), run_id,
                       as.integer(step_index), params_json,
                       payload_key_offset, "eager", dtype = "float32")
  
  # Only add report dataset definition if it wasn't written directly
  if (nzchar(payload_key_report)) {
    plan$add_dataset_def(report_path, "quant_report", as.character(type), run_id,
                         as.integer(step_index), params_json,
                         payload_key_report, "eager", dtype = "uint8")
  }

  handle$plan <- plan
  handle <- handle$update_stash(keys = input_key, new_values = list())
  return(handle)
}

#' Quantization Transform - Inverse Step
#'
#' Reads quantized values and reconstructs the original data using the
#' stored scale and offset. Validates the `quant_bits` attribute and
#' handles legacy files missing it.
#'
#' @param type,desc,handle Internal arguments used by the transform
#'   dispatcher.
#' @return None. Places the reconstructed array into `handle$stash`.
#' @keywords internal
invert_step.quant <- function(type, desc, handle) {
  run_id <- handle$current_run_id %||% "run-01"
  run_id <- sanitize_run_id(run_id)
  data_path <- paste0("/scans/", run_id, "/quantized")
  scale_path <- paste0("/scans/", run_id, "/quant_scale")
  offset_path <- paste0("/scans/", run_id, "/quant_offset")

  root <- handle$h5[["/"]]
  dset <- NULL
  q <- NULL
  attr_bits <- desc$params$bits %||% NA
  tryCatch({
    dset <- root[[data_path]]
    if (h5_attr_exists(dset, "quant_bits")) {
      attr_bits <- h5_attr_read(dset, "quant_bits")
    } else {
      warning("quant_bits HDF5 attribute missing; using descriptor value.",
              call. = FALSE)
    }
    q <- dset$read()
  }, error = function(e) {
    stop(paste0("Error reading dataset '", data_path, "': ",
                conditionMessage(e)), call. = FALSE)
  }, finally = {
    if (!is.null(dset) && inherits(dset, "H5D")) dset$close()
  })
  if (!is.null(desc$params$bits) && !is.na(desc$params$bits)) {
    if (!is.na(attr_bits) && attr_bits != desc$params$bits) {
      abort_lna(
        sprintf(
          "quant_bits attribute (%s) disagrees with descriptor bits (%s)",
          attr_bits, desc$params$bits
        ),
        .subclass = "lna_error_validation",
        location = "invert_step.quant:bits"
      )
    }
  }
  scale <- as.numeric(h5_read(root, scale_path))
  offset <- as.numeric(h5_read(root, offset_path))
  x <- q * scale + offset

  # Subsetting logic from original could be applied here if needed
  input_key <- if (!is.null(desc$inputs) && length(desc$inputs) > 0) desc$inputs[[1]] else "input"
  
  handle <- handle$update_stash(keys = character(), new_values = setNames(list(x), input_key))
  return(handle)
}

#' Compute quantization parameters globally
#' @keywords internal
.quantize_global <- function(x, bits, method, center) {
  stopifnot(is.numeric(x))
  if (any(!is.finite(x))) {
    abort_lna(
      "non-finite values found in input",
      .subclass = "lna_error_validation",
      location = ".quantize_global"
    )
  }
  rng <- range(x)
  m <- mean(x)
  s <- stats::sd(x)

  if (center) {
    if (identical(method, "sd")) {
      max_abs <- 3 * s
    } else {
      max_abs <- max(abs(rng - m))
    }
    scale <- (2 * max_abs) / (2^bits - 1)
    offset <- m - max_abs
  } else {
    if (identical(method, "sd")) {
      lo <- m - 3 * s
      hi <- m + 3 * s
    } else {
      lo <- rng[1]
      hi <- rng[2]
    }
    scale <- (hi - lo) / (2^bits - 1)
    offset <- lo
  }

  if (scale == 0) {
    # Zero variance input: all values are identical. Force scale to 1 so
    # quantized output is defined (all zeros).
    scale <- 1
    q_raw <- rep.int(0L, length(x))
    n_clipped_total <- 0L
  } else {
    q_raw <- round((x - offset) / scale)
    n_clipped_total <- sum(q_raw < 0 | q_raw > 2^bits - 1)
  }

  clip_pct <- if (length(x) > 0) {
    100 * n_clipped_total / length(x)
  } else {
    0
  }

  q <- q_raw
  q[q < 0] <- 0L
  q[q > 2^bits - 1] <- 2^bits - 1L

  list(q = q,
       scale = scale,
       offset = offset,
       n_clipped_total = as.integer(n_clipped_total),
       clip_pct = clip_pct)
}

#' Quantize voxel-wise time series
#'
#' Internal helper that performs quantization for each voxel. When
#' `collect_clip` is `TRUE`, the number of clipped samples is also returned.
#' @keywords internal
.quantize_voxel_core <- function(x, bits, method, center, collect_clip = FALSE) {
  if (any(!is.finite(x))) {
    abort_lna(
      "non-finite values found in input",
      .subclass = "lna_error_validation",
      location = if (collect_clip) ".quantize_voxel_block" else ".quantize_voxel"
    )
  }

  dims <- dim(x)
  vox <- prod(dims[1:3])
  time <- dims[4]
  mat <- matrix(as.numeric(x), vox, time)

  m <- rowMeans(mat)
  s <- apply(mat, 1, stats::sd)
  rng_lo <- apply(mat, 1, min)
  rng_hi <- apply(mat, 1, max)

  if (center) {
    max_abs <- if (identical(method, "sd")) 3 * s else pmax(abs(rng_hi - m), abs(rng_lo - m))
    scale <- (2 * max_abs) / (2^bits - 1)
    offset <- m - max_abs
  } else {
    if (identical(method, "sd")) {
      lo <- m - 3 * s
      hi <- m + 3 * s
    } else {
      lo <- rng_lo
      hi <- rng_hi
    }
    scale <- (hi - lo) / (2^bits - 1)
    offset <- lo
  }

  zero_idx <- scale == 0
  q <- matrix(0L, vox, time)
  nclip <- 0L
  if (any(!zero_idx)) {
    q_raw <- round((mat[!zero_idx, , drop = FALSE] - offset[!zero_idx]) / scale[!zero_idx])
    if (collect_clip) nclip <- sum(q_raw < 0 | q_raw > 2^bits - 1)
    q[!zero_idx, ] <- q_raw
    q[q < 0] <- 0L
    q[q > 2^bits - 1] <- 2^bits - 1L
  }
  scale[zero_idx] <- 1

  result <- list(
    q = array(as.integer(q), dim = dims),
    scale = array(scale, dim = dims[1:3]),
    offset = array(offset, dim = dims[1:3])
  )
  if (collect_clip) result$n_clipped_total <- as.integer(nclip)
  result
}

#' Compute quantization parameters per voxel (time series)
#' @keywords internal
.quantize_voxel <- function(x, bits, method, center) {
  .quantize_voxel_core(x, bits, method, center, collect_clip = FALSE)
}

#' Quantize a voxel block with clipping stats
#'
#' Helper used by block-wise processing to quantize a subset of voxels and
#' return clipping information along with scale/offset parameters.
#' @keywords internal
.quantize_voxel_block <- function(x, bits, method, center) {
  .quantize_voxel_core(x, bits, method, center, collect_clip = TRUE)
}
