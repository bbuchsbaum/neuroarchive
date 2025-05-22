#' Quantization Transform - Forward Step
#'
#' Applies simple min-max quantization to numeric data.
#' @keywords internal
forward_step.quant <- function(type, desc, handle) {
  p <- desc$params %||% list()
  bits <- p$bits %||% 8
  method <- p$method %||% "range"
  center <- p$center %||% TRUE
  scope <- p$scale_scope %||% "global"


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

  input_key <- if (!is.null(desc$inputs)) desc$inputs[[1]] else "input"
  x <- handle$get_inputs(input_key)[[1]]

  if (any(!is.finite(x))) {
    abort_lna(
      "quant cannot handle non-finite values – run an imputation/filtering transform first.",
      .subclass = "lna_error_validation",
      location = "forward_step.quant"
    )
  }

  if (identical(scope, "voxel")) {
    if (length(dim(x)) < 4) {
      warning("scale_scope='voxel' requires 4D data; falling back to global")
      scope <- "global"
    }
  }

  if (identical(scope, "voxel")) {
    params <- .quantize_voxel(x, bits, method, center)
  } else {
    params <- .quantize_global(as.numeric(x), bits, method, center)
    params$q <- array(params$q, dim = dim(x))
  }

  q <- params$q
  scale <- params$scale
  offset <- params$offset
  if (identical(scope, "global")) {
    clip_warn_pct <- lna_options("quant.clip_warn_pct")[[1]]
    if (is.null(clip_warn_pct)) clip_warn_pct <- 0.5
    clip_abort_pct <- lna_options("quant.clip_abort_pct")[[1]]
    if (is.null(clip_abort_pct)) clip_abort_pct <- 5.0
    allow_clip <- isTRUE(p$allow_clip)
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
    q[q < 0] <- 0L
    q[q > (2^bits - 1)] <- (2^bits - 1L)
  }

  storage.mode(q) <- "integer"
  storage_type_str <- if (bits <= 8) "uint8" else "uint16"

  run_id <- handle$current_run_id %||% "run-01"
  run_id <- sanitize_run_id(run_id)
  data_path <- paste0("/scans/", run_id, "/quantized")
  scale_path <- paste0("/scans/", run_id, "/quant_scale")
  offset_path <- paste0("/scans/", run_id, "/quant_offset")

  blockwise <- identical(scope, "voxel") && !is.null(handle$h5) && handle$h5$is_valid
  if (blockwise) {
    dims <- dim(x)
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
          block <- x[xi, yi, zi, , drop = FALSE]
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
    clip_pct <- if (length(x) > 0) 100 * n_clipped_total / length(x) else 0
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
    q <- NULL
    scale <- NULL
    offset <- NULL
  }

  input_range <- range(as.numeric(x))
  qs <- handle$meta$quant_stats
  if (!identical(scope, "voxel")) {
    var_x <- stats::var(as.numeric(x))
    snr_db <- 10 * log10(var_x / ((qs$scale_val)^2 / 12))
    hist_info <- hist(as.numeric(q), breaks = 64, plot = FALSE)
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
    var_x <- stats::var(as.numeric(x))
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
  params_json <- as.character(jsonlite::toJSON(p, auto_unbox = TRUE))

  plan$add_descriptor(plan$get_next_filename("quant"), list(type = "quant", params = p))
  payload_key_data <- data_path
  payload_key_scale <- scale_path
  payload_key_offset <- offset_path
  if (blockwise) {
    payload_key_data <- ""
    payload_key_scale <- ""
    payload_key_offset <- ""
  } else {
    plan$add_payload(data_path, q)
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

  handle$plan <- plan
  handle$update_stash(keys = input_key, new_values = list())
}

#' Quantization Transform - Inverse Step
#'
#' Reconstructs data from quantized representation stored in HDF5.
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

  subset <- handle$subset
  if (!is.null(subset$roi_mask)) {
    roi <- as.logical(subset$roi_mask)
    if (length(dim(x)) == 4 && length(dim(subset$roi_mask)) == 3 &&
        all(dim(subset$roi_mask) == dim(x)[1:3])) {
      vox_idx <- which(roi)
      mat <- matrix(as.numeric(x), prod(dim(x)[1:3]), dim(x)[4])
      x <- mat[vox_idx, , drop = FALSE]
    }
  }
  if (!is.null(subset$time_idx)) {
    idx <- as.integer(subset$time_idx)
    if (is.matrix(x)) {
      x <- x[, idx, drop = FALSE]
    } else if (length(dim(x)) == 4) {
      x <- x[,,, idx, drop = FALSE]
    } else {
      x <- x[idx]
    }
  }

  input_key <- if (!is.null(desc$inputs)) desc$inputs[[1]] else "input"
  handle$update_stash(keys = character(), new_values = setNames(list(x), input_key))
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

#' Compute quantization parameters per voxel (time series)
#' @keywords internal
.quantize_voxel <- function(x, bits, method, center) {
  if (any(!is.finite(x))) {
    abort_lna(
      "non-finite values found in input",
      .subclass = "lna_error_validation",
      location = ".quantize_voxel"
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
  if (any(!zero_idx)) {
    q[!zero_idx, ] <- round((mat[!zero_idx, , drop = FALSE] - offset[!zero_idx]) / scale[!zero_idx])
    q[q < 0] <- 0L
    q[q > 2^bits - 1] <- 2^bits - 1L
  }
  scale[zero_idx] <- 1

  arr_q <- array(as.integer(q), dim = dims)
  list(q = arr_q,
       scale = array(scale, dim = dims[1:3]),
       offset = array(offset, dim = dims[1:3]))
}

#' Quantize a voxel block with clipping stats
#'
#' Helper used by block-wise processing to quantize a subset of voxels and
#' return clipping information along with scale/offset parameters.
#' @keywords internal
.quantize_voxel_block <- function(x, bits, method, center) {
  if (any(!is.finite(x))) {
    abort_lna(
      "non-finite values found in input",
      .subclass = "lna_error_validation",
      location = ".quantize_voxel_block"
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
    nclip <- sum(q_raw < 0 | q_raw > 2^bits - 1)
    q[!zero_idx, ] <- q_raw
    q[q < 0] <- 0L
    q[q > 2^bits - 1] <- 2^bits - 1L
  }
  scale[zero_idx] <- 1

  arr_q <- array(as.integer(q), dim = dims)
  list(q = arr_q,
       scale = array(scale, dim = dims[1:3]),
       offset = array(offset, dim = dims[1:3]),
       n_clipped_total = as.integer(nclip))
}
