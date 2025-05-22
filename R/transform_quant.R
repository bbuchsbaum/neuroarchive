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

  storage.mode(q) <- "integer"
  storage_type_str <- if (bits <= 8) "uint8" else "uint16"

  run_id <- handle$current_run_id %||% "run-01"
  run_id <- sanitize_run_id(run_id)
  data_path <- paste0("/scans/", run_id, "/quantized")
  scale_path <- paste0("/scans/", run_id, "/quant_scale")
  offset_path <- paste0("/scans/", run_id, "/quant_offset")

  plan <- handle$plan
  step_index <- plan$next_index
  params_json <- as.character(jsonlite::toJSON(p, auto_unbox = TRUE))

  plan$add_descriptor(plan$get_next_filename("quant"), list(type = "quant", params = p))
  plan$add_payload(data_path, q)
  plan$add_dataset_def(data_path, "quantized", as.character(type), run_id,
                       as.integer(step_index), params_json,
                       data_path, "eager", dtype = storage_type_str)
  plan$add_payload(scale_path, scale)
  plan$add_dataset_def(scale_path, "scale", as.character(type), run_id,
                       as.integer(step_index), params_json,
                       scale_path, "eager", dtype = "float32")
  plan$add_payload(offset_path, offset)
  plan$add_dataset_def(offset_path, "offset", as.character(type), run_id,
                       as.integer(step_index), params_json,
                       offset_path, "eager", dtype = "float32")

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
    q <- rep.int(0L, length(x))
  } else {
    q <- round((x - offset) / scale)
    q[q < 0] <- 0L
    q[q > 2^bits - 1] <- 2^bits - 1L
  }

  list(q = q, scale = scale, offset = offset)
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
