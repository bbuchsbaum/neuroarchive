#' Quantization Transform - Forward Step
#'
#' Applies simple min-max quantization to numeric data.
#' @keywords internal
forward_step.quant <- function(type, desc, handle) {
  p <- desc$params %||% list()
  bits <- p$bits %||% 8

  input_key <- if (!is.null(desc$inputs)) desc$inputs[[1]] else "input"
  x <- handle$get_inputs(input_key)[[1]]
  rng <- range(as.numeric(x))
  scale <- (rng[2] - rng[1]) / (2^bits - 1)
  offset <- rng[1]
  q <- round((x - offset) / scale)
  storage.mode(q) <- "integer"

  run_id <- handle$current_run_id %||% "run-01"
  data_path <- paste0("/scans/", run_id, "/quantized")
  scale_path <- paste0("/scans/", run_id, "/quant_scale")
  offset_path <- paste0("/scans/", run_id, "/quant_offset")

  plan <- handle$plan
  step_index <- plan$next_index
  params_json <- jsonlite::toJSON(p, auto_unbox = TRUE)

  plan$add_descriptor(plan$get_next_filename("quant"), list(type = "quant", params = p))
  plan$add_payload(data_path, q)
  plan$add_dataset_def(data_path, "quantized", "quant", run_id,
                       as.integer(step_index), params_json,
                       data_path, "eager")
  plan$add_payload(scale_path, scale)
  plan$add_dataset_def(scale_path, "scale", "quant", run_id,
                       as.integer(step_index), params_json,
                       scale_path, "eager")
  plan$add_payload(offset_path, offset)
  plan$add_dataset_def(offset_path, "offset", "quant", run_id,
                       as.integer(step_index), params_json,
                       offset_path, "eager")

  handle$plan <- plan
  handle$update_stash(keys = input_key, new_values = list())
}

#' Quantization Transform - Inverse Step
#'
#' Reconstructs data from quantized representation stored in HDF5.
#' @keywords internal
invert_step.quant <- function(type, desc, handle) {
  run_id <- handle$current_run_id %||% "run-01"
  data_path <- paste0("/scans/", run_id, "/quantized")
  scale_path <- paste0("/scans/", run_id, "/quant_scale")
  offset_path <- paste0("/scans/", run_id, "/quant_offset")

  root <- handle$h5[["/"]]
  q <- h5_read(root, data_path)
  scale <- as.numeric(h5_read(root, scale_path))
  offset <- as.numeric(h5_read(root, offset_path))
  x <- q * scale + offset

  input_key <- if (!is.null(desc$inputs)) desc$inputs[[1]] else "input"
  handle$update_stash(keys = character(), new_values = setNames(list(x), input_key))
}
