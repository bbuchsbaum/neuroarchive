# Integer temporal lifting reference for the experimental multiresolution-run
# profile. This is deliberately clear R, not the production decoder.

.lift_safe_integer <- 2^53 - 1

.lift_check_integer <- function(value, name = "value") {
  if (!is.numeric(value) || any(!is.finite(value)) ||
      any(value != trunc(value)) || any(abs(value) > .lift_safe_integer)) {
    stop(name, " must contain exact integers in the binary64 safe domain",
         call. = FALSE)
  }
  invisible(TRUE)
}

.lift_check_result <- function(value, operation) {
  if (any(!is.finite(value)) || any(value != trunc(value)) ||
      any(abs(value) > .lift_safe_integer)) {
    stop("integer overflow in ", operation, call. = FALSE)
  }
  value
}

.lift53_forward_once <- function(value) {
  n <- length(value)
  if (n < 2L) return(as.numeric(value))
  coarse <- as.numeric(value[seq.int(1L, n, by = 2L)])
  odd <- as.numeric(value[seq.int(2L, n, by = 2L)])
  detail_count <- length(odd)
  left <- coarse[seq_len(detail_count)]
  right <- coarse[pmin(seq_len(detail_count) + 1L, length(coarse))]
  detail <- .lift_check_result(
    odd - floor((left + right) / 2), "CDF 5/3 predict"
  )
  left_detail <- detail[pmax(seq_along(coarse) - 1L, 1L)]
  right_detail <- detail[pmin(seq_along(coarse), detail_count)]
  coarse <- .lift_check_result(
    coarse + floor((left_detail + right_detail + 2) / 4),
    "CDF 5/3 update"
  )
  c(coarse, detail)
}

.lift53_inverse_once <- function(value) {
  n <- length(value)
  if (n < 2L) return(as.numeric(value))
  coarse_count <- (n + 1L) %/% 2L
  coarse <- as.numeric(value[seq_len(coarse_count)])
  detail <- as.numeric(value[(coarse_count + 1L):n])
  detail_count <- length(detail)
  left_detail <- detail[pmax(seq_along(coarse) - 1L, 1L)]
  right_detail <- detail[pmin(seq_along(coarse), detail_count)]
  even <- .lift_check_result(
    coarse - floor((left_detail + right_detail + 2) / 4),
    "inverse CDF 5/3 update"
  )
  left <- even[seq_len(detail_count)]
  right <- even[pmin(seq_len(detail_count) + 1L, length(even))]
  odd <- .lift_check_result(
    detail + floor((left + right) / 2), "inverse CDF 5/3 predict"
  )
  out <- numeric(n)
  out[seq.int(1L, n, by = 2L)] <- even
  out[seq.int(2L, n, by = 2L)] <- odd
  out
}

.lift_haar_forward_once <- function(value) {
  n <- length(value)
  if (n < 2L) return(as.numeric(value))
  even <- as.numeric(value[seq.int(1L, n, by = 2L)])
  odd <- as.numeric(value[seq.int(2L, n, by = 2L)])
  detail <- .lift_check_result(odd - even[seq_along(odd)], "Haar predict")
  coarse <- even
  coarse[seq_along(detail)] <- .lift_check_result(
    even[seq_along(detail)] + floor(detail / 2), "Haar update"
  )
  c(coarse, detail)
}

.lift_haar_inverse_once <- function(value) {
  n <- length(value)
  if (n < 2L) return(as.numeric(value))
  coarse_count <- (n + 1L) %/% 2L
  coarse <- as.numeric(value[seq_len(coarse_count)])
  detail <- as.numeric(value[(coarse_count + 1L):n])
  even <- coarse
  even[seq_along(detail)] <- .lift_check_result(
    coarse[seq_along(detail)] - floor(detail / 2), "inverse Haar update"
  )
  odd <- .lift_check_result(
    detail + even[seq_along(detail)], "inverse Haar predict"
  )
  out <- numeric(n)
  out[seq.int(1L, n, by = 2L)] <- even
  out[seq.int(2L, n, by = 2L)] <- odd
  out
}

lift_max_levels <- function(length) {
  length <- as.integer(length)
  levels <- 0L
  while (length >= 2L) {
    levels <- levels + 1L
    length <- (length + 1L) %/% 2L
  }
  levels
}

lift_band_layout <- function(length, levels) {
  length <- as.integer(length)
  levels <- min(as.integer(levels), lift_max_levels(length))
  if (levels < 0L) stop("levels must be non-negative", call. = FALSE)
  if (levels == 0L) {
    return(data.frame(
      kind = "coarse", level = 0L, start = 1L, end = length,
      length = length, temporal_support = length,
      stringsAsFactors = FALSE
    ))
  }
  sizes <- integer(levels)
  details <- integer(levels)
  current <- length
  for (level in seq_len(levels)) {
    sizes[[level]] <- current
    details[[level]] <- current %/% 2L
    current <- (current + 1L) %/% 2L
  }
  position <- current + 1L
  rows <- list(data.frame(
    kind = "coarse", level = levels, start = 1L, end = current,
    length = current, temporal_support = length,
    stringsAsFactors = FALSE
  ))
  for (level in rev(seq_len(levels))) {
    count <- details[[level]]
    rows[[length(rows) + 1L]] <- data.frame(
      kind = "detail", level = level, start = position,
      end = position + count - 1L, length = count,
      temporal_support = min(length, 2^(level + 1L)),
      stringsAsFactors = FALSE
    )
    position <- position + count
  }
  do.call(rbind, rows)
}

.lift_forward_vector <- function(value, levels, wavelet) {
  .lift_check_integer(value)
  levels <- min(as.integer(levels), lift_max_levels(length(value)))
  if (levels < 0L) stop("levels must be non-negative", call. = FALSE)
  out <- as.numeric(value)
  current <- length(out)
  forward_once <- switch(
    wavelet, cdf53 = .lift53_forward_once, haar = .lift_haar_forward_once,
    stop("wavelet must be 'cdf53' or 'haar'", call. = FALSE)
  )
  for (level in seq_len(levels)) {
    out[seq_len(current)] <- forward_once(out[seq_len(current)])
    current <- (current + 1L) %/% 2L
  }
  out
}

.lift_inverse_vector <- function(value, levels, wavelet) {
  .lift_check_integer(value, "coefficients")
  levels <- min(as.integer(levels), lift_max_levels(length(value)))
  if (levels < 0L) stop("levels must be non-negative", call. = FALSE)
  sizes <- integer(levels)
  current <- length(value)
  for (level in seq_len(levels)) {
    sizes[[level]] <- current
    current <- (current + 1L) %/% 2L
  }
  inverse_once <- switch(
    wavelet, cdf53 = .lift53_inverse_once, haar = .lift_haar_inverse_once,
    stop("wavelet must be 'cdf53' or 'haar'", call. = FALSE)
  )
  out <- as.numeric(value)
  for (size in rev(sizes)) {
    out[seq_len(size)] <- inverse_once(out[seq_len(size)])
  }
  out
}

lift_forward_unit <- function(series, levels = 3L,
                              wavelet = c("cdf53", "haar"),
                              unit_id = "unit-01") {
  wavelet <- match.arg(wavelet)
  series <- as.matrix(series)
  .lift_check_integer(series, "series")
  if (any(series < -2147483648) || any(series > 2147483647)) {
    stop("series must be in the signed-int32 QPP source domain", call. = FALSE)
  }
  applied_levels <- min(as.integer(levels), lift_max_levels(nrow(series)))
  coefficients <- apply(series, 2L, .lift_forward_vector,
                        levels = applied_levels, wavelet = wavelet)
  if (is.null(dim(coefficients))) {
    coefficients <- matrix(coefficients, nrow = nrow(series), ncol = ncol(series))
  }
  structure(list(
    coefficients = unname(coefficients), unit_id = as.character(unit_id),
    unit_length = nrow(series), parcel_count = ncol(series),
    levels = applied_levels, wavelet = wavelet,
    boundary_extension = "whole-sample symmetric endpoint duplication",
    rounding = if (identical(wavelet, "cdf53")) {
      "JPEG-2000 reversible 5/3 floor predict and floor((d_left+d_right+2)/4) update"
    } else {
      "integer Haar detail=odd-even and coarse=even+floor(detail/2)"
    },
    coefficient_order = "coarse-L, detail-L, ..., detail-1",
    bands = lift_band_layout(nrow(series), applied_levels),
    independence_unit = as.character(unit_id),
    lossy_support = "complete-independence-unit"
  ), class = c("lift_integer_unit", "list"))
}

lift_inverse_unit <- function(object, coefficients = object$coefficients) {
  if (!inherits(object, "lift_integer_unit")) {
    stop("object must be a lift_integer_unit", call. = FALSE)
  }
  coefficients <- as.matrix(coefficients)
  if (!identical(dim(coefficients),
                 c(object$unit_length, object$parcel_count))) {
    stop("coefficient dimensions contradict the unit descriptor", call. = FALSE)
  }
  reconstructed <- apply(
    coefficients, 2L, .lift_inverse_vector,
    levels = object$levels, wavelet = object$wavelet
  )
  if (is.null(dim(reconstructed))) {
    reconstructed <- matrix(
      reconstructed, nrow = object$unit_length, ncol = object$parcel_count
    )
  }
  unname(reconstructed)
}

lift_forward_partitioned <- function(units, levels = 3L,
                                     wavelet = c("cdf53", "haar")) {
  wavelet <- match.arg(wavelet)
  if (!is.list(units) || !length(units)) {
    stop("units must be a non-empty list of matrices", call. = FALSE)
  }
  ids <- names(units)
  if (is.null(ids) || any(!nzchar(ids))) ids <- sprintf("unit-%02d", seq_along(units))
  transformed <- Map(function(unit, id) {
    lift_forward_unit(unit, levels = levels, wavelet = wavelet, unit_id = id)
  }, units, ids)
  parcel_counts <- vapply(transformed, `[[`, integer(1), "parcel_count")
  if (length(unique(parcel_counts)) != 1L) {
    stop("all independence units must have the same parcel count", call. = FALSE)
  }
  structure(list(
    units = stats::setNames(transformed, ids), wavelet = wavelet,
    unit_ids = ids, independence_unit = "declared-list-element",
    cross_unit_operations = FALSE,
    lossy_support = "complete-independence-unit"
  ), class = c("lift_integer_partitioned", "list"))
}

lift_inverse_partitioned <- function(object) {
  if (!inherits(object, "lift_integer_partitioned")) {
    stop("object must be a lift_integer_partitioned", call. = FALSE)
  }
  lapply(object$units, lift_inverse_unit)
}

.lift_detail_step <- function(detail_steps, level) {
  if (length(detail_steps) == 1L) return(as.numeric(detail_steps))
  if (!is.null(names(detail_steps)) && as.character(level) %in% names(detail_steps)) {
    return(as.numeric(detail_steps[[as.character(level)]]))
  }
  if (level <= length(detail_steps)) return(as.numeric(detail_steps[[level]]))
  stop("detail_steps does not define every retained level", call. = FALSE)
}

lift_quantize_bands <- function(object, detail_steps = 1,
                                dead_zone = 0, drop_finest = 0L) {
  if (!inherits(object, "lift_integer_unit")) {
    stop("object must be a lift_integer_unit", call. = FALSE)
  }
  if (!is.numeric(detail_steps) || any(!is.finite(detail_steps)) ||
      any(detail_steps < 1) || any(detail_steps != trunc(detail_steps))) {
    stop("detail_steps must be positive integers", call. = FALSE)
  }
  if (!is.numeric(dead_zone) || length(dead_zone) != 1L ||
      !is.finite(dead_zone) || dead_zone < 0) {
    stop("dead_zone must be one finite non-negative value", call. = FALSE)
  }
  drop_finest <- as.integer(drop_finest)
  if (drop_finest < 0L || drop_finest > object$levels) {
    stop("drop_finest must be in 0..levels", call. = FALSE)
  }
  quantized <- object$coefficients
  reconstructed <- object$coefficients
  band_steps <- numeric(object$levels)
  for (level in seq_len(object$levels)) {
    band <- object$bands[
      object$bands$kind == "detail" & object$bands$level == level, , drop = FALSE
    ]
    index <- band$start:band$end
    step <- .lift_detail_step(detail_steps, level)
    band_steps[[level]] <- step
    if (level <= drop_finest) {
      quantized[index, ] <- 0
      reconstructed[index, ] <- 0
    } else {
      value <- object$coefficients[index, , drop = FALSE]
      q <- round(value / step)
      q[abs(value) <= dead_zone * step] <- 0
      q[q == 0] <- 0
      quantized[index, ] <- q
      reconstructed[index, ] <- q * step
    }
  }
  structure(list(
    transform = object, quantized_coefficients = quantized,
    reconstructed_coefficients = reconstructed,
    coarse_step = 1, detail_steps = band_steps,
    dead_zone = dead_zone, dropped_detail_levels = seq_len(drop_finest),
    coefficient_order = object$coefficient_order,
    lossy_support = "complete-independence-unit"
  ), class = c("lift_quantized_unit", "list"))
}

lift_progressive_coefficients <- function(object, detail_levels = integer()) {
  if (!inherits(object, "lift_quantized_unit")) {
    stop("object must be a lift_quantized_unit", call. = FALSE)
  }
  detail_levels <- sort(unique(as.integer(detail_levels)))
  if (anyNA(detail_levels) || any(detail_levels < 1L) ||
      any(detail_levels > object$transform$levels)) {
    stop("detail_levels must be a subset of 1..levels", call. = FALSE)
  }
  out <- matrix(0, object$transform$unit_length, object$transform$parcel_count)
  coarse <- object$transform$bands[object$transform$bands$kind == "coarse", ]
  out[coarse$start:coarse$end, ] <-
    object$reconstructed_coefficients[coarse$start:coarse$end, , drop = FALSE]
  for (level in detail_levels) {
    band <- object$transform$bands[
      object$transform$bands$kind == "detail" &
        object$transform$bands$level == level, , drop = FALSE
    ]
    out[band$start:band$end, ] <-
      object$reconstructed_coefficients[band$start:band$end, , drop = FALSE]
  }
  out
}

lift_inverse_quantized <- function(object,
                                   detail_levels = seq_len(object$transform$levels)) {
  coefficients <- lift_progressive_coefficients(object, detail_levels)
  lift_inverse_unit(object$transform, coefficients)
}
