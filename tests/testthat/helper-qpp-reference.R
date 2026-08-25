# Deliberately simple, independent QPP v1 numerical oracle.
# This is test code, not the production implementation.

.qpp_ref_q_min <- -(2^31)
.qpp_ref_q_max <- 2^31 - 1
.qpp_ref_residual_min <- -(2^33) + 2
.qpp_ref_residual_max <- 2^33 - 2

read_qpp_golden <- function() {
  path <- testthat::test_path("fixtures", "qpp-golden-v1.json")
  jsonlite::fromJSON(path, simplifyVector = FALSE)$fixtures
}

golden_q_matrix <- function(fixture) {
  columns <- lapply(fixture$q_by_column, as.numeric)
  do.call(cbind, columns)
}

.qpp_ref_assert_whole <- function(x, min = -Inf, max = Inf, what = "value") {
  if (!is.numeric(x) || any(!is.finite(x)) || any(x != trunc(x)) ||
      any(x < min) || any(x > max)) {
    stop(what, " must contain finite whole numbers in [", min, ", ", max, "]",
         call. = FALSE)
  }
  invisible(x)
}

.qpp_ref_quantize <- function(y, delta) {
  y <- as.matrix(y)
  if (!is.numeric(y) || any(!is.finite(y))) {
    stop("y must be a finite numeric matrix", call. = FALSE)
  }
  if (!is.numeric(delta) || length(delta) != ncol(y) ||
      any(!is.finite(delta)) || any(delta <= 0)) {
    stop("delta must contain one finite positive value per parcel", call. = FALSE)
  }
  scaled <- sweep(y, 2L, delta, "/")
  q <- round(scaled) # IEC 60559 round-to-nearest, ties-to-even in R.
  .qpp_ref_assert_whole(q, .qpp_ref_q_min, .qpp_ref_q_max, "quantized values")
  dim(q) <- dim(y)
  q
}

.qpp_ref_reconstruct <- function(q, delta, original = NULL) {
  q <- as.matrix(q)
  .qpp_ref_assert_whole(q, .qpp_ref_q_min, .qpp_ref_q_max, "q")
  if (!is.numeric(delta) || length(delta) != ncol(q) ||
      any(!is.finite(delta)) || any(delta <= 0)) {
    stop("delta must contain one finite positive value per parcel", call. = FALSE)
  }
  decoded <- sweep(q, 2L, delta, "*")
  out <- list(q = q, decoded = decoded, ideal_bound = delta / 2)
  if (!is.null(original)) {
    original <- as.matrix(original)
    if (!identical(dim(original), dim(decoded)) || any(!is.finite(original))) {
      stop("original must be finite and have the same dimensions as q", call. = FALSE)
    }
    err <- abs(decoded - original)
    out$maximum_absolute_error <- max(err)
    out$parcelwise_maximum_absolute_error <- apply(err, 2L, max)
    scale <- pmax(1, apply(abs(original), 2L, max), apply(abs(decoded), 2L, max))
    out$conservative_floating_bound <- out$ideal_bound + 8 * .Machine$double.eps * scale
    out$excess_over_ideal_bound <- pmax(
      out$parcelwise_maximum_absolute_error - out$ideal_bound,
      0
    )
  }
  out
}

.qpp_ref_residuals <- function(q, predictor) {
  q <- as.numeric(q)
  .qpp_ref_assert_whole(q, .qpp_ref_q_min, .qpp_ref_q_max, "q")
  predictor <- match.arg(as.character(predictor), c("P0", "P1", "P2"))
  n <- length(q)
  if (n < 1L) stop("q must not be empty", call. = FALSE)

  if (predictor == "P0") {
    anchors <- numeric()
    residuals <- q
  } else if (predictor == "P1") {
    anchors <- q[1L]
    residuals <- if (n == 1L) numeric() else q[-1L] - q[-n]
  } else {
    if (n < 2L) stop("P2 requires at least two samples", call. = FALSE)
    anchors <- q[1:2]
    residuals <- if (n == 2L) numeric() else {
      q[3:n] - 2 * q[2:(n - 1L)] + q[1:(n - 2L)]
    }
  }
  # Doubles exactly represent every residual produced by signed int32 q.
  .qpp_ref_assert_whole(
    residuals,
    .qpp_ref_residual_min,
    .qpp_ref_residual_max,
    "residuals"
  )
  list(anchors = anchors, residuals = residuals)
}

.qpp_ref_inverse_residuals <- function(anchors, residuals, predictor) {
  predictor <- match.arg(as.character(predictor), c("P0", "P1", "P2"))
  residuals <- as.numeric(residuals)
  if (predictor == "P0") return(residuals)
  if (predictor == "P1") {
    if (length(anchors) != 1L) stop("P1 requires one anchor", call. = FALSE)
    return(c(anchors, anchors + cumsum(residuals)))
  }
  if (length(anchors) != 2L) stop("P2 requires two anchors", call. = FALSE)
  q <- c(as.numeric(anchors), numeric(length(residuals)))
  if (length(residuals)) {
    for (i in seq_along(residuals)) {
      j <- i + 2L
      q[j] <- residuals[i] + 2 * q[j - 1L] - q[j - 2L]
    }
  }
  q
}

.qpp_ref_zigzag <- function(x) {
  .qpp_ref_assert_whole(x, .qpp_ref_residual_min, .qpp_ref_residual_max,
                        "signed values")
  if (!length(x)) return(numeric())
  ifelse(x >= 0, 2 * x, -2 * x - 1)
}

.qpp_ref_unzigzag <- function(x) {
  .qpp_ref_assert_whole(x, 0, 2^34 - 1, "zig-zag values")
  if (!length(x)) return(numeric())
  ifelse(x %% 2 == 0, x / 2, -(x + 1) / 2)
}

.qpp_ref_width <- function(unsigned) {
  .qpp_ref_assert_whole(unsigned, 0, 2^34 - 1, "unsigned values")
  if (!length(unsigned) || max(unsigned) == 0) return(0L)
  width <- 0L
  limit <- 1
  maximum <- max(unsigned)
  while (limit <= maximum) {
    width <- width + 1L
    limit <- limit * 2
  }
  width
}

.qpp_ref_pack <- function(unsigned, width) {
  .qpp_ref_assert_whole(unsigned, 0, 2^34 - 1, "unsigned values")
  .qpp_ref_assert_whole(width, 0, 34, "width")
  width <- as.integer(width)
  if (!length(unsigned) || width == 0L) {
    if (length(unsigned) && any(unsigned != 0)) stop("width zero only encodes zero", call. = FALSE)
    return(raw())
  }
  if (any(unsigned >= 2^width)) stop("value does not fit width", call. = FALSE)
  out <- integer(ceiling(length(unsigned) * width / 8))
  bit_cursor <- 0L
  for (value in unsigned) {
    for (bit in 0:(width - 1L)) {
      byte <- bit_cursor %/% 8L + 1L
      offset <- bit_cursor %% 8L
      out[byte] <- out[byte] + ((value %/% (2^bit)) %% 2) * (2^offset)
      bit_cursor <- bit_cursor + 1L
    }
  }
  as.raw(out)
}

.qpp_ref_unpack <- function(bytes, count, width) {
  .qpp_ref_assert_whole(count, 0, 2^31 - 1, "count")
  .qpp_ref_assert_whole(width, 0, 34, "width")
  count <- as.integer(count)
  width <- as.integer(width)
  needed <- ceiling(count * width / 8)
  if (length(bytes) != needed) stop("packed byte length mismatch", call. = FALSE)
  if (width == 0L) return(numeric(count))
  ans <- numeric(count)
  bit_cursor <- 0L
  b <- as.integer(bytes)
  for (i in seq_len(count)) {
    value <- 0
    for (bit in 0:(width - 1L)) {
      byte <- bit_cursor %/% 8L + 1L
      offset <- bit_cursor %% 8L
      value <- value + ((b[byte] %/% (2^offset)) %% 2) * (2^bit)
      bit_cursor <- bit_cursor + 1L
    }
    ans[i] <- value
  }
  ans
}

.qpp_ref_u16le <- function(x) {
  .qpp_ref_assert_whole(x, 0, 2^16 - 1, "uint16")
  as.raw(c(x %% 256, (x %/% 256) %% 256))
}

.qpp_ref_u32le <- function(x) {
  .qpp_ref_assert_whole(x, 0, 2^32 - 1, "uint32")
  as.raw(c(x %% 256, (x %/% 256) %% 256, (x %/% 2^16) %% 256,
           (x %/% 2^24) %% 256))
}

.qpp_ref_i32le <- function(x) {
  .qpp_ref_assert_whole(x, .qpp_ref_q_min, .qpp_ref_q_max, "int32")
  .qpp_ref_u32le(ifelse(x < 0, x + 2^32, x))
}

.qpp_ref_choose_predictor <- function(q) {
  candidates <- c("P0", "P1", if (length(q) >= 2L) "P2")
  plans <- lapply(candidates, function(predictor) {
    rr <- .qpp_ref_residuals(q, predictor)
    unsigned <- .qpp_ref_zigzag(rr$residuals)
    width <- .qpp_ref_width(unsigned)
    # Predictor record headers are equal size, so only anchors and payload vary.
    cost <- length(rr$anchors) * 32 + length(unsigned) * width
    list(predictor = predictor, predictor_id = match(predictor, c("P0", "P1", "P2")) - 1L,
         anchors = rr$anchors, unsigned = unsigned, width = width, cost_bits = cost)
  })
  costs <- vapply(plans, `[[`, numeric(1), "cost_bits")
  plans[[which.min(costs)]] # which.min gives canonical P0, P1, P2 tie order.
}

.qpp_ref_encode_tile <- function(q) {
  q <- as.matrix(q)
  if (nrow(q) < 1L || nrow(q) > 2^16 - 1 ||
      ncol(q) < 1L || ncol(q) > 2^16 - 1) {
    stop("tile dimensions must be in 1..65535", call. = FALSE)
  }
  .qpp_ref_assert_whole(q, .qpp_ref_q_min, .qpp_ref_q_max, "q")
  bytes <- c(charToRaw("QPP1"), .qpp_ref_u16le(nrow(q)), .qpp_ref_u16le(ncol(q)))
  plans <- vector("list", ncol(q))
  payload_bits <- 0
  payload_bytes <- 0

  for (k in seq_len(ncol(q))) {
    plan <- .qpp_ref_choose_predictor(q[, k])
    payload <- .qpp_ref_pack(plan$unsigned, plan$width)
    record <- c(
      as.raw(c(plan$predictor_id, plan$width, length(plan$anchors), 0L)),
      .qpp_ref_u16le(length(plan$unsigned)),
      .qpp_ref_u32le(length(payload))
    )
    if (length(plan$anchors)) {
      for (anchor in plan$anchors) record <- c(record, .qpp_ref_i32le(anchor))
    }
    bytes <- c(bytes, record, payload)
    plan$payload_nbyte <- length(payload)
    plan$residual_count <- length(plan$unsigned)
    plans[[k]] <- plan
    payload_bits <- payload_bits + length(plan$unsigned) * plan$width
    payload_bytes <- payload_bytes + length(payload)
  }
  total_bits <- length(bytes) * 8
  list(
    bytes = bytes,
    plan = plans,
    metrics = list(
      payload_bits = payload_bits,
      payload_storage_bits = payload_bytes * 8,
      header_and_anchor_bits = total_bits - payload_bytes * 8,
      total_bits = total_bits,
      effective_bits_per_parcel_tr = total_bits / length(q)
    )
  )
}

.qpp_ref_decode_tile <- function(encoded) {
  bytes <- if (is.list(encoded)) encoded$bytes else encoded
  if (!is.raw(bytes)) stop("encoded tile must be raw", call. = FALSE)
  pos <- 1L
  take <- function(n) {
    if (n < 0 || pos + n - 1L > length(bytes)) stop("truncated QPP tile", call. = FALSE)
    ans <- if (n == 0L) raw() else bytes[pos:(pos + n - 1L)]
    pos <<- pos + n
    ans
  }
  u8 <- function() as.integer(take(1L))
  u16 <- function() { b <- as.integer(take(2L)); b[1L] + 256 * b[2L] }
  u32 <- function() { b <- as.integer(take(4L)); sum(b * c(1, 2^8, 2^16, 2^24)) }
  i32 <- function() { x <- u32(); if (x >= 2^31) x - 2^32 else x }

  if (!identical(take(4L), charToRaw("QPP1"))) stop("invalid QPP magic", call. = FALSE)
  n_time <- u16()
  n_parcel <- u16()
  if (n_time < 1 || n_parcel < 1) stop("impossible QPP dimensions", call. = FALSE)
  q <- matrix(0, nrow = n_time, ncol = n_parcel)
  for (k in seq_len(n_parcel)) {
    predictor_id <- u8()
    width <- u8()
    n_anchor <- u8()
    reserved <- u8()
    n_residual <- u16()
    payload_nbyte <- u32()
    if (predictor_id > 2L || reserved != 0L || width > 34L) stop("invalid QPP header", call. = FALSE)
    predictor <- c("P0", "P1", "P2")[[predictor_id + 1L]]
    expected_anchor <- c(P0 = 0L, P1 = 1L, P2 = 2L)[[predictor]]
    if (n_anchor != expected_anchor || n_residual != n_time - n_anchor ||
        payload_nbyte != ceiling(n_residual * width / 8)) {
      stop("inconsistent QPP record lengths", call. = FALSE)
    }
    anchors <- if (n_anchor) vapply(seq_len(n_anchor), function(i) i32(), numeric(1)) else numeric()
    unsigned <- .qpp_ref_unpack(take(payload_nbyte), n_residual, width)
    values <- .qpp_ref_inverse_residuals(anchors, .qpp_ref_unzigzag(unsigned), predictor)
    .qpp_ref_assert_whole(values, .qpp_ref_q_min, .qpp_ref_q_max, "decoded q")
    q[, k] <- values
  }
  if (pos != length(bytes) + 1L) stop("trailing bytes after QPP tile", call. = FALSE)
  q
}

.qpp_ref_hex <- function(bytes) paste(sprintf("%02x", as.integer(bytes)), collapse = "")

.qpp_ref_unhex <- function(text) {
  if (!is.character(text) || length(text) != 1L || nchar(text) %% 2L != 0L ||
      grepl("[^0-9A-Fa-f]", text)) {
    stop("hex input must contain complete hexadecimal bytes", call. = FALSE)
  }
  starts <- seq.int(1L, nchar(text), by = 2L)
  as.raw(strtoi(substring(text, starts, starts + 1L), base = 16L))
}

.qpp_ref_synthetic_ar <- function(n_time = 129L, n_parcel = 4L, rho = 0.8) {
  innovations <- outer(seq_len(n_time), seq_len(n_parcel), function(t, k) {
    sin(t * (k + 1) * 0.37) + cos(t * (k + 2) * 0.11) + ((t * k) %% 7 - 3) / 5
  })
  y <- matrix(0, n_time, n_parcel)
  y[1L, ] <- innovations[1L, ]
  if (n_time > 1L) for (t in 2:n_time) y[t, ] <- rho * y[t - 1L, ] + innovations[t, ]
  y
}

.qpp_ref_chunk_coverage <- function(index, time_index, parcel_index) {
  required <- c("time_start", "time_end", "parcel_start", "parcel_end", "byte_length")
  if (!is.data.frame(index) || !all(required %in% names(index))) {
    stop("index is missing required coverage columns", call. = FALSE)
  }
  hit <- index$time_start <= max(time_index) & index$time_end >= min(time_index) &
    index$parcel_start <= max(parcel_index) & index$parcel_end >= min(parcel_index)
  list(tile_count = sum(hit), byte_count = sum(index$byte_length[hit]), rows = which(hit))
}
