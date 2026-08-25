# QPP rate-speed-fidelity court.
#
# This script is intentionally outside the package API. It selects candidate
# settings on a deterministic calibration corpus, freezes that choice, and only
# then evaluates it on a disjoint corpus. Synthetic data make the court
# redistributable and reproducible; a production promotion still requires an
# external-data receipt following inst/extdata/qpp-calibration/README.md.

.qpp_bench_version <- "1.0.0"
.qpp_bench_index_record_bytes <- 67L

.qpp_bench_hrf <- function(length = 24L) {
  time <- 0:(as.integer(length) - 1L)
  value <- stats::dgamma(time, shape = 6, scale = 1) -
    0.35 * stats::dgamma(time, shape = 12, scale = 1)
  value / max(abs(value))
}

.qpp_bench_convolve <- function(onsets, amplitudes, n_time, hrf) {
  out <- numeric(n_time)
  for (i in seq_along(onsets)) {
    index <- onsets[[i]] + seq_along(hrf) - 1L
    keep <- index <= n_time
    out[index[keep]] <- out[index[keep]] + amplitudes[[i]] * hrf[keep]
  }
  out
}

qpp_benchmark_simulate_run <- function(seed, split, run_index,
                                       n_time = 192L, n_parcel = 64L) {
  stopifnot(split %in% c("calibration", "evaluation"))
  n_time <- as.integer(n_time)
  n_parcel <- as.integer(n_parcel)
  if (n_time < 64L || n_parcel < 4L) {
    stop("benchmark runs require at least 64 TRs and 4 parcels", call. = FALSE)
  }
  set.seed(as.integer(seed))
  time <- seq_len(n_time)
  drift <- as.numeric(scale(seq(-1, 1, length.out = n_time)))
  block <- rep(c(-1, 1), each = 12L, length.out = n_time)
  block <- .qpp_bench_convolve(which(diff(c(0, block)) != 0),
                               rep(1, sum(diff(c(0, block)) != 0)),
                               n_time, .qpp_bench_hrf(18L))
  block <- as.numeric(scale(block))
  design <- cbind(intercept = 1, block = block, drift = drift)

  trial_onsets <- unique(as.integer(round(seq(10, n_time - 28, length.out = 12L))))
  trial_labels <- rep(c(-1, 1), length.out = length(trial_onsets))
  hrf <- .qpp_bench_hrf()
  trial_design <- vapply(trial_onsets, function(onset) {
    .qpp_bench_convolve(onset, 1, n_time, hrf)
  }, numeric(n_time))
  colnames(trial_design) <- sprintf("trial-%02d", seq_len(ncol(trial_design)))

  task_pattern <- as.numeric(scale(seq(-1, 1, length.out = n_parcel)))
  trial_amplitudes <- outer(trial_labels * 0.16, task_pattern)
  trial_amplitudes <- trial_amplitudes +
    matrix(stats::rnorm(length(trial_amplitudes), sd = 0.015),
           nrow = nrow(trial_amplitudes))
  beta <- rbind(
    intercept = stats::rnorm(n_parcel, 0, 0.03),
    block = stats::rnorm(n_parcel, 0.12, 0.04) * sign(task_pattern),
    drift = stats::rnorm(n_parcel, 0, 0.015)
  )

  factor_count <- 4L
  factors <- matrix(0, n_time, factor_count)
  innovations <- matrix(stats::rnorm(n_time * factor_count, sd = 0.04),
                        nrow = n_time)
  for (i in 2:n_time) factors[i, ] <- 0.65 * factors[i - 1L, ] + innovations[i, ]
  loadings <- matrix(stats::rnorm(n_parcel * factor_count, sd = 0.2),
                     nrow = factor_count)
  groups <- rep(seq_len(factor_count), length.out = n_parcel)
  loadings[cbind(groups, seq_len(n_parcel))] <-
    loadings[cbind(groups, seq_len(n_parcel))] + 0.75

  innovation_sd <- seq(0.085, 0.115, length.out = n_parcel)
  noise_innov <- sweep(matrix(stats::rnorm(n_time * n_parcel), nrow = n_time),
                       2L, innovation_sd, "*")
  noise <- matrix(0, n_time, n_parcel)
  noise[1L, ] <- noise_innov[1L, ] / sqrt(1 - 0.35^2)
  for (i in 2:n_time) noise[i, ] <- 0.35 * noise[i - 1L, ] + noise_innov[i, ]

  series <- design %*% beta + trial_design %*% trial_amplitudes +
    factors %*% loadings + noise
  colnames(series) <- sprintf("parcel-%04d", seq_len(n_parcel))
  list(
    id = sprintf("synthetic-%s-run-%02d", split, run_index),
    split = split,
    seed = as.integer(seed),
    series = unname(series),
    design = design,
    contrast = c(intercept = 0, block = 1, drift = 0),
    trial_design = trial_design,
    trial_labels = as.integer(trial_labels),
    simulation_innovation_sd = innovation_sd,
    units = "percent-signal-change"
  )
}

qpp_benchmark_corpus <- function(split = c("calibration", "evaluation"),
                                 n_runs = 2L, n_time = 192L,
                                 n_parcel = 64L) {
  split <- match.arg(split)
  base_seed <- if (identical(split, "calibration")) 202608250L else 202608750L
  lapply(seq_len(as.integer(n_runs)), function(i) {
    qpp_benchmark_simulate_run(
      seed = base_seed + i, split = split, run_index = i,
      n_time = n_time, n_parcel = n_parcel
    )
  })
}

qpp_benchmark_edge_cases <- function(n_time = 129L) {
  n_time <- as.integer(n_time)
  list(
    constant = matrix(7, n_time, 2L),
    ramp = cbind(seq_len(n_time), -seq_len(n_time)),
    impulse = cbind(replace(integer(n_time), n_time %/% 2L, 1000L),
                    replace(integer(n_time), n_time, -1000L)),
    alternating = cbind(rep(c(-1L, 1L), length.out = n_time),
                        rep(c(-32768L, 32767L), length.out = n_time))
  )
}

.qpp_bench_fit <- function(design, series) {
  qr_design <- qr(design)
  beta <- qr.coef(qr_design, series)
  fitted <- design %*% beta
  residual <- series - fitted
  list(
    beta = beta,
    residual = residual,
    df = nrow(design) - qr_design$rank,
    xtx_inverse = chol2inv(qr.R(qr_design))
  )
}

.qpp_bench_trial_betas <- function(run, series) {
  design <- cbind(run$design, run$trial_design)
  fit <- .qpp_bench_fit(design, series)
  start <- ncol(run$design) + 1L
  fit$beta[start:nrow(fit$beta), , drop = FALSE]
}

.qpp_bench_classifier <- function(corpus) {
  betas <- do.call(rbind, lapply(corpus, function(run) {
    .qpp_bench_trial_betas(run, run$series)
  }))
  labels <- unlist(lapply(corpus, `[[`, "trial_labels"), use.names = FALSE)
  positive <- colMeans(betas[labels == 1L, , drop = FALSE])
  negative <- colMeans(betas[labels == -1L, , drop = FALSE])
  pooled <- apply(betas, 2L, stats::var)
  weight <- (positive - negative) / pmax(pooled, 1e-8)
  # Decision-value agreement needs a fixed scale. Unit-normalizing the
  # calibration-fitted discriminant makes its tolerance interpretable in the
  # trial-beta units instead of depending on parcel count or covariance scale.
  weight <- weight / pmax(sqrt(sum(weight^2)), .Machine$double.eps)
  intercept <- -sum(weight * (positive + negative) / 2)
  list(weight = weight, intercept = intercept,
       fit_split = "calibration", n_trials = length(labels))
}

.qpp_bench_acf <- function(series, max_lag = 8L) {
  vapply(seq_len(max_lag), function(lag) {
    mean(vapply(seq_len(ncol(series)), function(k) {
      stats::cor(series[seq_len(nrow(series) - lag), k],
                 series[(lag + 1L):nrow(series), k])
    }, numeric(1)))
  }, numeric(1))
}

.qpp_bench_spectrum <- function(series) {
  centered <- sweep(series, 2L, colMeans(series), "-")
  power <- apply(centered, 2L, function(x) Mod(stats::fft(x))^2)
  power[seq_len(floor(nrow(series) / 2L) + 1L), , drop = FALSE]
}

.qpp_bench_connectivity <- function(series) {
  value <- stats::cor(series)
  value[upper.tri(value)]
}

.qpp_bench_scientific_metrics <- function(run, reconstructed, classifier) {
  original <- run$series
  error <- reconstructed - original
  original_fit <- .qpp_bench_fit(run$design, original)
  reconstructed_fit <- .qpp_bench_fit(run$design, reconstructed)
  contrast <- as.numeric(run$contrast)
  original_contrast <- as.numeric(crossprod(contrast, original_fit$beta))
  reconstructed_contrast <- as.numeric(crossprod(contrast, reconstructed_fit$beta))
  contrast_variance <- as.numeric(t(contrast) %*% original_fit$xtx_inverse %*% contrast)
  original_sigma <- colSums(original_fit$residual^2) / original_fit$df
  reconstructed_sigma <- colSums(reconstructed_fit$residual^2) / reconstructed_fit$df
  original_t <- original_contrast / sqrt(contrast_variance * original_sigma)
  reconstructed_t <- reconstructed_contrast /
    sqrt(contrast_variance * reconstructed_sigma)
  original_trial <- .qpp_bench_trial_betas(run, original)
  reconstructed_trial <- .qpp_bench_trial_betas(run, reconstructed)
  original_decision <- as.numeric(original_trial %*% classifier$weight +
                                    classifier$intercept)
  reconstructed_decision <- as.numeric(reconstructed_trial %*% classifier$weight +
                                         classifier$intercept)
  original_spectrum <- .qpp_bench_spectrum(original)
  reconstructed_spectrum <- .qpp_bench_spectrum(reconstructed)
  original_connectivity <- .qpp_bench_connectivity(original)
  reconstructed_connectivity <- .qpp_bench_connectivity(reconstructed)
  list(
    maximum_absolute_error = max(abs(error)),
    RMSE = sqrt(mean(error^2)),
    spectrum_relative_RMSE = sqrt(
      sum((reconstructed_spectrum - original_spectrum)^2) /
        pmax(sum(original_spectrum^2), .Machine$double.eps)
    ),
    autocorrelation_max_absolute_error = max(abs(
      .qpp_bench_acf(reconstructed) - .qpp_bench_acf(original)
    )),
    glm_beta_relative_RMSE = sqrt(
      sum((reconstructed_fit$beta - original_fit$beta)^2) /
        pmax(sum(original_fit$beta^2), .Machine$double.eps)
    ),
    glm_contrast_max_absolute_error = max(abs(
      reconstructed_contrast - original_contrast
    )),
    glm_t_max_absolute_error = max(abs(reconstructed_t - original_t)),
    connectivity_edge_RMSE = sqrt(mean(
      (reconstructed_connectivity - original_connectivity)^2
    )),
    connectivity_edge_max_absolute_error = max(abs(
      reconstructed_connectivity - original_connectivity
    )),
    single_trial_beta_RMSE = sqrt(mean(
      (reconstructed_trial - original_trial)^2
    )),
    mvpa_prediction_agreement = mean(
      sign(reconstructed_decision) == sign(original_decision)
    ),
    mvpa_decision_value_RMSE = sqrt(mean(
      (reconstructed_decision - original_decision)^2
    ))
  )
}

.qpp_bench_noise_sd <- function(corpus) {
  estimates <- lapply(corpus, function(run) {
    design <- cbind(run$design, run$trial_design)
    residual <- .qpp_bench_fit(design, run$series)$residual
    previous <- residual[-nrow(residual), , drop = FALSE]
    current <- residual[-1L, , drop = FALSE]
    rho <- colSums(previous * current) /
      pmax(colSums(previous^2), .Machine$double.eps)
    rho <- pmin(0.95, pmax(-0.95, rho))
    innovation <- current - sweep(previous, 2L, rho, "*")
    sqrt(colMeans(innovation^2))
  })
  apply(do.call(rbind, estimates), 2L, stats::median)
}

.qpp_bench_delta <- function(calibration, alpha) {
  sqrt(12 * as.numeric(alpha)) * .qpp_bench_noise_sd(calibration)
}

.qpp_bench_encode_once <- function(series, delta, time_chunk, parcel_chunk) {
  time_chunk <- as.integer(time_chunk)
  parcel_chunk <- as.integer(parcel_chunk)
  payloads <- list()
  decoded <- matrix(NA_real_, nrow(series), ncol(series))
  row <- 0L
  for (time_start in seq.int(1L, nrow(series), by = time_chunk)) {
    time_end <- min(nrow(series), time_start + time_chunk - 1L)
    for (parcel_start in seq.int(1L, ncol(series), by = parcel_chunk)) {
      parcel_end <- min(ncol(series), parcel_start + parcel_chunk - 1L)
      parcel_index <- parcel_start:parcel_end
      q <- round(sweep(series[time_start:time_end, parcel_index, drop = FALSE],
                       2L, delta[parcel_index], "/"))
      bytes <- neuroarchive:::qpp_encode_tile(q)
      roundtrip <- neuroarchive:::qpp_decode_tile(bytes)
      if (!identical(unname(roundtrip), unname(q * 1))) {
        stop("QPP q-roundtrip failed inside benchmark", call. = FALSE)
      }
      decoded[time_start:time_end, parcel_index] <-
        sweep(roundtrip, 2L, delta[parcel_index], "*")
      row <- row + 1L
      payloads[[row]] <- bytes
    }
  }
  list(
    payloads = payloads,
    reconstructed = decoded,
    payload_bytes = sum(lengths(payloads)),
    tile_count = length(payloads)
  )
}

.qpp_bench_decode_once <- function(encoded, delta, n_time, n_parcel,
                                   time_chunk, parcel_chunk) {
  decoded <- matrix(NA_real_, n_time, n_parcel)
  row <- 0L
  for (time_start in seq.int(1L, n_time, by = time_chunk)) {
    time_end <- min(n_time, time_start + time_chunk - 1L)
    for (parcel_start in seq.int(1L, n_parcel, by = parcel_chunk)) {
      parcel_end <- min(n_parcel, parcel_start + parcel_chunk - 1L)
      parcel_index <- parcel_start:parcel_end
      row <- row + 1L
      q <- neuroarchive:::qpp_decode_tile(encoded$payloads[[row]])
      decoded[time_start:time_end, parcel_index] <-
        sweep(q, 2L, delta[parcel_index], "*")
    }
  }
  decoded
}

.qpp_bench_elapsed <- function(fun, repetitions) {
  elapsed <- numeric(repetitions)
  fun()
  for (i in seq_len(repetitions)) {
    elapsed[[i]] <- system.time(fun())[["elapsed"]]
  }
  list(
    median_seconds = stats::median(elapsed),
    min_seconds = min(elapsed),
    max_seconds = max(elapsed),
    repetitions = as.integer(repetitions)
  )
}

.qpp_bench_allocations <- function(fun) {
  path <- tempfile("qpp-rate-rprofmem-")
  on.exit(unlink(path), add = TRUE)
  Rprofmem(path)
  on.exit(Rprofmem(NULL), add = TRUE)
  fun()
  Rprofmem(NULL)
  lines <- readLines(path, warn = FALSE)
  sizes <- suppressWarnings(as.numeric(sub(" .*", "", lines)))
  sizes <- sizes[is.finite(sizes)]
  list(
    method = "Rprofmem-visible R allocations",
    allocation_count = length(sizes),
    cumulative_allocation_bytes = sum(sizes),
    maximum_single_allocation_bytes = if (length(sizes)) max(sizes) else 0
  )
}

.qpp_bench_summarize_metrics <- function(metrics) {
  names <- names(metrics[[1L]])
  values <- lapply(names, function(name) {
    current <- vapply(metrics, `[[`, numeric(1), name)
    if (identical(name, "mvpa_prediction_agreement")) min(current) else max(current)
  })
  stats::setNames(values, names)
}

qpp_benchmark_candidate <- function(corpus, delta, alpha, time_chunk,
                                    parcel_chunk, classifier,
                                    repetitions = 3L) {
  encoded <- lapply(corpus, function(run) {
    .qpp_bench_encode_once(run$series, delta, time_chunk, parcel_chunk)
  })
  reconstructed <- Map(function(run, value) {
    .qpp_bench_decode_once(
      value, delta, nrow(run$series), ncol(run$series),
      time_chunk, parcel_chunk
    )
  }, corpus, encoded)
  metrics <- Map(function(run, value) {
    .qpp_bench_scientific_metrics(run, value, classifier)
  }, corpus, reconstructed)
  logical_values <- sum(vapply(corpus, function(x) length(x$series), numeric(1)))
  payload_bytes <- sum(vapply(encoded, `[[`, numeric(1), "payload_bytes"))
  tile_count <- sum(vapply(encoded, `[[`, numeric(1), "tile_count"))
  index_bytes <- tile_count * .qpp_bench_index_record_bytes
  step_bytes <- length(delta) * 8 * length(corpus)
  settings <- list(
    alpha = alpha, time_chunk = as.integer(time_chunk),
    parcel_chunk = as.integer(parcel_chunk), outer_codec = "raw",
    predictor_cost = "anchors32-plus-residual-width-v1",
    predictor_tie_order = c("P0", "P1", "P2")
  )
  metadata_bytes <- nchar(as.character(jsonlite::toJSON(
    settings, auto_unbox = TRUE, null = "null", digits = NA
  )), type = "bytes") * length(corpus)
  total_transform_bytes <- payload_bytes + index_bytes + step_bytes + metadata_bytes
  encode_timing <- .qpp_bench_elapsed(function() {
    invisible(lapply(corpus, function(run) {
      .qpp_bench_encode_once(run$series, delta, time_chunk, parcel_chunk)
    }))
  }, repetitions)
  decode_timing <- .qpp_bench_elapsed(function() {
    invisible(Map(function(run, value) {
      .qpp_bench_decode_once(value, delta, nrow(run$series), ncol(run$series),
                             time_chunk, parcel_chunk)
    }, corpus, encoded))
  }, repetitions)
  first <- encoded[[1L]]
  sample_rows <- unique(as.integer(round(seq(
    1, length(first$payloads), length.out = min(8L, length(first$payloads))
  ))))
  random_timing <- .qpp_bench_elapsed(function() {
    invisible(lapply(first$payloads[sample_rows], neuroarchive:::qpp_decode_tile))
  }, repetitions)
  gzip_bytes <- sum(vapply(encoded, function(value) {
    length(memCompress(do.call(c, value$payloads), type = "gzip"))
  }, numeric(1)))
  concatenated <- lapply(encoded, function(value) do.call(c, value$payloads))
  gzip_payloads <- lapply(concatenated, memCompress, type = "gzip")
  gzip_encode_timing <- .qpp_bench_elapsed(function() {
    invisible(lapply(concatenated, memCompress, type = "gzip"))
  }, repetitions)
  gzip_decode_timing <- .qpp_bench_elapsed(function() {
    invisible(lapply(gzip_payloads, memDecompress, type = "gzip"))
  }, repetitions)
  parcel_float32_bytes <- logical_values * 4
  metric_summary <- .qpp_bench_summarize_metrics(metrics)
  list(
    candidate = settings,
    calibration_delta = as.numeric(delta),
    observed_quantization_variance_fraction = mean(
      (delta^2 / 12) / pmax(.qpp_bench_noise_sd(corpus)^2, .Machine$double.eps)
    ),
    rate = list(
      byte_boundary = paste(
        "measured QPP payload plus exact v1 index records, float64 steps,",
        "and serialized candidate settings; excludes HDF5 container overhead"
      ),
      logical_values = logical_values,
      payload_bytes = payload_bytes,
      index_bytes = index_bytes,
      step_bytes = step_bytes,
      metadata_bytes = metadata_bytes,
      total_transform_bytes = total_transform_bytes,
      effective_payload_bits_per_parcel_TR = payload_bytes * 8 / logical_values,
      effective_total_bits_per_parcel_TR = total_transform_bytes * 8 / logical_values,
      parcel_float32_bytes = parcel_float32_bytes,
      measured_parcel_float32_compression_ratio =
        parcel_float32_bytes / total_transform_bytes,
      gzip_exploratory_payload_bytes = gzip_bytes,
      gzip_payload_size_fraction = gzip_bytes / payload_bytes,
      tile_count = tile_count
    ),
    runtime = list(
      encode = c(encode_timing, list(
        million_values_per_second = logical_values /
          pmax(encode_timing$median_seconds, .Machine$double.eps) / 1e6
      )),
      sequential_decode = c(decode_timing, list(
        million_values_per_second = logical_values /
          pmax(decode_timing$median_seconds, .Machine$double.eps) / 1e6
      )),
      random_tile_decode = c(random_timing, list(
        tiles = length(sample_rows),
        touched_payload_bytes = sum(lengths(first$payloads[sample_rows]))
      )),
      exploratory_gzip_outer = list(
        encode = gzip_encode_timing, decode = gzip_decode_timing,
        note = "not a v1 decoder-supported outer codec"
      ),
      allocation_profile = .qpp_bench_allocations(function() {
        .qpp_bench_decode_once(first, delta, nrow(corpus[[1L]]$series),
                               ncol(corpus[[1L]]$series), time_chunk,
                               parcel_chunk)
      }),
      peak_memory = list(
        method = "deterministic working-set upper bound",
        bytes = max(lengths(first$payloads)) +
          min(time_chunk, nrow(corpus[[1L]]$series)) *
          min(parcel_chunk, ncol(corpus[[1L]]$series)) * 16
      )
    ),
    fidelity = metric_summary,
    per_run_fidelity = stats::setNames(metrics, vapply(corpus, `[[`, "", "id"))
  )
}

qpp_benchmark_gate_definition <- function() {
  list(
    version = "1.0.0",
    estimand = paste(
      "Agreement between original and independently pointwise-quantized",
      "parcel-by-TR series; predictor coding must be exactly q-reversible"
    ),
    aggregation = paste(
      "Worst run for error metrics and minimum run for prediction agreement;",
      "candidate selection uses calibration runs only"
    ),
    tolerances = list(
      observed_quantization_variance_fraction = 0.00105,
      RMSE_multiplier_over_uniform_oracle = 1.20,
      spectrum_relative_RMSE = 0.01,
      autocorrelation_max_absolute_error = 0.005,
      glm_beta_relative_RMSE = 0.01,
      glm_contrast_max_absolute_error = 0.003,
      glm_t_max_absolute_error = 0.12,
      connectivity_edge_RMSE = 0.002,
      connectivity_edge_max_absolute_error = 0.008,
      single_trial_beta_RMSE = 0.004,
      mvpa_prediction_agreement = 1,
      mvpa_decision_value_RMSE = 0.01
    ),
    critical_checks = c(
      "signed-int32 q domain", "exact q decode", "delta-derived max-error bound",
      "no evaluation data in selection", "installed-package validation"
    ),
    uncertainty = paste(
      "Runtime reports median and min/max over repetitions; scientific gates",
      "retain the worst individual run rather than a standard-error summary"
    )
  )
}

qpp_benchmark_apply_gate <- function(result, gate = qpp_benchmark_gate_definition()) {
  tolerance <- gate$tolerances
  oracle_rmse <- sqrt(mean(result$calibration_delta^2 / 12))
  checks <- c(
    variance_fraction =
      result$observed_quantization_variance_fraction <=
        tolerance$observed_quantization_variance_fraction,
    max_error = result$fidelity$maximum_absolute_error <=
      max(result$calibration_delta) / 2 + 64 * .Machine$double.eps,
    RMSE = result$fidelity$RMSE <=
      tolerance$RMSE_multiplier_over_uniform_oracle * oracle_rmse,
    spectrum = result$fidelity$spectrum_relative_RMSE <=
      tolerance$spectrum_relative_RMSE,
    autocorrelation = result$fidelity$autocorrelation_max_absolute_error <=
      tolerance$autocorrelation_max_absolute_error,
    glm_beta = result$fidelity$glm_beta_relative_RMSE <=
      tolerance$glm_beta_relative_RMSE,
    glm_contrast = result$fidelity$glm_contrast_max_absolute_error <=
      tolerance$glm_contrast_max_absolute_error,
    glm_t = result$fidelity$glm_t_max_absolute_error <=
      tolerance$glm_t_max_absolute_error,
    connectivity_RMSE = result$fidelity$connectivity_edge_RMSE <=
      tolerance$connectivity_edge_RMSE,
    connectivity_max = result$fidelity$connectivity_edge_max_absolute_error <=
      tolerance$connectivity_edge_max_absolute_error,
    single_trial = result$fidelity$single_trial_beta_RMSE <=
      tolerance$single_trial_beta_RMSE,
    mvpa_prediction = result$fidelity$mvpa_prediction_agreement >=
      tolerance$mvpa_prediction_agreement,
    mvpa_decision = result$fidelity$mvpa_decision_value_RMSE <=
      tolerance$mvpa_decision_value_RMSE
  )
  list(pass = all(checks), checks = as.list(checks))
}

.qpp_bench_pack_unsigned <- function(values, bits) {
  values <- as.integer(values)
  bits <- as.integer(bits)
  out <- raw(ceiling(length(values) * bits / 8))
  bit_position <- 0L
  for (value in values) {
    remaining <- bits
    current <- value
    while (remaining > 0L) {
      byte <- bit_position %/% 8L + 1L
      offset <- bit_position %% 8L
      take <- min(remaining, 8L - offset)
      mask <- bitwShiftL(bitwShiftL(1L, take) - 1L, offset)
      part <- bitwShiftL(bitwAnd(current, bitwShiftL(1L, take) - 1L), offset)
      out[[byte]] <- as.raw(bitwOr(as.integer(out[[byte]]), bitwAnd(part, mask)))
      current <- bitwShiftR(current, take)
      remaining <- remaining - take
      bit_position <- bit_position + take
    }
  }
  out
}

.qpp_bench_write_numeric_raw <- function(value, size) {
  connection <- rawConnection(raw(), open = "wb")
  on.exit(close(connection), add = TRUE)
  writeBin(as.numeric(value), connection, size = size, endian = "little")
  rawConnectionValue(connection)
}

qpp_benchmark_baselines <- function(corpus, selected) {
  series <- lapply(corpus, `[[`, "series")
  logical_values <- sum(vapply(series, length, numeric(1)))
  float32 <- lapply(series, .qpp_bench_write_numeric_raw, size = 4L)
  gzip <- lapply(float32, memCompress, type = "gzip")
  quant_baseline <- function(bits, packed = FALSE) {
    bytes <- 0
    for (value in series) {
      q <- neuroarchive:::.quantize_global(as.numeric(value), bits,
                                           method = "range", center = TRUE)
      code <- as.integer(q$q)
      payload <- if (isTRUE(packed)) {
        .qpp_bench_pack_unsigned(code, bits)
      } else if (bits <= 8L) {
        as.raw(code)
      } else {
        .qpp_bench_write_numeric_raw(code, 2L)
      }
      bytes <- bytes + length(payload) + 8L
    }
    bytes
  }
  delta_rle_bytes <- 0
  for (value in series) {
    q <- round(sweep(value, 2L, selected$calibration_delta, "/"))
    delta <- apply(q, 2L, diff)
    if (is.null(dim(delta))) delta <- matrix(delta, ncol = 1L)
    encoded <- lapply(seq_len(ncol(delta)), function(k) rle(delta[, k]))
    delta_rle_bytes <- delta_rle_bytes + ncol(q) * 4L +
      sum(vapply(encoded, function(x) 8L * length(x$lengths), numeric(1)))
  }
  qpp_bytes <- selected$rate$total_transform_bytes
  values <- list(
    list(name = "uncompressed parcel float32", bytes = sum(lengths(float32)),
         boundary = "actual little-endian float32 serialization", loss = "float32 cast"),
    list(name = "parcel float32 plus gzip", bytes = sum(lengths(gzip)),
         boundary = "actual base gzip serialization of float32 bytes", loss = "float32 cast"),
    list(name = "current quant uint8", bytes = quant_baseline(8L, FALSE),
         boundary = "actual uint8 codes plus float64 global scale and offset",
         loss = "whole-corpus fitted range quantizer; not analysis-safe"),
    list(name = "N-bit 12", bytes = quant_baseline(12L, TRUE),
         boundary = "actual 12-bit packed codes plus float64 global scale and offset",
         loss = "whole-corpus fitted range quantizer; not analysis-safe"),
    list(name = "delta plus RLE", bytes = delta_rle_bytes,
         boundary = "actual signed-int32 anchor and RLE length/value accounting",
         loss = "inherits selected QPP pointwise quantization"),
    list(name = "QPP selected", bytes = qpp_bytes,
         boundary = selected$rate$byte_boundary,
         loss = "external-calibration pointwise quantization")
  )
  lapply(values, function(value) {
    value$logical_values <- logical_values
    value$effective_bits_per_parcel_TR <- value$bytes * 8 / logical_values
    value$compression_ratio_vs_parcel_float32 <- logical_values * 4 / value$bytes
    value$byte_kind <- "measured_serialized_or_exact_wire_components"
    value
  })
}

.qpp_bench_archive_measurement <- function(run, selected, calibration_checksum) {
  n_parcel <- ncol(run$series)
  array <- array(t(run$series), dim = c(n_parcel, 1L, 1L, nrow(run$series)))
  mask <- array(TRUE, dim = c(n_parcel, 1L, 1L))
  labels <- array(seq_len(n_parcel), dim = dim(mask))
  affine <- diag(4)
  grid_id <- neuroarchive:::.lna_grid_id(
    dim(mask), affine, "synthetic-index", "first-axis-fastest"
  )
  path <- tempfile(fileext = ".h5")
  on.exit(unlink(path), add = TRUE)
  neuroarchive::write_lna(
    array, path, transforms = c("parcel.mean", "qpp"), mask = mask,
    transform_params = list(
      "parcel.mean" = list(
        labels = labels, fit_scope = "fixed",
        source_id = "qpp-benchmark-identity-parcels-v1",
        signal_units = run$units, spatial_units = "synthetic-index",
        voxel_ordering = "first-axis-fastest", map_grid_id = grid_id
      ),
      qpp = list(
        delta = selected$calibration_delta,
        quantizer_policy = "external-calibration", units = run$units,
        calibration_id = "qpp-synthetic-calibration-v1",
        calibration_checksum = calibration_checksum,
        time_chunk = selected$candidate$time_chunk,
        parcel_chunk = selected$candidate$parcel_chunk
      )
    ),
    header = list(dims = c(n_parcel, 1L, 1L), spacing = c(1, 1, 1),
                  origin = c(0, 0, 0), transform = affine,
                  spatial_units = "synthetic-index",
                  voxel_ordering = "first-axis-fastest"), checksum = "none"
  )
  receipt <- neuroarchive::validate_lna(
    path, checksum = FALSE, qpp = "all", receipt = TRUE
  )
  list(
    boundary = "complete HDF5 LNA file with identity parcel map",
    file_bytes = unname(file.size(path)),
    qpp_validation_receipt_checksum = receipt$qpp[[1L]]$checksums$receipt,
    installed_decoder = receipt$qpp[[1L]]$transform,
    tiles_checked = receipt$qpp[[1L]]$tiles_checked,
    tile_count = receipt$qpp[[1L]]$tile_count
  )
}

.qpp_bench_hash_corpus <- function(corpus) {
  connection <- rawConnection(raw(), open = "wb")
  on.exit(close(connection), add = TRUE)
  for (run in corpus) {
    writeBin(charToRaw(run$id), connection)
    writeBin(as.numeric(run$series), connection, size = 8L, endian = "little")
  }
  paste0("sha256:", digest::digest(rawConnectionValue(connection),
                                    algo = "sha256", serialize = FALSE))
}

qpp_run_rate_distortion_court <- function(
    n_runs = 2L, n_time = 192L, n_parcel = 64L,
    alpha = c(0.00025, 0.001, 0.004),
    time_chunk = c(64L, 128L), parcel_chunk = c(64L, 128L),
    repetitions = 3L) {
  calibration <- qpp_benchmark_corpus(
    "calibration", n_runs, n_time, n_parcel
  )
  evaluation <- qpp_benchmark_corpus(
    "evaluation", n_runs, n_time, n_parcel
  )
  if (length(intersect(vapply(calibration, `[[`, "", "id"),
                       vapply(evaluation, `[[`, "", "id")))) {
    stop("calibration and evaluation corpus IDs overlap", call. = FALSE)
  }
  classifier <- .qpp_bench_classifier(calibration)
  grid <- expand.grid(
    alpha = alpha, time_chunk = as.integer(time_chunk),
    parcel_chunk = as.integer(parcel_chunk), KEEP.OUT.ATTRS = FALSE
  )
  calibration_results <- vector("list", nrow(grid))
  gates <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    delta <- .qpp_bench_delta(calibration, grid$alpha[[i]])
    calibration_results[[i]] <- qpp_benchmark_candidate(
      calibration, delta, grid$alpha[[i]], grid$time_chunk[[i]],
      grid$parcel_chunk[[i]], classifier, repetitions
    )
    gates[[i]] <- qpp_benchmark_apply_gate(calibration_results[[i]])
  }
  passing <- which(vapply(gates, `[[`, logical(1), "pass"))
  if (!length(passing)) stop("no QPP candidate passed calibration gates", call. = FALSE)
  compression <- vapply(calibration_results[passing], function(x) {
    x$rate$measured_parcel_float32_compression_ratio
  }, numeric(1))
  selected_index <- passing[[which.max(compression)]]
  selected_calibration <- calibration_results[[selected_index]]
  evaluation_result <- qpp_benchmark_candidate(
    evaluation, selected_calibration$calibration_delta,
    selected_calibration$candidate$alpha,
    selected_calibration$candidate$time_chunk,
    selected_calibration$candidate$parcel_chunk,
    classifier, repetitions
  )
  evaluation_gate <- qpp_benchmark_apply_gate(evaluation_result)
  calibration_checksum <- .qpp_bench_hash_corpus(calibration)
  evaluation_checksum <- .qpp_bench_hash_corpus(evaluation)
  edge_roundtrip <- vapply(qpp_benchmark_edge_cases(), function(q) {
    identical(unname(neuroarchive:::qpp_decode_tile(
      neuroarchive:::qpp_encode_tile(q)
    )), unname(q * 1))
  }, logical(1))
  candidate_summary <- lapply(seq_along(calibration_results), function(i) {
    current <- calibration_results[[i]]
    list(
      candidate = current$candidate,
      pass = gates[[i]]$pass,
      checks = gates[[i]]$checks,
      total_transform_bytes = current$rate$total_transform_bytes,
      effective_total_bits_per_parcel_TR =
        current$rate$effective_total_bits_per_parcel_TR,
      parcel_float32_compression_ratio =
        current$rate$measured_parcel_float32_compression_ratio,
      encode_million_values_per_second =
        current$runtime$encode$million_values_per_second,
      decode_million_values_per_second =
        current$runtime$sequential_decode$million_values_per_second,
      fidelity = current$fidelity
    )
  })
  manifest <- list(
    manifest_version = "1.0.0",
    corpus_kind = "redistributable deterministic synthetic parcellated fMRI",
    calibration = list(
      ids = vapply(calibration, `[[`, "", "id"),
      checksum = calibration_checksum,
      selection_role = "candidate fitting and selection"
    ),
    evaluation = list(
      ids = vapply(evaluation, `[[`, "", "id"),
      checksum = evaluation_checksum,
      selection_role = "untouched until one candidate was frozen"
    ),
    disjoint = TRUE,
    dimensions = list(runs_per_split = as.integer(n_runs),
                      time = as.integer(n_time), parcels = as.integer(n_parcel)),
    units = "percent-signal-change",
    restricted_real_data_in_repository = FALSE,
    external_data_promotion_required = TRUE
  )
  list(
    receipt_version = .qpp_bench_version,
    profile_status = "experimental",
    stable_profile_claim = FALSE,
    reason = paste(
      "Synthetic calibration/evaluation court passed or failed as reported;",
      "representative external fMRI and clean package/release gates remain required"
    ),
    selection_rule = paste(
      "Among calibration-pass candidates, maximize measured compression versus",
      "parcel float32; evaluation metrics are not consulted during selection"
    ),
    selected_candidate_index = selected_index,
    selected_candidate = selected_calibration$candidate,
    outer_codec_decision = list(
      compared = c("raw", "exploratory base gzip"),
      selected = "raw",
      gzip_payload_size_fraction =
        selected_calibration$rate$gzip_payload_size_fraction,
      reason = paste(
        "QPP v1 has an exact raw outer-codec allow-list and bounded random",
        "tile reads; gzip is size/speed evidence only and is not silently",
        "promoted into the wire format"
      )
    ),
    frozen_delta = selected_calibration$calibration_delta,
    gate_definition = qpp_benchmark_gate_definition(),
    calibration_gate = gates[[selected_index]],
    evaluation_gate = evaluation_gate,
    calibration_candidates = candidate_summary,
    selected_calibration = selected_calibration,
    untouched_evaluation = evaluation_result,
    baselines = qpp_benchmark_baselines(evaluation, evaluation_result),
    archive_measurement = .qpp_bench_archive_measurement(
      evaluation[[1L]], evaluation_result, calibration_checksum
    ),
    edge_case_q_roundtrip = as.list(edge_roundtrip),
    corpus_manifest = manifest,
    accounting = list(
      measured = c(
        "serialized baseline bytes", "QPP tile payload bytes",
        "QPP index/step/settings bytes", "complete selected LNA file bytes",
        "elapsed runtime", "touched random-access payload bytes"
      ),
      derived_identity = paste(
        "voxel-float32 ratio = (V/K) times parcel-float32 ratio; no voxel ratio",
        "is reported without an explicit measured V/K"
      )
    ),
    execution = list(
      repetitions = as.integer(repetitions),
      uncertainty = qpp_benchmark_gate_definition()$uncertainty,
      R = R.version.string,
      neuroarchive = as.character(utils::packageVersion("neuroarchive")),
      platform = R.version$platform,
      system = as.list(Sys.info())
    )
  )
}

qpp_write_rate_distortion_receipt <- function(receipt, path) {
  json <- jsonlite::toJSON(
    receipt, auto_unbox = TRUE, null = "null", digits = NA, pretty = TRUE
  )
  writeLines(as.character(json), path, useBytes = TRUE)
  invisible(path)
}
