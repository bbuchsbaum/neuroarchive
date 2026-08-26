# Sparse/bit-plane benchmark for the partitioned integer lifting prototype.

.qpp_wavelet_source <- function(path) {
  if (!file.exists(path)) stop("missing benchmark dependency: ", path, call. = FALSE)
  source(path, local = parent.frame())
}

if (!exists("lift_forward_unit", mode = "function")) {
  .qpp_wavelet_source(file.path("inst", "prototypes",
                                 "temporal-lift53-reference.R"))
}
if (!exists("qpp_benchmark_corpus", mode = "function")) {
  .qpp_wavelet_source(file.path("inst", "benchmarks",
                                 "qpp-rate-distortion.R"))
}

.qpp_wavelet_zigzag <- function(value) {
  .lift_check_integer(value, "wavelet coefficient")
  ifelse(value >= 0, 2 * value, -2 * value - 1)
}

.qpp_wavelet_band_bytes <- function(values, step, level, kind) {
  values <- as.numeric(values)
  nonzero <- values != 0
  padded <- c(as.raw(nonzero), raw((-length(nonzero)) %% 8L))
  bitmap <- if (length(padded)) packBits(padded, type = "raw") else raw()
  unsigned <- .qpp_wavelet_zigzag(values[nonzero])
  width <- if (length(unsigned) && max(unsigned) > 0) {
    ceiling(log2(max(unsigned) + 1))
  } else {
    0L
  }
  packed <- if (width > 0L) {
    .qpp_bench_pack_unsigned(unsigned, width)
  } else {
    raw()
  }
  header_bytes <- 24L
  list(
    kind = kind, level = as.integer(level), step = as.numeric(step),
    logical_coefficients = length(values), nonzero_coefficients = sum(nonzero),
    zero_fraction = mean(!nonzero), bit_width = as.integer(width),
    bitmap_bytes = length(bitmap), packed_value_bytes = length(packed),
    header_bytes = header_bytes,
    total_bytes = length(bitmap) + length(packed) + header_bytes
  )
}

.qpp_wavelet_storage <- function(quantized) {
  transform <- quantized$transform
  payloads <- list()
  for (i in seq_len(nrow(transform$bands))) {
    band <- transform$bands[i, , drop = FALSE]
    index <- band$start:band$end
    step <- if (identical(band$kind, "coarse")) {
      1
    } else {
      quantized$detail_steps[[band$level]]
    }
    payloads[[i]] <- .qpp_wavelet_band_bytes(
      quantized$quantized_coefficients[index, , drop = FALSE],
      step, band$level, band$kind
    )
  }
  names(payloads) <- ifelse(
    vapply(payloads, `[[`, "", "kind") == "coarse", "coarse",
    paste0("detail-", vapply(payloads, `[[`, integer(1), "level"))
  )
  cumulative <- numeric(length(payloads))
  cumulative[[1L]] <- payloads[[1L]]$total_bytes
  if (length(payloads) > 1L) {
    for (i in 2:length(payloads)) {
      cumulative[[i]] <- cumulative[[i - 1L]] + payloads[[i]]$total_bytes
    }
  }
  list(
    encoding = paste(
      "one actual presence bitmap and one LSB-first packed unsigned zig-zag",
      "value stream per named band; 24-byte fixed band header"
    ),
    bands = payloads,
    progressive_stages = stats::setNames(as.list(cumulative), names(payloads)),
    total_bytes = sum(vapply(payloads, `[[`, numeric(1), "total_bytes")),
    nonzero_coefficients = sum(vapply(
      payloads, `[[`, numeric(1), "nonzero_coefficients"
    )),
    logical_coefficients = sum(vapply(
      payloads, `[[`, numeric(1), "logical_coefficients"
    ))
  )
}

qpp_wavelet_candidate_definitions <- function(levels = 3L) {
  levels <- as.integer(levels)
  exact <- rep(1, levels)
  names(exact) <- as.character(seq_len(levels))
  gentle <- pmax(1, 2^(pmax(0, 2L - seq_len(levels))))
  names(gentle) <- names(exact)
  compact <- pmax(1, 2^(pmax(0, 3L - seq_len(levels))))
  names(compact) <- names(exact)
  unlist(lapply(c("cdf53", "haar"), function(wavelet) {
    list(
      list(id = paste0(wavelet, "-lossless"), wavelet = wavelet,
           levels = levels, detail_steps = exact, dead_zone = 0,
           drop_finest = 0L),
      list(id = paste0(wavelet, "-gentle"), wavelet = wavelet,
           levels = levels, detail_steps = gentle, dead_zone = 0.5,
           drop_finest = 0L),
      list(id = paste0(wavelet, "-drop-finest"), wavelet = wavelet,
           levels = levels, detail_steps = exact, dead_zone = 0,
           drop_finest = 1L),
      list(id = paste0(wavelet, "-compact"), wavelet = wavelet,
           levels = levels, detail_steps = compact, dead_zone = 0.5,
           drop_finest = 1L)
    )
  }), recursive = FALSE)
}

qpp_wavelet_gate_definition <- function() {
  list(
    version = "1.0.0",
    estimand = paste(
      "Worst-run scientific agreement after run-partitioned integer lifting,",
      "declared detail quantization/band dropping, and exact inverse"
    ),
    tolerances = list(
      maximum_absolute_error = 0.04,
      RMSE = 0.01,
      spectrum_relative_RMSE = 0.03,
      autocorrelation_max_absolute_error = 0.012,
      glm_beta_relative_RMSE = 0.025,
      glm_contrast_max_absolute_error = 0.008,
      glm_t_max_absolute_error = 0.25,
      connectivity_edge_RMSE = 0.005,
      connectivity_edge_max_absolute_error = 0.02,
      single_trial_beta_RMSE = 0.01,
      mvpa_prediction_agreement = 1,
      mvpa_decision_value_RMSE = 0.012
    ),
    admission_rate_fraction_of_qpp = 0.90,
    aggregation = "maximum per-run error; minimum prediction agreement",
    selection = paste(
      "select the smallest calibration-pass wavelet payload; freeze before",
      "evaluation; admit only if evaluation passes and bytes <= 90% of QPP"
    )
  )
}

qpp_wavelet_apply_gate <- function(result,
                                   gate = qpp_wavelet_gate_definition()) {
  tolerance <- gate$tolerances
  metrics <- result$fidelity
  checks <- vapply(names(tolerance), function(name) {
    if (identical(name, "mvpa_prediction_agreement")) {
      metrics[[name]] >= tolerance[[name]]
    } else {
      metrics[[name]] <= tolerance[[name]]
    }
  }, logical(1))
  list(pass = all(checks), checks = as.list(checks))
}

.qpp_wavelet_quantized_series <- function(run, delta) {
  q <- round(sweep(run$series, 2L, delta, "/"))
  q[q == 0] <- 0
  q
}

qpp_wavelet_benchmark_candidate <- function(corpus, delta, candidate,
                                             classifier, repetitions = 3L) {
  encoded <- lapply(corpus, function(run) {
    q <- .qpp_wavelet_quantized_series(run, delta)
    transformed <- lift_forward_unit(
      q, candidate$levels, candidate$wavelet, unit_id = run$id
    )
    quantized <- lift_quantize_bands(
      transformed, candidate$detail_steps, candidate$dead_zone,
      candidate$drop_finest
    )
    reconstructed_q <- lift_inverse_quantized(quantized)
    list(
      q = q, transformed = transformed, quantized = quantized,
      reconstructed_q = reconstructed_q,
      storage = .qpp_wavelet_storage(quantized)
    )
  })
  reconstructed <- Map(function(value) {
    sweep(value$reconstructed_q, 2L, delta, "*")
  }, encoded)
  metrics <- Map(function(run, value) {
    .qpp_bench_scientific_metrics(run, value, classifier)
  }, corpus, reconstructed)
  summary <- .qpp_bench_summarize_metrics(metrics)
  logical_values <- sum(vapply(corpus, function(run) length(run$series), numeric(1)))
  total_bytes <- sum(vapply(encoded, function(value) {
    value$storage$total_bytes
  }, numeric(1)))
  forward_timing <- .qpp_bench_elapsed(function() {
    invisible(lapply(corpus, function(run) {
      lift_forward_unit(
        .qpp_wavelet_quantized_series(run, delta), candidate$levels,
        candidate$wavelet, run$id
      )
    }))
  }, repetitions)
  inverse_timing <- .qpp_bench_elapsed(function() {
    invisible(lapply(encoded, function(value) {
      lift_inverse_quantized(value$quantized)
    }))
  }, repetitions)
  exact_q <- all(vapply(encoded, function(value) {
    identical(value$q, value$reconstructed_q)
  }, logical(1)))
  list(
    candidate = candidate,
    partition_contract = list(
      independence_unit = "run", cross_unit_operations = FALSE,
      lossy_support = "complete-run"
    ),
    full_coefficient_q_roundtrip_exact =
      candidate$drop_finest == 0L && all(candidate$detail_steps == 1) &&
        candidate$dead_zone == 0 && exact_q,
    rate = list(
      byte_boundary = paste(
        "measured sparse bitmap plus packed zig-zag band streams and fixed",
        "band headers; excludes future LNA container/index overhead"
      ),
      logical_values = logical_values, total_bytes = total_bytes,
      effective_bits_per_parcel_TR = total_bytes * 8 / logical_values,
      nonzero_fraction = sum(vapply(encoded, function(value) {
        value$storage$nonzero_coefficients
      }, numeric(1))) / logical_values,
      progressive_stages_first_run = encoded[[1L]]$storage$progressive_stages
    ),
    runtime = list(forward = forward_timing, inverse = inverse_timing),
    fidelity = summary,
    per_run_fidelity = stats::setNames(metrics, vapply(corpus, `[[`, "", "id"))
  )
}

qpp_run_wavelet_prototype_court <- function(
    n_runs = 2L, n_time = 192L, n_parcel = 128L, levels = 3L,
    alpha = 0.001, qpp_time_chunk = 128L, qpp_parcel_chunk = 128L,
    repetitions = 3L) {
  calibration <- qpp_benchmark_corpus("calibration", n_runs, n_time, n_parcel)
  evaluation <- qpp_benchmark_corpus("evaluation", n_runs, n_time, n_parcel)
  classifier <- .qpp_bench_classifier(calibration)
  delta <- .qpp_bench_delta(calibration, alpha)
  candidates <- qpp_wavelet_candidate_definitions(levels)
  calibration_results <- lapply(candidates, function(candidate) {
    qpp_wavelet_benchmark_candidate(
      calibration, delta, candidate, classifier, repetitions
    )
  })
  calibration_gates <- lapply(calibration_results, qpp_wavelet_apply_gate)
  qpp_calibration <- qpp_benchmark_candidate(
    calibration, delta, alpha, qpp_time_chunk, qpp_parcel_chunk,
    classifier, repetitions
  )
  passing <- which(vapply(calibration_gates, `[[`, logical(1), "pass"))
  if (!length(passing)) stop("no wavelet candidate passed calibration", call. = FALSE)
  candidate_bytes <- vapply(calibration_results[passing], function(value) {
    value$rate$total_bytes
  }, numeric(1))
  selected_index <- passing[[which.min(candidate_bytes)]]
  selected <- calibration_results[[selected_index]]
  evaluation_result <- qpp_wavelet_benchmark_candidate(
    evaluation, delta, candidates[[selected_index]], classifier, repetitions
  )
  evaluation_gate <- qpp_wavelet_apply_gate(evaluation_result)
  qpp_evaluation <- qpp_benchmark_candidate(
    evaluation, delta, alpha, qpp_time_chunk, qpp_parcel_chunk,
    classifier, repetitions
  )
  rate_fraction <- evaluation_result$rate$total_bytes /
    qpp_evaluation$rate$payload_bytes
  lossless <- Filter(function(value) {
    isTRUE(value$full_coefficient_q_roundtrip_exact)
  }, calibration_results)
  list(
    receipt_version = "1.0.0",
    profile_status = "prototype",
    selected_candidate_index = selected_index,
    selected_candidate = candidates[[selected_index]],
    selection_rule = qpp_wavelet_gate_definition()$selection,
    gate_definition = qpp_wavelet_gate_definition(),
    calibration_candidates = Map(function(result, gate) {
      list(candidate = result$candidate, gate = gate, rate = result$rate,
           fidelity = result$fidelity,
           full_coefficient_q_roundtrip_exact =
             result$full_coefficient_q_roundtrip_exact)
    }, calibration_results, calibration_gates),
    selected_calibration = selected,
    untouched_evaluation = evaluation_result,
    evaluation_gate = evaluation_gate,
    qpp_comparator = list(
      settings = qpp_evaluation$candidate,
      byte_boundary = "measured QPP tile payload only; excludes index, steps, and container",
      calibration_payload_bytes = qpp_calibration$rate$payload_bytes,
      evaluation_payload_bytes = qpp_evaluation$rate$payload_bytes,
      evaluation_effective_payload_bits_per_parcel_TR =
        qpp_evaluation$rate$effective_payload_bits_per_parcel_TR
    ),
    wavelet_rate_fraction_of_qpp = rate_fraction,
    lossless_exact_for_both_wavelets =
      length(lossless) == 2L &&
        identical(sort(vapply(lossless, function(value) {
          value$candidate$wavelet
        }, "")), c("cdf53", "haar")),
    admit_to_qpp12 = isTRUE(evaluation_gate$pass) &&
      rate_fraction <= qpp_wavelet_gate_definition()$admission_rate_fraction_of_qpp,
    admission_reason = if (isTRUE(evaluation_gate$pass) &&
                           rate_fraction <= 0.90) {
      "calibration-selected profile passed untouched evaluation and used <=90% of QPP bytes"
    } else {
      paste(
        "do not integrate: requires both untouched scientific pass and <=90%",
        "of QPP bytes at the declared measured boundaries"
      )
    },
    package_boundary = list(
      production_decoder_owner = "neuroarchive-core",
      fmrilatent_role = "fitting, calibration, and experimental transforms",
      fmrilatent_companion_mote = "bd-01M0WJX0QVS0ZE07EF3VCN3AG4",
      fmrilatent_source_edited = FALSE
    ),
    corpus_manifest = list(
      calibration_ids = vapply(calibration, `[[`, "", "id"),
      evaluation_ids = vapply(evaluation, `[[`, "", "id"),
      disjoint = !length(intersect(
        vapply(calibration, `[[`, "", "id"),
        vapply(evaluation, `[[`, "", "id")
      )),
      delta_fit_scope = "external-calibration"
    ),
    execution = list(
      repetitions = as.integer(repetitions), R = R.version.string,
      neuroarchive = as.character(utils::packageVersion("neuroarchive")),
      platform = R.version$platform, system = as.list(Sys.info())
    )
  )
}

qpp_write_wavelet_prototype_receipt <- function(receipt, path) {
  json <- jsonlite::toJSON(
    receipt, auto_unbox = TRUE, null = "null", digits = NA, pretty = TRUE
  )
  writeLines(as.character(json), path, useBytes = TRUE)
  invisible(path)
}
