# Provenance-asserted reference identity and spatial compatibility policy.

.lna_reference_schema_version <- "1.0.0"
.lna_reference_origins <- c(
  "bids-space-entity", "spatial-reference-sidecar", "user-receipt", "unknown"
)
.lna_reference_kinds <- c("template", "named-reference", "unknown-native")
.lna_reference_statuses <- c("asserted", "overridden", "unknown")

.reference_abort <- function(message, field) {
  abort_lna(
    message,
    field = field,
    .subclass = "lna_error_reference_assertion",
    location = paste0("lna_reference_assertion:", field)
  )
}

.reference_scalar_string <- function(x, field, allow_null = FALSE) {
  if (allow_null && is.null(x)) return(NULL)
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    .reference_abort(paste0(field, " must be one non-empty string"), field)
  }
  enc2utf8(x)
}

.reference_named_list <- function(x, field, allow_empty = TRUE) {
  if (!is.list(x) || is.data.frame(x)) {
    .reference_abort(paste0(field, " must be a named list"), field)
  }
  if (!length(x) && allow_empty) return(x)
  nms <- names(x)
  if (is.null(nms) || anyNA(nms) || any(!nzchar(nms)) || anyDuplicated(nms)) {
    .reference_abort(
      paste0(field, " must have unique, non-empty names"),
      field
    )
  }
  x
}

.reference_json_value <- function(x, field) {
  tryCatch(
    .source_record_json_value(x, field),
    lna_error_source_record = function(e) {
      .reference_abort(conditionMessage(e), field)
    }
  )
  x
}

.reference_receipt <- function(x, field = "receipt") {
  x <- .reference_named_list(x, field, allow_empty = FALSE)
  required <- c("receipt_id", "asserted_by", "asserted_at", "reason")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    .reference_abort(
      paste0(field, " is missing field(s): ", paste(missing, collapse = ", ")),
      field
    )
  }
  for (nm in required) {
    x[[nm]] <- .reference_scalar_string(x[[nm]], paste0(field, ".", nm))
  }
  .reference_json_value(x, field)
  x
}

.reference_history_entry <- function(x, index) {
  field <- paste0("history[[", index, "]]")
  x <- .reference_named_list(x, field, allow_empty = FALSE)
  required <- c(
    "schema_version", "kind", "identifier", "origin", "evidence", "status",
    "receipt"
  )
  unknown <- setdiff(names(x), required)
  missing <- setdiff(required, names(x))
  if (length(unknown) || length(missing)) {
    .reference_abort(
      paste0(
        field, " must contain exactly: ", paste(required, collapse = ", ")
      ),
      field
    )
  }
  .reference_assertion_from_raw(
    c(x, list(history = list())),
    record_id = NULL,
    history_entry = TRUE
  )[required]
}

.reference_snapshot <- function(x) {
  unclass(x)[c(
    "schema_version", "kind", "identifier", "origin", "evidence", "status",
    "receipt"
  )]
}

.reference_assertion_from_raw <- function(
    x,
    record_id = NULL,
    history_entry = FALSE) {
  x <- .reference_named_list(x, "reference_assertion", allow_empty = FALSE)
  required <- c(
    "schema_version", "kind", "identifier", "origin", "evidence", "status",
    "receipt", "history"
  )
  unknown <- setdiff(names(x), required)
  missing <- setdiff(required, names(x))
  if (length(unknown)) {
    .reference_abort(
      paste0(
        "reference_assertion has unknown field(s): ",
        paste(unknown, collapse = ", ")
      ),
      "reference_assertion"
    )
  }
  if (length(missing)) {
    .reference_abort(
      paste0(
        "reference_assertion is missing field(s): ",
        paste(missing, collapse = ", ")
      ),
      "reference_assertion"
    )
  }

  schema_version <- .reference_scalar_string(
    x$schema_version, "schema_version"
  )
  if (!identical(schema_version, .lna_reference_schema_version)) {
    .reference_abort(
      paste0("unsupported reference schema_version: ", schema_version),
      "schema_version"
    )
  }
  kind <- .reference_scalar_string(x$kind, "kind")
  identifier <- .reference_scalar_string(x$identifier, "identifier")
  origin <- .reference_scalar_string(x$origin, "origin")
  status <- .reference_scalar_string(x$status, "status")
  if (!kind %in% .lna_reference_kinds) {
    .reference_abort(
      paste0("kind must be one of: ", paste(.lna_reference_kinds, collapse = ", ")),
      "kind"
    )
  }
  if (!origin %in% .lna_reference_origins) {
    .reference_abort(
      paste0(
        "origin must be one of: ",
        paste(.lna_reference_origins, collapse = ", ")
      ),
      "origin"
    )
  }
  if (!status %in% .lna_reference_statuses) {
    .reference_abort(
      paste0(
        "status must be one of: ",
        paste(.lna_reference_statuses, collapse = ", ")
      ),
      "status"
    )
  }

  evidence <- x$evidence %||% list()
  .reference_named_list(evidence, "evidence")
  .reference_json_value(evidence, "evidence")
  history <- x$history %||% list()
  if (!is.list(history) || is.data.frame(history)) {
    .reference_abort("history must be a list", "history")
  }
  normalized_history <- if (length(history)) {
    lapply(seq_along(history), function(i) {
      .reference_history_entry(history[[i]], i)
    })
  } else {
    list()
  }
  receipt <- if (is.null(x$receipt)) NULL else .reference_receipt(x$receipt)

  if (identical(origin, "unknown")) {
    if (!identical(kind, "unknown-native") || !identical(status, "unknown")) {
      .reference_abort(
        "origin='unknown' requires kind='unknown-native' and status='unknown'",
        "origin"
      )
    }
    if (!history_entry) {
      record_id <- .reference_scalar_string(record_id, "record_id")
      expected <- paste0("unknown-native:", record_id)
      if (!identical(identifier, expected)) {
        .reference_abort(
          paste0("unknown-native identifier must equal ", expected),
          "identifier"
        )
      }
    } else if (!startsWith(identifier, "unknown-native:")) {
      .reference_abort(
        "unknown-native history identifiers must start with 'unknown-native:'",
        "identifier"
      )
    }
    if (!is.null(receipt) || length(normalized_history)) {
      .reference_abort(
        "an unknown-native assertion cannot carry a receipt or history",
        "receipt"
      )
    }
  } else {
    if (identical(kind, "unknown-native")) {
      .reference_abort(
        "known assertion origins cannot use kind='unknown-native'",
        "kind"
      )
    }
    if (identical(status, "unknown")) {
      .reference_abort(
        "known assertion origins cannot use status='unknown'",
        "status"
      )
    }
  }

  if (identical(origin, "bids-space-entity")) {
    if (!identical(evidence$space_entity, identifier)) {
      .reference_abort(
        "BIDS-space assertions require evidence.space_entity matching identifier",
        "evidence.space_entity"
      )
    }
  }
  if (identical(origin, "spatial-reference-sidecar")) {
    if (!identical(evidence$SpatialReference, identifier)) {
      .reference_abort(
        paste0(
          "sidecar assertions require evidence.SpatialReference matching ",
          "identifier"
        ),
        "evidence.SpatialReference"
      )
    }
  }
  if (identical(origin, "user-receipt")) {
    if (is.null(receipt)) {
      .reference_abort(
        "user reference overrides require an explicit receipt",
        "receipt"
      )
    }
    if (!identical(status, "overridden") ||
        (!history_entry && !length(normalized_history))) {
      .reference_abort(
        "user reference overrides require status='overridden' and prior history",
        "history"
      )
    }
  } else if (!is.null(receipt)) {
    .reference_abort(
      "receipt is reserved for origin='user-receipt'",
      "receipt"
    )
  }
  if (!identical(origin, "user-receipt") && identical(status, "overridden")) {
    .reference_abort(
      "status='overridden' requires origin='user-receipt'",
      "status"
    )
  }

  structure(
    list(
      schema_version = schema_version,
      kind = kind,
      identifier = identifier,
      origin = origin,
      evidence = evidence,
      status = status,
      receipt = receipt,
      history = normalized_history
    ),
    class = "lna_reference_assertion"
  )
}

#' Create a provenance-asserted reference identity
#'
#' Reference identity is an assertion about spatial meaning, not a hash of
#' geometry. Template names are therefore never inferred from dimensions or an
#' affine. Missing information becomes a per-record unknown-native assertion.
#'
#' @param record_id Record identifier that owns the assertion.
#' @param identifier Asserted reference identifier. Leave `NULL` only for an
#'   unknown-native assertion.
#' @param kind One of `"template"`, `"named-reference"`, or
#'   `"unknown-native"`.
#' @param origin One of `"bids-space-entity"`,
#'   `"spatial-reference-sidecar"`, `"user-receipt"`, or `"unknown"`.
#' @param evidence Named JSON-safe evidence. BIDS and sidecar origins require
#'   matching `space_entity` and `SpatialReference` evidence, respectively.
#' @param status Assertion status. Usually inferred from `origin`.
#' @param receipt Explicit user-override receipt.
#' @param history Prior assertion snapshots. User-facing overrides should use
#'   [override_lna_reference()] so history cannot be dropped accidentally.
#' @return An object of class `lna_reference_assertion`.
#' @export
lna_reference_assertion <- function(
    record_id,
    identifier = NULL,
    kind = NULL,
    origin = "unknown",
    evidence = list(),
    status = NULL,
    receipt = NULL,
    history = list()) {
  origin <- .reference_scalar_string(origin, "origin")
  if (identical(origin, "unknown")) {
    kind <- kind %||% "unknown-native"
    identifier <- identifier %||% paste0("unknown-native:", record_id)
    status <- status %||% "unknown"
  } else {
    kind <- kind %||% "named-reference"
    status <- status %||% if (identical(origin, "user-receipt")) {
      "overridden"
    } else {
      "asserted"
    }
  }
  .reference_assertion_from_raw(
    list(
      schema_version = .lna_reference_schema_version,
      kind = kind,
      identifier = identifier,
      origin = origin,
      evidence = evidence,
      status = status,
      receipt = receipt,
      history = history
    ),
    record_id = record_id
  )
}

.as_lna_reference_assertion <- function(x, record_id) {
  if (inherits(x, "lna_reference_assertion")) x <- unclass(x)
  .reference_assertion_from_raw(x, record_id = record_id)
}

#' Override a source record's reference assertion with a receipt
#'
#' The returned source record contains the prior assertion and all earlier
#' history as immutable snapshots. The input record is not modified.
#'
#' @param record An `lna_source_record`.
#' @param identifier New asserted reference identifier.
#' @param receipt Named list with `receipt_id`, `asserted_by`, `asserted_at`,
#'   and `reason`.
#' @param kind New reference kind, normally `"template"`.
#' @param evidence Additional JSON-safe user evidence.
#' @return A new `lna_source_record`.
#' @export
override_lna_reference <- function(
    record,
    identifier,
    receipt,
    kind = "template",
    evidence = list()) {
  record <- as_lna_source_record(record)
  previous <- record$reference_assertion
  history <- c(previous$history, list(.reference_snapshot(previous)))
  assertion <- lna_reference_assertion(
    record_id = record$record_id,
    identifier = identifier,
    kind = kind,
    origin = "user-receipt",
    evidence = evidence,
    status = "overridden",
    receipt = receipt,
    history = history
  )
  raw <- as.list(record)
  raw$reference_assertion <- unclass(assertion)
  as_lna_source_record(raw)
}

.spatial_compatibility_reason <- function(reference, grid, support, policy) {
  reasons <- character()
  if (!reference) reasons <- c(reasons, "reference-assertion-mismatch")
  if (!grid) reasons <- c(reasons, "computed-grid-mismatch")
  if (identical(policy, "same-support") && !support) {
    reasons <- c(reasons, "computed-support-missing-or-mismatch")
  }
  reasons
}

#' Compare asserted reference and computed spatial identities
#'
#' This receipt keeps scientific assertions separate from byte-derived
#' identities. Compatibility always requires matching reference assertions and
#' matching computed grids; `support_policy = "same-support"` additionally
#' requires two present, equal support identities.
#'
#' @param x,y LNA source records.
#' @param support_policy Either `"same-support"` (fail closed when support is
#'   absent) or `"grid-only"`.
#' @return An `lna_spatial_compatibility` receipt.
#' @export
lna_spatial_compatibility <- function(
    x,
    y,
    support_policy = c("same-support", "grid-only")) {
  support_policy <- match.arg(support_policy)
  x <- as_lna_source_record(x)
  y <- as_lna_source_record(y)

  x_ref <- x$reference_assertion
  y_ref <- y$reference_assertion
  both_unknown <- identical(x_ref$kind, "unknown-native") &&
    identical(y_ref$kind, "unknown-native")
  reference_match <- if (both_unknown) {
    identical(x$record_id, y$record_id) &&
      identical(x_ref$identifier, y_ref$identifier)
  } else {
    !identical(x_ref$kind, "unknown-native") &&
      !identical(y_ref$kind, "unknown-native") &&
      identical(x_ref$kind, y_ref$kind) &&
      identical(x_ref$identifier, y_ref$identifier)
  }

  grid_match <- identical(x$geometry$grid_id, y$geometry$grid_id)
  support_present <- !is.null(x$geometry$support_id) &&
    !is.null(y$geometry$support_id)
  support_match <- support_present &&
    identical(x$geometry$support_id, y$geometry$support_id)
  support_compatible <- if (identical(support_policy, "same-support")) {
    support_match
  } else {
    TRUE
  }
  compatible <- reference_match && grid_match && support_compatible

  structure(
    list(
      compatible = compatible,
      policy = list(
        reference = "matching-assertion-required",
        grid = "same-grid-required",
        support = support_policy
      ),
      reference = list(
        compatible = reference_match,
        identity_kind = "provenance-asserted",
        verification = "assertion-validated",
        left = x_ref$identifier,
        right = y_ref$identifier,
        left_origin = x_ref$origin,
        right_origin = y_ref$origin
      ),
      grid = list(
        compatible = grid_match,
        identity_kind = "computed-verified",
        verification = "recomputed-from-geometry",
        left = x$geometry$grid_id,
        right = y$geometry$grid_id
      ),
      support = list(
        compatible = support_compatible,
        identities_match = support_match,
        identity_kind = "computed-identifier",
        verification = if (support_present) {
          "identifier-only"
        } else {
          "missing"
        },
        left = x$geometry$support_id,
        right = y$geometry$support_id,
        required = identical(support_policy, "same-support")
      ),
      reasons = .spatial_compatibility_reason(
        reference_match, grid_match, support_compatible, support_policy
      )
    ),
    class = "lna_spatial_compatibility"
  )
}
