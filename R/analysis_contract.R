# Machine-checkable analysis and leakage contracts.

.analysis_fit_scopes <- c(
  "not-applicable", "unknown", "fixed", "external-calibration", "fold",
  "pointwise", "causal-prefix", "chunk", "trial-epoch", "run", "session",
  "subject"
)

.analysis_unit_scopes <- c(
  "pointwise", "causal-prefix", "chunk", "trial-epoch", "run", "session",
  "subject", "fold"
)

.analysis_lossy_supports <- c(.analysis_unit_scopes, "parcel-at-same-TR")
.analysis_coding_dependencies <- c("none", .analysis_unit_scopes)

.analysis_scalar_string <- function(x, field, allow_empty = FALSE) {
  ok <- is.character(x) && length(x) == 1L && !is.na(x)
  if (ok && !allow_empty) ok <- nzchar(x)
  if (!ok) {
    abort_lna(
      paste0(field, " must be one non-empty string"),
      .subclass = "lna_error_analysis_contract",
      location = paste0("analysis_contract:", field)
    )
  }
  x
}

.analysis_normalization_record <- function(x, canonicalize_checksums = FALSE) {
  if (is.null(x)) {
    return(list(
      method = "unspecified", fit_scope = "unknown", source_id = "unrecorded",
      calibration_identity = NULL, fit_id = NULL
    ))
  }
  if (!is.list(x)) {
    abort_lna(
      "normalization_provenance must be a named list",
      .subclass = "lna_error_analysis_contract",
      location = "analysis_contract:normalization_provenance"
    )
  }
  out <- list(
    method = x$method %||% "unspecified",
    fit_scope = x$fit_scope %||% "unknown",
    source_id = x$source_id %||% "unrecorded",
    calibration_identity = x$calibration_identity %||% NULL,
    fit_id = x$fit_id %||% NULL
  )
  if (!is.null(out$calibration_identity)) {
    out$calibration_identity <- .analysis_validate_calibration(
      out$calibration_identity,
      "normalization_provenance.calibration_identity",
      canonicalize = canonicalize_checksums
    )
  }
  .analysis_validate_normalization(out)
  out
}

.analysis_validate_calibration <- function(x, field = "calibration_identity",
                                           canonicalize = FALSE) {
  if (!is.list(x)) {
    abort_lna(
      paste0(field, " is required and must be a named list"),
      .subclass = "lna_error_analysis_contract",
      location = paste0("analysis_contract:", field)
    )
  }
  id <- .analysis_scalar_string(x$id, paste0(field, ".id"))
  checksum <- .lna_canonical_sha256(
    x$checksum, paste0(field, ".checksum"),
    canonicalize = canonicalize,
    location = paste0("analysis_contract:", field, ".checksum")
  )
  relationship <- .analysis_scalar_string(
    x$relationship, paste0(field, ".relationship")
  )
  if (!identical(relationship, "external-disjoint")) {
    abort_lna(
      paste0(field, " must declare relationship='external-disjoint'"),
      .subclass = "lna_error_analysis_contract",
      location = paste0("analysis_contract:", field)
    )
  }
  list(id = id, checksum = checksum, relationship = relationship)
}

.analysis_validate_normalization <- function(x) {
  if (!is.list(x)) {
    abort_lna(
      "normalization_provenance must be a named list",
      .subclass = "lna_error_analysis_contract",
      location = "analysis_contract:normalization_provenance"
    )
  }
  .analysis_scalar_string(x$method, "normalization_provenance.method")
  scope <- .analysis_scalar_string(
    x$fit_scope, "normalization_provenance.fit_scope"
  )
  if (!scope %in% .analysis_fit_scopes) {
    abort_lna(
      paste0("normalization_provenance.fit_scope is invalid: ", scope),
      .subclass = "lna_error_analysis_contract",
      location = "analysis_contract:normalization_provenance.fit_scope"
    )
  }
  .analysis_scalar_string(x$source_id, "normalization_provenance.source_id")
  if (identical(scope, "external-calibration")) {
    .analysis_validate_calibration(
      x$calibration_identity, "normalization_provenance.calibration_identity"
    )
  } else if (!is.null(x$calibration_identity)) {
    abort_lna(
      "normalization calibration_identity requires external-calibration fit_scope",
      .subclass = "lna_error_analysis_contract",
      location = "analysis_contract:normalization_provenance.calibration_identity"
    )
  }
  if (scope %in% setdiff(.analysis_fit_scopes, c(
    "not-applicable", "unknown", "fixed", "external-calibration", "pointwise"
  ))) {
    .analysis_scalar_string(x$fit_id, "normalization_provenance.fit_id")
  }
  invisible(TRUE)
}

.analysis_external_identity <- function(id, checksum) {
  list(
    id = .analysis_scalar_string(id, "calibration_id"),
    checksum = .lna_canonical_sha256(
      checksum, "calibration_checksum",
      location = "analysis_contract:calibration_checksum"
    ),
    relationship = "external-disjoint"
  )
}

.parcel_mean_analysis_contract <- function(p) {
  fit_scope <- .analysis_scalar_string(p$fit_scope, "parcellation_fit_scope")
  if (!fit_scope %in% c(
    "fixed", "external-calibration", "fold", "run", "session", "subject"
  )) {
    abort_lna(
      paste0("invalid parcel.mean fit_scope: ", fit_scope),
      .subclass = "lna_error_analysis_contract",
      location = "analysis_contract:parcel.mean_fit_scope"
    )
  }
  if (fit_scope %in% c("fold", "run", "session", "subject") &&
      (!is.list(p$fit_provenance) || !length(p$fit_provenance) ||
       is.null(names(p$fit_provenance)) || any(!nzchar(names(p$fit_provenance))))) {
    abort_lna(
      "fitted parcel.mean declarations require named fit_provenance",
      .subclass = "lna_error_analysis_contract",
      location = "analysis_contract:parcel.mean_fit_provenance"
    )
  }
  calibration <- if (identical(fit_scope, "external-calibration")) {
    .analysis_external_identity(p$calibration_id, p$calibration_checksum)
  } else {
    NULL
  }
  list(
    schema_version = "1.0.0",
    parcellation_fit_scope = fit_scope,
    quantizer_fit_scope = "not-applicable",
    lossy_support = "parcel-at-same-TR",
    coding_dependency = "none",
    independence_unit = "pointwise",
    calibration_identity = calibration,
    normalization_provenance = .analysis_normalization_record(
      p$normalization_provenance
    ),
    loss_profile = list(
      kind = "spatial-parcel-mean",
      map_checksum = p$map_checksum,
      weight_policy = p$weight_policy
    ),
    transformed_units = list(),
    fitted_operations = list()
  )
}

.qpp_analysis_contract <- function(p) {
  policy <- .analysis_scalar_string(p$quantizer_policy, "quantizer_policy")
  expected_scope <- switch(
    policy,
    fixed = "fixed",
    `external-calibration` = "external-calibration",
    `fitted-run-maxabs` = "run",
    abort_lna(
      paste0("unknown qpp quantizer policy in analysis contract: ", policy),
      .subclass = "lna_error_analysis_contract",
      location = "analysis_contract:qpp_policy"
    )
  )
  if (!identical(p$fit_scope, expected_scope)) {
    abort_lna(
      paste0(
        "qpp quantizer policy '", policy, "' contradicts fit_scope '",
        p$fit_scope %||% "<missing>", "'"
      ),
      .subclass = "lna_error_analysis_contract",
      location = "analysis_contract:qpp_fit_scope"
    )
  }
  expected_profile <- if (policy %in% c("fixed", "external-calibration")) {
    declared_profile <- p$profile_id %||% ""
    if (identical(declared_profile, "analysis-safe")) {
      .qpp_validate_parcel_binding_fields(
        p$parcel_binding %||% NULL,
        location = "analysis_contract:qpp_parcel_binding"
      )
      "analysis-safe"
    } else if (identical(declared_profile, "parcel-series-qpp") &&
               is.null(p$parcel_binding)) {
      "parcel-series-qpp"
    } else {
      abort_lna(
        "fixed/external qpp must declare either bound analysis-safe or standalone parcel-series-qpp",
        .subclass = "lna_error_analysis_contract",
        location = "analysis_contract:qpp_profile"
      )
    }
  } else {
    "fitted-run"
  }
  if (!identical(p$profile_id, expected_profile)) {
    abort_lna(
      paste0("qpp profile_id contradicts quantizer policy '", policy, "'"),
      .subclass = "lna_error_analysis_contract",
      location = "analysis_contract:qpp_profile"
    )
  }
  if (identical(policy, "fitted-run-maxabs")) {
    .analysis_scalar_string(p$fit_id, "qpp.fit_id")
  }
  calibration <- if (identical(policy, "external-calibration")) {
    .analysis_external_identity(p$calibration_id, p$calibration_checksum)
  } else {
    NULL
  }
  list(
    schema_version = "1.0.0",
    parcellation_fit_scope = "not-applicable",
    quantizer_fit_scope = expected_scope,
    lossy_support = "pointwise",
    coding_dependency = "chunk",
    independence_unit = "pointwise",
    calibration_identity = calibration,
    normalization_provenance = .analysis_normalization_record(
      p$normalization_provenance
    ),
    loss_profile = list(
      kind = "pointwise-error-bounded-quantization",
      error_bound = "one-half-parcel-step",
      parcel_steps_checksum = p$parcel_steps_checksum,
      clipping = FALSE
    ),
    transformed_units = list(),
    fitted_operations = list()
  )
}

.analysis_validate_transformed_units <- function(x) {
  if (!is.list(x)) {
    abort_lna(
      "transformed_units must be an array of unit declarations",
      .subclass = "lna_error_analysis_contract",
      location = "analysis_contract:transformed_units"
    )
  }
  for (entry in x) {
    if (!is.list(entry)) {
      abort_lna(
        "each transformed_units entry must be a named list",
        .subclass = "lna_error_analysis_contract",
        location = "analysis_contract:transformed_units"
      )
    }
    scope <- .analysis_scalar_string(entry$scope, "transformed_units.scope")
    if (!scope %in% .analysis_unit_scopes || identical(scope, "pointwise")) {
      abort_lna(
        paste0("transformed_units scope is invalid: ", scope),
        .subclass = "lna_error_analysis_contract",
        location = "analysis_contract:transformed_units.scope"
      )
    }
    ids <- entry$ids
    if (!is.character(ids) || !length(ids) || anyNA(ids) || any(!nzchar(ids)) ||
        anyDuplicated(ids)) {
      abort_lna(
        "transformed_units ids must be unique non-empty strings",
        .subclass = "lna_error_analysis_contract",
        location = "analysis_contract:transformed_units.ids"
      )
    }
  }
  invisible(TRUE)
}

.analysis_validate_fitted_operations <- function(x) {
  if (!is.list(x)) {
    abort_lna(
      "fitted_operations must be an array",
      .subclass = "lna_error_analysis_contract",
      location = "analysis_contract:fitted_operations"
    )
  }
  for (op in x) {
    if (!is.list(op)) {
      abort_lna(
        "each fitted operation must be a named list",
        .subclass = "lna_error_analysis_contract",
        location = "analysis_contract:fitted_operations"
      )
    }
    .analysis_scalar_string(op$operation_id, "fitted_operations.operation_id")
    scope <- .analysis_scalar_string(op$fit_scope, "fitted_operations.fit_scope")
    if (!scope %in% .analysis_fit_scopes || scope %in% c("not-applicable", "unknown")) {
      abort_lna(
        paste0("fitted operation has invalid fit_scope: ", scope),
        .subclass = "lna_error_analysis_contract",
        location = "analysis_contract:fitted_operations.fit_scope"
      )
    }
    if (identical(scope, "external-calibration")) {
      .analysis_validate_calibration(
        op$calibration_identity, "fitted_operations.calibration_identity"
      )
    } else if (!scope %in% c("fixed", "pointwise")) {
      .analysis_scalar_string(op$fit_id, "fitted_operations.fit_id")
    }
  }
  invisible(TRUE)
}

.analysis_validate_declaration <- function(x, transform_id) {
  if (!is.list(x)) {
    abort_lna(
      paste0("transform ", transform_id, " is missing analysis_contract"),
      .subclass = "lna_error_analysis_contract",
      location = paste0("analysis_contract:", transform_id)
    )
  }
  if (!is.null(x$allowed_split_unit)) {
    abort_lna(
      paste0(transform_id, " must not assert allowed_split_unit; it is derived"),
      .subclass = "lna_error_analysis_contract",
      location = paste0("analysis_contract:", transform_id)
    )
  }
  required <- c(
    "schema_version", "parcellation_fit_scope", "quantizer_fit_scope",
    "lossy_support", "coding_dependency", "independence_unit",
    "normalization_provenance", "loss_profile", "transformed_units",
    "fitted_operations"
  )
  missing <- required[vapply(required, function(nm) is.null(x[[nm]]), logical(1))]
  if (length(missing)) {
    abort_lna(
      paste0(transform_id, " analysis contract is missing: ",
             paste(missing, collapse = ", ")),
      .subclass = "lna_error_analysis_contract",
      location = paste0("analysis_contract:", transform_id)
    )
  }
  if (!identical(x$schema_version, "1.0.0")) {
    abort_lna(
      paste0(transform_id, " has unsupported analysis contract version"),
      .subclass = "lna_error_analysis_contract",
      location = paste0("analysis_contract:", transform_id)
    )
  }
  for (field in c("parcellation_fit_scope", "quantizer_fit_scope")) {
    value <- .analysis_scalar_string(x[[field]], field)
    if (!value %in% .analysis_fit_scopes) {
      abort_lna(
        paste0(transform_id, " has invalid ", field, ": ", value),
        .subclass = "lna_error_analysis_contract",
        location = paste0("analysis_contract:", transform_id, ":", field)
      )
    }
  }
  if (!.analysis_scalar_string(x$lossy_support, "lossy_support") %in%
      .analysis_lossy_supports) {
    abort_lna(
      paste0(transform_id, " has invalid lossy_support"),
      .subclass = "lna_error_analysis_contract",
      location = paste0("analysis_contract:", transform_id, ":lossy_support")
    )
  }
  if (!.analysis_scalar_string(x$coding_dependency, "coding_dependency") %in%
      .analysis_coding_dependencies) {
    abort_lna(
      paste0(transform_id, " has invalid coding_dependency"),
      .subclass = "lna_error_analysis_contract",
      location = paste0("analysis_contract:", transform_id, ":coding_dependency")
    )
  }
  if (!.analysis_scalar_string(x$independence_unit, "independence_unit") %in%
      .analysis_unit_scopes) {
    abort_lna(
      paste0(transform_id, " has invalid independence_unit"),
      .subclass = "lna_error_analysis_contract",
      location = paste0("analysis_contract:", transform_id, ":independence_unit")
    )
  }
  .analysis_validate_normalization(x$normalization_provenance)
  if (!is.list(x$loss_profile) ||
      !is.character(x$loss_profile$kind) || length(x$loss_profile$kind) != 1L ||
      !nzchar(x$loss_profile$kind)) {
    abort_lna(
      paste0(transform_id, " loss_profile requires a non-empty kind"),
      .subclass = "lna_error_analysis_contract",
      location = paste0("analysis_contract:", transform_id, ":loss_profile")
    )
  }
  .analysis_validate_transformed_units(x$transformed_units)
  .analysis_validate_fitted_operations(x$fitted_operations)
  invisible(TRUE)
}

.analysis_assert_contract_match <- function(actual, expected, transform_id) {
  fields <- names(expected)
  mismatched <- fields[!vapply(fields, function(field) {
    identical(actual[[field]], expected[[field]])
  }, logical(1))]
  if (length(mismatched)) {
    abort_lna(
      paste0(
        transform_id, " analysis declaration contradicts descriptor fields: ",
        paste(mismatched, collapse = ", ")
      ),
      .subclass = "lna_error_analysis_contract",
      location = paste0("analysis_contract:", transform_id, ":contradiction")
    )
  }
  invisible(TRUE)
}

.analysis_descriptors_from_h5 <- function(h5) {
  if (!inherits(h5, "H5File") || !h5$exists("transforms")) {
    abort_lna(
      "archive has no readable /transforms group",
      .subclass = "lna_error_analysis_contract",
      location = "analysis_contract:archive"
    )
  }
  group <- h5[["transforms"]]
  on.exit(safe_h5_close(group), add = TRUE)
  discovered <- discover_transforms(group)
  out <- lapply(discovered$name, function(nm) read_json_descriptor(group, nm))
  names(out) <- discovered$name
  out
}

.analysis_descriptors <- function(x) {
  if (inherits(x, "Plan")) return(x$descriptors)
  if (inherits(x, "DataHandle")) {
    if (!is.null(x$plan)) return(x$plan$descriptors)
    if (!is.null(x$h5)) return(.analysis_descriptors_from_h5(x$h5))
  }
  if (inherits(x, "H5File")) return(.analysis_descriptors_from_h5(x))
  if (is.character(x) && length(x) == 1L && file.exists(x)) {
    h5 <- open_h5(x, mode = "r")
    on.exit(close_h5_safely(h5), add = TRUE)
    return(.analysis_descriptors_from_h5(h5))
  }
  if (is.list(x) && is.character(x$type) && length(x$type) == 1L) {
    return(list(transform_00 = x))
  }
  if (is.list(x) && length(x) && all(vapply(x, function(desc) {
    is.list(desc) && is.character(desc$type) && length(desc$type) == 1L
  }, logical(1)))) {
    return(x)
  }
  abort_lna(
    "x must be an LNA path, H5File, DataHandle, Plan, or descriptor list",
    .subclass = "lna_error_analysis_contract",
    location = "analysis_contract:x"
  )
}

.analysis_contract_from_descriptor <- function(desc, name) {
  type <- .analysis_scalar_string(desc$type, paste0(name, ".type"))
  version <- .analysis_scalar_string(desc$version, paste0(name, ".version"))
  transform_id <- paste0(type, "@", version)
  declared <- desc$analysis_contract
  .analysis_validate_declaration(declared, transform_id)
  if (identical(type, "qpp")) {
    expected <- .qpp_analysis_contract(desc$params %||% list())
    .analysis_assert_contract_match(declared, expected, transform_id)
  } else if (identical(type, "parcel.mean")) {
    expected <- .parcel_mean_analysis_contract(desc$params %||% list())
    .analysis_assert_contract_match(declared, expected, transform_id)
  }
  list(name = name, transform_id = transform_id, type = type,
       declaration = declared)
}

.analysis_add_scope_constraint <- function(out, scope, transform, source,
                                           calibration = NULL,
                                           prospective = TRUE) {
  if (scope %in% c("not-applicable", "fixed", "pointwise")) return(out)
  if (identical(scope, "unknown")) {
    abort_lna(
      paste0(transform, " has unknown provenance for ", source),
      .subclass = "lna_error_analysis_contract",
      location = paste0("analysis_contract:", transform, ":", source)
    )
  }
  if (identical(scope, "external-calibration")) {
    .analysis_validate_calibration(calibration, paste0(source, ".calibration_identity"))
    return(out)
  }
  out[[length(out) + 1L]] <- list(
    transform = transform, source = source, unit = scope,
    fitted = TRUE,
    strict_prospective_incompatible = isTRUE(prospective) &&
      !identical(scope, "fold"),
    transformed_ids = NULL
  )
  out
}

.analysis_transformed_ids <- function(contract, scope, transform) {
  matches <- Filter(function(x) identical(x$scope, scope), contract$transformed_units)
  if (length(matches) != 1L) {
    abort_lna(
      paste0(
        transform, " lossy_support '", scope,
        "' requires exactly one transformed_units declaration"
      ),
      .subclass = "lna_error_analysis_contract",
      location = paste0("analysis_contract:", transform, ":transformed_units")
    )
  }
  matches[[1L]]$ids
}

.analysis_constraints <- function(entries) {
  out <- list()
  for (entry in entries) {
    d <- entry$declaration
    transform <- entry$transform_id
    out <- .analysis_add_scope_constraint(
      out, d$parcellation_fit_scope, transform, "parcellation_fit_scope",
      d$calibration_identity
    )
    out <- .analysis_add_scope_constraint(
      out, d$quantizer_fit_scope, transform, "quantizer_fit_scope",
      d$calibration_identity
    )
    out <- .analysis_add_scope_constraint(
      out, d$normalization_provenance$fit_scope, transform,
      "normalization_provenance.fit_scope",
      d$normalization_provenance$calibration_identity
    )
    for (op in d$fitted_operations) {
      out <- .analysis_add_scope_constraint(
        out, op$fit_scope, transform,
        paste0("fitted_operation:", op$operation_id),
        op$calibration_identity
      )
    }
    lossy_unit <- if (d$lossy_support %in% c("pointwise", "parcel-at-same-TR")) {
      "pointwise"
    } else {
      d$lossy_support
    }
    if (!identical(lossy_unit, "pointwise")) {
      out[[length(out) + 1L]] <- list(
        transform = transform, source = "lossy_support", unit = lossy_unit,
        fitted = FALSE, strict_prospective_incompatible = FALSE,
        transformed_ids = .analysis_transformed_ids(d, lossy_unit, transform)
      )
    }
  }
  out
}

.analysis_unit_dominates <- function(a, b) {
  hierarchy <- c(
    pointwise = 0L, `trial-epoch` = 1L, run = 2L, session = 3L,
    subject = 4L
  )
  a %in% names(hierarchy) && b %in% names(hierarchy) && hierarchy[[a]] > hierarchy[[b]]
}

.analysis_allowed_units <- function(constraints) {
  units <- unique(vapply(constraints, `[[`, character(1), "unit"))
  if (!length(units)) return("pointwise")
  keep <- !vapply(units, function(unit) {
    any(vapply(setdiff(units, unit), function(other) {
      .analysis_unit_dominates(other, unit)
    }, logical(1)))
  }, logical(1))
  order <- c("pointwise", "trial-epoch", "run", "session", "subject",
             "chunk", "causal-prefix", "fold")
  units <- units[keep]
  units[order(match(units, order), units)]
}

#' Derive the analysis-safety contract for an LNA archive
#'
#' The result is derived from typed per-transform declarations. In particular,
#' `allowed_split_unit` is never accepted from archive metadata.
#'
#' @param x An LNA path, open HDF5 file, `DataHandle`, `Plan`, transform
#'   descriptor, or list of transform descriptors.
#' @return An object of class `lna_analysis_contract`. Its
#'   `allowed_split_unit` is a machine term: `"pointwise"` means individual
#'   observations (TRs for the QPP profile) may be assigned independently;
#'   `"run"`, `"session"`, or `"subject"` means each such unit must remain
#'   whole. Exact coding dependencies are reported separately and do not widen
#'   lossy or fitted support.
#' @export
analysis_contract <- function(x) {
  if (inherits(x, "lna_analysis_contract")) return(x)
  descriptors <- .analysis_descriptors(x)
  if (!length(descriptors)) {
    abort_lna(
      "archive contains no transform descriptors",
      .subclass = "lna_error_analysis_contract",
      location = "analysis_contract:descriptors"
    )
  }
  entries <- Map(
    function(desc, name) .analysis_contract_from_descriptor(desc, name),
    descriptors, names(descriptors) %||% paste0("transform_", seq_along(descriptors))
  )
  constraints <- .analysis_constraints(entries)
  coding <- lapply(entries, function(x) list(
    transform = x$transform_id,
    dependency = x$declaration$coding_dependency
  ))
  structure(list(
    schema_version = "1.0.0",
    transforms = entries,
    constraints = constraints,
    coding_dependencies = coding,
    allowed_split_unit = .analysis_allowed_units(constraints)
  ), class = "lna_analysis_contract")
}

.analysis_fold_column <- function(unit) {
  switch(unit,
         `trial-epoch` = "trial_epoch",
         `causal-prefix` = "causal_prefix",
         unit)
}

.analysis_fold_keys <- function(folds, unit) {
  column <- .analysis_fold_column(unit)
  if (!column %in% names(folds)) {
    abort_lna(
      paste0("fold table requires column '", column, "' for ", unit, " support"),
      .subclass = "lna_error_split_incompatible",
      location = paste0("assert_split_compatible:", unit),
      support_unit = unit
    )
  }
  if (anyNA(folds[[column]]) || any(!nzchar(as.character(folds[[column]])))) {
    abort_lna(
      paste0("fold column '", column, "' contains missing or empty unit ids"),
      .subclass = "lna_error_split_incompatible",
      location = paste0("assert_split_compatible:", unit),
      support_unit = unit
    )
  }
  ancestors <- switch(
    unit,
    subject = "subject",
    session = c("subject", "session"),
    run = c("subject", "session", "run"),
    `trial-epoch` = c("subject", "session", "run", "trial_epoch"),
    chunk = c("subject", "session", "run", "chunk"),
    `causal-prefix` = c("subject", "session", "run", "causal_prefix"),
    column
  )
  ancestors <- ancestors[ancestors %in% names(folds)]
  values <- lapply(folds[ancestors], as.character)
  keys <- do.call(paste, c(values, sep = "\u001f"))
  list(keys = keys, labels = as.character(folds[[column]]), column = column)
}

.analysis_abort_split <- function(constraint, units, detail = NULL) {
  units <- unique(as.character(units))
  message <- paste0(
    constraint$transform, " ", constraint$source,
    " requires complete ", constraint$unit, " units; incompatible units: ",
    paste(units, collapse = ", ")
  )
  if (!is.null(detail)) message <- paste0(message, " (", detail, ")")
  abort_lna(
    message,
    .subclass = "lna_error_split_incompatible",
    location = paste0("assert_split_compatible:", constraint$transform),
    transform = constraint$transform,
    support_unit = constraint$unit,
    incompatible_units = units,
    source = constraint$source
  )
}

#' Assert that proposed folds are compatible with an LNA analysis contract
#'
#' Coding dependencies are reported by `analysis_contract()` but deliberately
#' do not constrain folds when their inverse is exact.
#'
#' @param x Anything accepted by [analysis_contract()].
#' @param folds A data frame with one row per observation, a non-missing `fold`
#'   column, and columns naming any required units (`subject`, `session`, `run`,
#'   `trial_epoch`, `chunk`, `causal_prefix`, or `fit_fold`).
#' @param strict_prospective If `TRUE`, reject transforms fitted on the
#'   evaluation subject/session/run even when the fitted unit is kept whole.
#' @return Invisibly, a compatibility receipt containing the derived machine
#'   `allowed_split_unit`, the checked constraints, and `compatible = TRUE`.
#' @export
assert_split_compatible <- function(x, folds, strict_prospective = FALSE) {
  contract <- analysis_contract(x)
  if (!is.data.frame(folds) || !nrow(folds) || !"fold" %in% names(folds) ||
      anyNA(folds$fold) || any(!nzchar(as.character(folds$fold)))) {
    abort_lna(
      "folds must be a non-empty data frame with a non-missing fold column",
      .subclass = "lna_error_split_incompatible",
      location = "assert_split_compatible:folds"
    )
  }
  if (!is.logical(strict_prospective) || length(strict_prospective) != 1L ||
      is.na(strict_prospective)) {
    abort_lna(
      "strict_prospective must be one TRUE or FALSE value",
      .subclass = "lna_error_split_incompatible",
      location = "assert_split_compatible:strict_prospective"
    )
  }
  if ("observation_id" %in% names(folds) &&
      (anyNA(folds$observation_id) || anyDuplicated(folds$observation_id))) {
    abort_lna(
      "observation_id must be complete and unique",
      .subclass = "lna_error_split_incompatible",
      location = "assert_split_compatible:observation_id"
    )
  }

  checked <- list()
  for (constraint in contract$constraints) {
    if (isTRUE(strict_prospective) &&
        isTRUE(constraint$strict_prospective_incompatible)) {
      .analysis_abort_split(
        constraint, constraint$unit,
        "strict prospective evaluation forbids same-data fitting"
      )
    }
    if (identical(constraint$unit, "fold")) {
      if (!isTRUE(constraint$fitted)) {
        missing_units <- setdiff(
          constraint$transformed_ids %||% character(),
          unique(as.character(folds$fold))
        )
        if (length(missing_units)) {
          .analysis_abort_split(
            constraint, missing_units,
            "declared transformed fold units are absent from the fold table"
          )
        }
        checked[[length(checked) + 1L]] <- constraint
        next
      }
      if (!"fit_fold" %in% names(folds) || anyNA(folds$fit_fold) ||
          !identical(as.character(folds$fit_fold), as.character(folds$fold))) {
        bad <- if ("fit_fold" %in% names(folds)) {
          unique(as.character(folds$fold[
            is.na(folds$fit_fold) |
              as.character(folds$fit_fold) != as.character(folds$fold)
          ]))
        } else {
          unique(as.character(folds$fold))
        }
        .analysis_abort_split(
          constraint, bad,
          "fit_fold must identify the training fit applied to each held-out fold"
        )
      }
      checked[[length(checked) + 1L]] <- constraint
      next
    }
    unit <- .analysis_fold_keys(folds, constraint$unit)
    if (!is.null(constraint$transformed_ids)) {
      observed <- unique(unit$labels)
      mismatch <- setdiff(constraint$transformed_ids, observed)
      if (length(mismatch)) {
        .analysis_abort_split(
          constraint, mismatch,
          "declared transformed units are absent from the fold table"
        )
      }
    }
    folds_per_unit <- split(as.character(folds$fold), unit$keys)
    split_keys <- names(folds_per_unit)[vapply(
      folds_per_unit, function(value) length(unique(value)) > 1L, logical(1)
    )]
    if (length(split_keys)) {
      bad_labels <- unique(unit$labels[unit$keys %in% split_keys])
      .analysis_abort_split(constraint, bad_labels)
    }
    checked[[length(checked) + 1L]] <- constraint
  }
  invisible(structure(list(
    compatible = TRUE,
    allowed_split_unit = contract$allowed_split_unit,
    checked_constraints = checked,
    ignored_coding_dependencies = contract$coding_dependencies
  ), class = "lna_split_compatibility_receipt"))
}
