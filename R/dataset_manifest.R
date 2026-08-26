# Typed, versioned manifest contract for portable LNA exchange datasets.

.lna_manifest_schema_version <- "1.1.0"
.lna_manifest_legacy_schema_version <- "1.0.0"
.lna_manifest_null_encoding <- "n/a"
.lna_manifest_supported_major <- 1L
.lna_manifest_column_types <- c(
  "string", "integer", "number", "boolean", "json", "path", "uri",
  "checksum"
)
.lna_manifest_status_values <- c(
  "planned", "writing", "complete", "failed", "stale"
)

.manifest_abort <- function(message, field, ...) {
  abort_lna(
    message,
    ...,
    field = field,
    .subclass = "lna_error_manifest",
    location = paste0("lna_manifest:", field)
  )
}

.manifest_scalar_string <- function(x, field) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    .manifest_abort(paste0(field, " must be one non-empty string"), field)
  }
  enc2utf8(x)
}

.manifest_schema_major <- function(version) {
  version <- .manifest_scalar_string(version, "manifest_schema_version")
  if (!grepl("^[0-9]+[.][0-9]+[.][0-9]+$", version)) {
    .manifest_abort(
      "manifest_schema_version must use semantic versioning (major.minor.patch)",
      "manifest_schema_version"
    )
  }
  as.integer(strsplit(version, ".", fixed = TRUE)[[1L]][[1L]])
}

.manifest_schema_version_supported <- function(version) {
  identical(.manifest_schema_major(version), .lna_manifest_supported_major)
}

.manifest_column <- function(
    type,
    nullable,
    unit = NULL,
    description,
    role = "data",
    values = NULL) {
  column <- list(
    type = type,
    nullable = isTRUE(nullable),
    unit = unit,
    description = description,
    role = role
  )
  if (!is.null(values)) column$values <- values
  column
}

.lna_manifest_core_columns <- function() {
  list(
    record_id = .manifest_column(
      "string", FALSE, description = "Stable source-image record identity",
      role = "identity"
    ),
    source_path = .manifest_column(
      "path", TRUE, description = "Source path relative to its BIDS dataset",
      role = "source"
    ),
    source_uri = .manifest_column(
      "uri", TRUE, description = "Absolute source URI when no local path exists",
      role = "source"
    ),
    source_checksum = .manifest_column(
      "checksum", FALSE, description = "Checksum of the represented source image",
      role = "integrity"
    ),
    source_extension = .manifest_column(
      "string", FALSE, description = "Original source imaging extension",
      role = "source"
    ),
    source_bytes = .manifest_column(
      "integer", TRUE, unit = "byte",
      description = "Byte size of the represented source image", role = "size"
    ),
    lna_path = .manifest_column(
      "path", FALSE, description = "LNA artifact path relative to archive root",
      role = "artifact"
    ),
    lna_checksum = .manifest_column(
      "checksum", TRUE, description = "Checksum of the finalized LNA artifact",
      role = "integrity"
    ),
    lna_bytes = .manifest_column(
      "integer", TRUE, unit = "byte",
      description = "Byte size of the finalized LNA artifact", role = "size"
    ),
    datatype = .manifest_column(
      "string", FALSE, description = "BIDS datatype directory", role = "bids"
    ),
    suffix = .manifest_column(
      "string", FALSE, description = "BIDS suffix before the LNA extension",
      role = "bids"
    ),
    scope = .manifest_column(
      "string", FALSE, description = "BIDS raw or derivative scope", role = "bids",
      values = c("raw", "derivative")
    ),
    pipeline = .manifest_column(
      "string", TRUE, description = "Derivative pipeline identity", role = "bids"
    ),
    entity_sub = .manifest_column(
      "string", FALSE, description = "Exact BIDS subject entity value",
      role = "bids-entity"
    ),
    codec_profile_id = .manifest_column(
      "string", FALSE, description = "Declared codec profile identity",
      role = "codec"
    ),
    loss_mode = .manifest_column(
      "string", FALSE, description = "Loss declaration for the codec profile",
      role = "codec", values = c("lossless", "lossy", "unspecified")
    ),
    codec_profile_json = .manifest_column(
      "json", FALSE, description = "Complete codec profile and parameters",
      role = "codec"
    ),
    shared_basis_json = .manifest_column(
      "json", FALSE,
      description = "Exact shared-basis identity and fit provenance, or status none",
      role = "codec-resource"
    ),
    associations_json = .manifest_column(
      "json", FALSE, description = "Portable companion and association closure",
      role = "association"
    ),
    metadata_json = .manifest_column(
      "json", FALSE, description = "Resolved source metadata", role = "metadata"
    ),
    metadata_provenance_json = .manifest_column(
      "json", FALSE,
      description = "Field-level provenance for resolved metadata",
      role = "provenance"
    ),
    n_volumes = .manifest_column(
      "integer", TRUE, unit = "volume",
      description = "Number of temporal volumes", role = "temporal"
    ),
    repetition_time_seconds = .manifest_column(
      "number", TRUE, unit = "second",
      description = "Resolved repetition time", role = "temporal"
    ),
    reference_identifier = .manifest_column(
      "string", FALSE, description = "Provenance-asserted reference identity",
      role = "spatial"
    ),
    reference_kind = .manifest_column(
      "string", FALSE, description = "Reference assertion kind", role = "spatial"
    ),
    reference_origin = .manifest_column(
      "string", FALSE, description = "Origin of the reference assertion",
      role = "provenance"
    ),
    reference_status = .manifest_column(
      "string", FALSE, description = "Status of the reference assertion",
      role = "provenance"
    ),
    reference_json = .manifest_column(
      "json", FALSE, description = "Complete reference assertion and receipts",
      role = "provenance"
    ),
    grid_id = .manifest_column(
      "checksum", FALSE, description = "Canonical stored-grid identity",
      role = "spatial"
    ),
    support_id = .manifest_column(
      "checksum", TRUE, description = "Canonical stored-support identity",
      role = "spatial"
    ),
    status = .manifest_column(
      "string", FALSE, description = "Artifact packaging state", role = "state",
      values = .lna_manifest_status_values
    )
  )
}

.manifest_insert_named_after <- function(values, name, value, after) {
  position <- match(after, names(values))
  if (is.na(position)) return(values)
  append(values, stats::setNames(list(value), name), after = position)
}

.manifest_migrate_legacy_schema <- function(schema) {
  result <- list(schema = schema, added_columns = character())
  if (!is.list(schema) || is.data.frame(schema) ||
      !identical(schema$manifest_schema_version, .lna_manifest_legacy_schema_version) ||
      !is.list(schema$columns) || is.data.frame(schema$columns) ||
      "shared_basis_json" %in% names(schema$columns)) {
    return(result)
  }
  order <- schema$column_order
  if (is.list(order)) order <- unlist(order, use.names = FALSE)
  legacy_core <- setdiff(names(.lna_manifest_core_columns()), "shared_basis_json")
  if (!is.character(order) || !identical(order, names(schema$columns)) ||
      length(setdiff(legacy_core, names(schema$columns)))) {
    return(result)
  }
  definition <- .lna_manifest_core_columns()$shared_basis_json
  schema$columns <- .manifest_insert_named_after(
    schema$columns, "shared_basis_json", definition, "codec_profile_json"
  )
  position <- match("codec_profile_json", order)
  schema$column_order <- append(order, "shared_basis_json", after = position)
  result$schema <- schema
  result$added_columns <- "shared_basis_json"
  result
}

.manifest_migrate_legacy_table <- function(raw, schema, added_columns) {
  if (!length(added_columns)) return(raw)
  expected_legacy <- setdiff(schema$column_order, added_columns)
  if (!identical(names(raw), expected_legacy)) {
    .manifest_abort(
      "legacy TSV columns must exactly match its declared pre-migration schema",
      "column_order"
    )
  }
  if ("shared_basis_json" %in% added_columns) {
    raw$shared_basis_json <- rep(
      .manifest_canonical_json(list(status = "none"), "shared_basis_json"),
      nrow(raw)
    )
  }
  raw[schema$column_order]
}

.manifest_entity_column <- function(entity) {
  .manifest_column(
    "string", TRUE,
    description = paste0("Exact BIDS entity value for '", entity, "'"),
    role = "bids-entity"
  )
}

.manifest_is_relative_path <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    return(FALSE)
  }
  if (grepl("[\\\\]", path) || startsWith(path, "/") ||
      grepl("^[A-Za-z]:", path) || grepl("[[:cntrl:]]", path)) {
    return(FALSE)
  }
  components <- strsplit(path, "/", fixed = TRUE)[[1L]]
  length(components) > 0L &&
    all(nzchar(components)) &&
    !any(components %in% c(".", ".."))
}

.manifest_validate_path <- function(path, field) {
  if (!.manifest_is_relative_path(path)) {
    .manifest_abort(
      paste0(field, " must be a traversal-safe archive-relative path"), field
    )
  }
  enc2utf8(path)
}

.manifest_validate_json_value <- function(value, field) {
  tryCatch(
    .source_record_json_value(value, field),
    lna_error_source_record = function(e) {
      .manifest_abort(conditionMessage(e), field)
    }
  )
  value
}

.manifest_canonical_json <- function(value, field) {
  .manifest_validate_json_value(value, field)
  jsonlite::toJSON(
    value, auto_unbox = TRUE, null = "null", digits = NA,
    matrix = "rowmajor", pretty = FALSE
  )
}

.manifest_column_definition <- function(definition, name) {
  field <- paste0("columns.", name)
  if (!is.list(definition) || is.data.frame(definition)) {
    .manifest_abort(paste0(field, " must be an object"), field)
  }
  required <- c("type", "nullable", "unit", "description", "role")
  missing <- setdiff(required, names(definition))
  if (length(missing)) {
    .manifest_abort(
      paste0(field, " is missing: ", paste(missing, collapse = ", ")), field
    )
  }
  type <- .manifest_scalar_string(definition$type, paste0(field, ".type"))
  if (!type %in% .lna_manifest_column_types) {
    .manifest_abort(
      paste0(field, ".type is unsupported: ", type), paste0(field, ".type")
    )
  }
  if (!is.logical(definition$nullable) || length(definition$nullable) != 1L ||
      is.na(definition$nullable)) {
    .manifest_abort(
      paste0(field, ".nullable must be TRUE or FALSE"),
      paste0(field, ".nullable")
    )
  }
  unit <- definition$unit
  if (!is.null(unit)) unit <- .manifest_scalar_string(unit, paste0(field, ".unit"))
  description <- .manifest_scalar_string(
    definition$description, paste0(field, ".description")
  )
  role <- .manifest_scalar_string(definition$role, paste0(field, ".role"))
  values <- definition$values %||% NULL
  if (!is.null(values)) {
    if (is.list(values)) values <- unlist(values, use.names = FALSE)
    if (!is.character(values) || !length(values) || anyNA(values) ||
        any(!nzchar(values)) || anyDuplicated(values)) {
      .manifest_abort(
        paste0(field, ".values must be unique non-empty strings"),
        paste0(field, ".values")
      )
    }
  }
  normalized <- list(
    type = type,
    nullable = definition$nullable,
    unit = unit,
    description = description,
    role = role
  )
  if (!is.null(values)) normalized$values <- values
  extras <- definition[setdiff(
    names(definition), c("type", "nullable", "unit", "description", "role", "values")
  )]
  c(normalized, extras)
}

.manifest_validate_schema <- function(schema) {
  if (!is.list(schema) || is.data.frame(schema)) {
    .manifest_abort("manifest JSON sidecar must contain an object", "schema")
  }
  required <- c(
    "manifest_schema_version", "table", "null_encoding", "row_semantics",
    "column_order", "columns", "entity_semantics", "identity_algorithms",
    "checksum_algorithms"
  )
  missing <- setdiff(required, names(schema))
  if (length(missing)) {
    .manifest_abort(
      paste0("manifest schema is missing: ", paste(missing, collapse = ", ")),
      "schema"
    )
  }
  major <- .manifest_schema_major(schema$manifest_schema_version)
  if (!identical(major, .lna_manifest_supported_major)) {
    .manifest_abort(
      paste0(
        "unsupported manifest schema major version ", major,
        "; supported major version is ", .lna_manifest_supported_major
      ),
      "manifest_schema_version"
    )
  }
  if (!identical(schema$table, "lna_manifest.tsv")) {
    .manifest_abort("schema.table must equal lna_manifest.tsv", "table")
  }
  if (!identical(schema$null_encoding, .lna_manifest_null_encoding)) {
    .manifest_abort("schema.null_encoding must equal 'n/a'", "null_encoding")
  }
  if (!identical(schema$row_semantics, "one-row-per-source-image")) {
    .manifest_abort(
      "schema.row_semantics must equal one-row-per-source-image",
      "row_semantics"
    )
  }
  order <- schema$column_order
  if (is.list(order)) order <- unlist(order, use.names = FALSE)
  if (!is.character(order) || !length(order) || anyNA(order) ||
      any(!nzchar(order)) || anyDuplicated(order)) {
    .manifest_abort(
      "schema.column_order must contain unique non-empty column names",
      "column_order"
    )
  }
  if (!is.list(schema$columns) || is.data.frame(schema$columns) ||
      is.null(names(schema$columns)) || anyDuplicated(names(schema$columns))) {
    .manifest_abort("schema.columns must be a uniquely named object", "columns")
  }
  if (!identical(order, names(schema$columns))) {
    .manifest_abort(
      "schema.column_order must exactly match the order of schema.columns",
      "column_order"
    )
  }
  core <- .lna_manifest_core_columns()
  absent_core <- setdiff(names(core), names(schema$columns))
  if (length(absent_core)) {
    .manifest_abort(
      paste0("manifest is missing core column(s): ", paste(absent_core, collapse = ", ")),
      "columns"
    )
  }
  normalized <- lapply(names(schema$columns), function(name) {
    .manifest_column_definition(schema$columns[[name]], name)
  })
  names(normalized) <- names(schema$columns)
  for (name in names(core)) {
    expected <- core[[name]]
    actual <- normalized[[name]]
    for (property in c("type", "nullable", "unit")) {
      if (!identical(actual[[property]], expected[[property]])) {
        .manifest_abort(
          paste0("core column ", name, " has incompatible ", property),
          paste0("columns.", name, ".", property)
        )
      }
    }
    if (!is.null(expected$values) && !identical(actual$values, expected$values)) {
      .manifest_abort(
        paste0("core column ", name, " has incompatible allowed values"),
        paste0("columns.", name, ".values")
      )
    }
  }
  entity_names <- names(normalized)[startsWith(names(normalized), "entity_")]
  for (name in entity_names) {
    if (!grepl("^entity_[a-z][a-z0-9]*$", name) ||
        !identical(normalized[[name]]$type, "string") ||
        !identical(normalized[[name]]$role, "bids-entity")) {
      .manifest_abort(
        paste0(name, " is not a valid optional BIDS entity column"),
        paste0("columns.", name)
      )
    }
  }
  if (!is.list(schema$entity_semantics) ||
      !identical(schema$entity_semantics$prefix, "entity_") ||
      !identical(schema$entity_semantics$comparison, "exact")) {
    .manifest_abort(
      "schema.entity_semantics must declare entity_ prefix and exact comparison",
      "entity_semantics"
    )
  }
  for (field in c("identity_algorithms", "checksum_algorithms")) {
    if (!is.list(schema[[field]]) || is.data.frame(schema[[field]])) {
      .manifest_abort(paste0("schema.", field, " must be an object"), field)
    }
    .manifest_validate_json_value(schema[[field]], field)
  }
  schema$column_order <- order
  schema$columns <- normalized
  schema
}

.manifest_schema <- function(column_definitions) {
  schema <- list(
    manifest_schema_version = .lna_manifest_schema_version,
    table = "lna_manifest.tsv",
    null_encoding = .lna_manifest_null_encoding,
    row_semantics = "one-row-per-source-image",
    column_order = names(column_definitions),
    columns = column_definitions,
    entity_semantics = list(
      prefix = "entity_",
      optional = TRUE,
      comparison = "exact",
      missing = .lna_manifest_null_encoding
    ),
    identity_algorithms = list(
      grid_id = "lna-grid-v1-sha256",
      support_id = "lna-support-v1-sha256",
      reference_id = "provenance-asserted-not-inferred"
    ),
    checksum_algorithms = list(
      source_checksum = "sha256",
      lna_checksum = "sha256"
    )
  )
  .manifest_validate_schema(schema)
}

.manifest_value_is_missing <- function(value, index, type) {
  if (identical(type, "json")) return(is.null(value[[index]]))
  is.na(value[[index]])
}

.manifest_validate_scalar <- function(value, definition, field) {
  type <- definition$type
  if (type == "json") {
    json <- .manifest_canonical_json(value, field)
    return(jsonlite::fromJSON(json, simplifyVector = FALSE))
  }
  if (type == "integer") {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value != floor(value)) {
      .manifest_abort(paste0(field, " must be an integer"), field)
    }
    return(as.double(value))
  }
  if (type == "number") {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value)) {
      .manifest_abort(paste0(field, " must be a finite number"), field)
    }
    return(as.double(value))
  }
  if (type == "boolean") {
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      .manifest_abort(paste0(field, " must be TRUE or FALSE"), field)
    }
    return(value)
  }
  value <- .manifest_scalar_string(value, field)
  if (identical(value, .lna_manifest_null_encoding) ||
      grepl("[\t\r\n]", value)) {
    .manifest_abort(
      paste0(field, " contains the reserved null value or a TSV delimiter"), field
    )
  }
  if (type == "path") value <- .manifest_validate_path(value, field)
  if (type == "uri" && !grepl("^[A-Za-z][A-Za-z0-9+.-]*:", value)) {
    .manifest_abort(paste0(field, " must be an absolute URI"), field)
  }
  if (type == "checksum" && !grepl(.lna_source_checksum_pattern, value)) {
    .manifest_abort(
      paste0(field, " must be sha256:<64 lowercase hexadecimal digits>"), field
    )
  }
  if (!is.null(definition$values) && !value %in% definition$values) {
    .manifest_abort(
      paste0(field, " must be one of: ", paste(definition$values, collapse = ", ")),
      field
    )
  }
  value
}

.manifest_validate_column <- function(values, definition, name, n_rows) {
  if (identical(definition$type, "json")) {
    if (!is.list(values) || length(values) != n_rows) {
      .manifest_abort(
        paste0(name, " must be a JSON list-column with ", n_rows, " values"), name
      )
    }
  } else if (length(values) != n_rows || is.list(values)) {
    .manifest_abort(paste0(name, " must contain ", n_rows, " scalar values"), name)
  }
  normalized <- vector("list", n_rows)
  for (i in seq_len(n_rows)) {
    missing <- .manifest_value_is_missing(values, i, definition$type)
    field <- paste0(name, "[", i, "]")
    if (missing) {
      if (!isTRUE(definition$nullable)) {
        .manifest_abort(paste0(field, " may not be null"), field)
      }
      normalized[i] <- list(NULL)
    } else {
      value <- if (identical(definition$type, "json")) values[[i]] else values[[i]]
      normalized[i] <- list(
        .manifest_validate_scalar(value, definition, field)
      )
    }
  }
  if (identical(definition$type, "json")) return(normalized)
  if (definition$type %in% c("integer", "number")) {
    return(vapply(normalized, function(x) x %||% NA_real_, numeric(1)))
  }
  if (identical(definition$type, "boolean")) {
    return(vapply(normalized, function(x) x %||% NA, logical(1)))
  }
  vapply(normalized, function(x) x %||% NA_character_, character(1))
}

.manifest_validate_associations <- function(associations, row, record_ids) {
  field <- paste0("associations_json[", row, "]")
  if (!is.list(associations) || (!is.null(names(associations)) && length(associations))) {
    .manifest_abort(paste0(field, " must be a JSON array"), field)
  }
  targets <- character()
  for (i in seq_along(associations)) {
    association <- associations[[i]]
    item_field <- paste0(field, "[[", i, "]]")
    if (!is.list(association) || is.null(names(association))) {
      .manifest_abort(paste0(item_field, " must be a JSON object"), item_field)
    }
    for (name in c("source_path", "destination_relative")) {
      value <- association[[name]] %||% NULL
      if (!is.null(value)) {
        .manifest_validate_path(value, paste0(item_field, ".", name))
      }
    }
    uri <- association$source_uri %||% NULL
    if (!is.null(uri)) {
      .manifest_validate_scalar(
        uri,
        .manifest_column(
          "uri", FALSE, description = "Association URI", role = "association"
        ),
        paste0(item_field, ".source_uri")
      )
    }
    checksum <- association$source_checksum %||% NULL
    if (!is.null(checksum)) {
      .manifest_validate_scalar(
        checksum,
        .manifest_column(
          "checksum", FALSE, description = "Association checksum",
          role = "association"
        ),
        paste0(item_field, ".source_checksum")
      )
    }
    target <- association$target_record_id %||% NULL
    if (!is.null(target)) targets <- c(targets, target)
  }
  missing <- setdiff(targets, record_ids)
  if (length(missing)) {
    .manifest_abort(
      paste0(field, " targets missing record(s): ", paste(unique(missing), collapse = ", ")),
      field
    )
  }
  invisible(TRUE)
}

.manifest_validate_semantics <- function(data) {
  if (!nrow(data)) .manifest_abort("manifest must contain at least one row", "data")
  if (any(!grepl("^[A-Za-z0-9][A-Za-z0-9._:-]*$", data$record_id))) {
    .manifest_abort("record_id contains an invalid identifier", "record_id")
  }
  if (any(!endsWith(data$lna_path, ".lna.h5"))) {
    .manifest_abort("lna_path values must end with .lna.h5", "lna_path")
  }
  for (name in c("datatype", "suffix", "entity_sub")) {
    if (any(!grepl("^[A-Za-z0-9][A-Za-z0-9+.-]*$", data[[name]]))) {
      .manifest_abort(paste0(name, " contains an invalid BIDS label"), name)
    }
  }
  pipeline <- data$pipeline[!is.na(data$pipeline)]
  if (length(pipeline) && any(!grepl("^[A-Za-z0-9][A-Za-z0-9+.-]*$", pipeline))) {
    .manifest_abort("pipeline contains an invalid path label", "pipeline")
  }
  entity_columns <- names(data)[startsWith(names(data), "entity_")]
  for (name in entity_columns) {
    values <- data[[name]][!is.na(data[[name]])]
    if (length(values) && any(
      !grepl("^[A-Za-z0-9][A-Za-z0-9+.-]*$", values) |
        grepl("[/\\\\_]", values)
    )) {
      .manifest_abort(paste0(name, " contains an invalid BIDS entity value"), name)
    }
  }
  if (any(data$source_bytes < 0, na.rm = TRUE) ||
      any(data$lna_bytes < 0, na.rm = TRUE)) {
    .manifest_abort("byte sizes must be non-negative", "bytes")
  }
  if (any(data$n_volumes < 1, na.rm = TRUE)) {
    .manifest_abort("n_volumes must be positive when present", "n_volumes")
  }
  record_ids <- data$record_id
  for (i in seq_len(nrow(data))) {
    .manifest_validate_associations(data$associations_json[[i]], i, record_ids)
    profile <- data$codec_profile_json[[i]]
    if (!is.list(profile) ||
        !identical(profile$profile_id, data$codec_profile_id[[i]]) ||
        !identical(profile$loss_mode, data$loss_mode[[i]])) {
      .manifest_abort(
        paste0("codec_profile_json[", i, "] disagrees with flattened codec fields"),
        paste0("codec_profile_json[", i, "]")
      )
    }
    assertion <- tryCatch(
      .reference_assertion_from_raw(
        data$reference_json[[i]], record_id = data$record_id[[i]]
      ),
      lna_error_reference_assertion = function(e) {
        .manifest_abort(
          paste0("reference_json[", i, "] is invalid: ", conditionMessage(e)),
          paste0("reference_json[", i, "]")
        )
      }
    )
    flattened <- c("identifier", "kind", "origin", "status")
    manifest_names <- paste0("reference_", flattened)
    if (any(vapply(seq_along(flattened), function(j) {
      !identical(assertion[[flattened[[j]]]], data[[manifest_names[[j]]]][[i]])
    }, logical(1)))) {
      .manifest_abort(
        paste0("reference_json[", i, "] disagrees with flattened reference fields"),
        paste0("reference_json[", i, "]")
      )
    }
    metadata <- data$metadata_json[[i]]
    metadata_tr <- metadata$RepetitionTime %||% NULL
    flat_tr <- data$repetition_time_seconds[[i]]
    if (is.null(metadata_tr)) {
      if (!is.na(flat_tr)) {
        .manifest_abort(
          paste0("metadata_json[", i, "] lacks the flattened RepetitionTime"),
          paste0("metadata_json[", i, "]")
        )
      }
    } else if (!is.numeric(metadata_tr) || length(metadata_tr) != 1L ||
        !is.finite(metadata_tr) || is.na(flat_tr) ||
        !isTRUE(all.equal(as.double(metadata_tr), flat_tr, tolerance = 0))) {
      .manifest_abort(
        paste0("metadata_json[", i, "] disagrees with repetition_time_seconds"),
        paste0("metadata_json[", i, "]")
      )
    }
  }
  invisible(TRUE)
}

#' Validate an LNA dataset manifest
#'
#' Validates the JSON schema contract before interpreting any table values.
#' Unknown optional columns are accepted when their type, nullability, unit,
#' description, and role are declared in the sidecar.
#'
#' @param manifest An `lna_manifest` object.
#' @return The normalized manifest, invisibly.
#' @export
validate_lna_manifest <- function(manifest) {
  if (!is.list(manifest) || is.data.frame(manifest) ||
      is.null(manifest$schema) || is.null(manifest$data)) {
    .manifest_abort("manifest must contain schema and data", "manifest")
  }
  schema <- .manifest_validate_schema(manifest$schema)
  data <- manifest$data
  if (!is.data.frame(data)) {
    .manifest_abort("manifest$data must be a data frame", "data")
  }
  if (!identical(names(data), schema$column_order)) {
    .manifest_abort(
      "manifest data columns must exactly match schema.column_order", "data"
    )
  }
  n_rows <- nrow(data)
  normalized <- data
  for (name in names(schema$columns)) {
    normalized[[name]] <- .manifest_validate_column(
      data[[name]], schema$columns[[name]], name, n_rows
    )
  }
  if (anyDuplicated(normalized$record_id)) {
    duplicate <- unique(normalized$record_id[duplicated(normalized$record_id)])
    .manifest_abort(
      paste0("duplicate record_id: ", paste(duplicate, collapse = ", ")),
      "record_id", duplicate_values = duplicate
    )
  }
  if (anyDuplicated(normalized$lna_path)) {
    duplicate <- unique(normalized$lna_path[duplicated(normalized$lna_path)])
    .manifest_abort(
      paste0("duplicate lna_path: ", paste(duplicate, collapse = ", ")),
      "lna_path", duplicate_values = duplicate
    )
  }
  if (any(is.na(normalized$source_path) & is.na(normalized$source_uri))) {
    rows <- which(is.na(normalized$source_path) & is.na(normalized$source_uri))
    .manifest_abort(
      paste0("rows require source_path or source_uri: ", paste(rows, collapse = ", ")),
      "source"
    )
  }
  derivative_missing <- normalized$scope == "derivative" & is.na(normalized$pipeline)
  raw_pipeline <- normalized$scope == "raw" & !is.na(normalized$pipeline)
  if (any(derivative_missing) || any(raw_pipeline)) {
    .manifest_abort(
      "pipeline is required for derivative rows and null for raw rows", "pipeline"
    )
  }
  complete <- normalized$status == "complete"
  if (any(complete & (is.na(normalized$lna_checksum) | is.na(normalized$lna_bytes)))) {
    .manifest_abort(
      "complete rows require lna_checksum and lna_bytes", "status"
    )
  }
  .manifest_validate_semantics(normalized)
  normalized <- normalized[order(normalized$lna_path, normalized$record_id), , drop = FALSE]
  rownames(normalized) <- NULL
  structure(
    list(schema = schema, data = normalized),
    class = c("lna_manifest", "list")
  )
}

.manifest_extra_columns <- function(extra_columns, n_rows) {
  if (!is.list(extra_columns) || is.data.frame(extra_columns)) {
    .manifest_abort("extra_columns must be a named list", "extra_columns")
  }
  if (!length(extra_columns)) return(list(definitions = list(), values = list()))
  if (is.null(names(extra_columns)) || any(!nzchar(names(extra_columns))) ||
      anyDuplicated(names(extra_columns))) {
    .manifest_abort("extra_columns must have unique names", "extra_columns")
  }
  reserved <- names(.lna_manifest_core_columns())
  if (length(intersect(names(extra_columns), reserved)) ||
      any(startsWith(names(extra_columns), "entity_"))) {
    .manifest_abort(
      "extra_columns may not replace core or entity columns", "extra_columns"
    )
  }
  definitions <- list()
  values <- list()
  for (name in sort(names(extra_columns))) {
    item <- extra_columns[[name]]
    if (!is.list(item) || is.null(item$values)) {
      .manifest_abort(
        paste0("extra_columns.", name, " requires values and a column definition"),
        paste0("extra_columns.", name)
      )
    }
    definition <- item$definition %||% item[setdiff(names(item), "values")]
    definitions[[name]] <- .manifest_column_definition(definition, name)
    values[[name]] <- .manifest_validate_column(
      item$values, definitions[[name]], name, n_rows
    )
  }
  list(definitions = definitions, values = values)
}

.manifest_state_for <- function(artifact_state, record_id) {
  if (!length(artifact_state)) return(list())
  artifact_state[[record_id]] %||% list()
}

#' Construct a typed manifest from a dataset plan
#'
#' @param plan A validated `lna_dataset_plan`.
#' @param artifact_state Optional named list keyed by record ID. Each item may
#'   supply `source_bytes`, `lna_checksum`, `lna_bytes`, and `status`.
#' @param extra_columns Optional named future columns, each containing `values`
#'   and a `definition` with type, nullability, unit, description, and role.
#' @return A validated `lna_manifest`.
#' @export
lna_manifest <- function(plan, artifact_state = list(), extra_columns = list()) {
  if (!inherits(plan, "lna_dataset_plan")) {
    .manifest_abort("plan must be an lna_dataset_plan", "plan")
  }
  if (!is.list(artifact_state) || is.data.frame(artifact_state) ||
      (length(artifact_state) && (is.null(names(artifact_state)) ||
        any(!nzchar(names(artifact_state))) || anyDuplicated(names(artifact_state))))) {
    .manifest_abort("artifact_state must be uniquely named by record_id", "artifact_state")
  }
  record_ids <- vapply(plan$artifacts, `[[`, character(1), "record_id")
  unknown_state <- setdiff(names(artifact_state), record_ids)
  if (length(unknown_state)) {
    .manifest_abort(
      paste0("artifact_state has unknown record_id(s): ", paste(unknown_state, collapse = ", ")),
      "artifact_state"
    )
  }
  entity_names <- unique(unlist(lapply(plan$artifacts, function(x) names(x$entities))))
  entity_names <- .dataset_plan_entity_names(stats::setNames(as.list(entity_names), entity_names))
  entity_names <- setdiff(entity_names, "sub")
  entity_columns <- paste0("entity_", entity_names)
  definitions <- c(
    .lna_manifest_core_columns(),
    stats::setNames(lapply(entity_names, .manifest_entity_column), entity_columns)
  )
  n_rows <- length(plan$artifacts)
  extras <- .manifest_extra_columns(extra_columns, n_rows)
  definitions <- c(definitions, extras$definitions)

  rows <- lapply(plan$artifacts, function(entry) {
    state <- .manifest_state_for(artifact_state, entry$record_id)
    unknown <- setdiff(
      names(state), c("source_bytes", "lna_checksum", "lna_bytes", "status")
    )
    if (length(unknown)) {
      .manifest_abort(
        paste0("artifact_state.", entry$record_id, " has unknown field(s): ",
          paste(unknown, collapse = ", ")),
        paste0("artifact_state.", entry$record_id)
      )
    }
    tr <- entry$resolved_metadata$RepetitionTime %||% NA_real_
    reference <- entry$reference_assertion
    row <- list(
      record_id = entry$record_id,
      source_path = entry$source$path %||% NA_character_,
      source_uri = entry$source$uri %||% NA_character_,
      source_checksum = entry$source$checksum,
      source_extension = entry$bids$extension,
      source_bytes = state$source_bytes %||% NA_real_,
      lna_path = entry$destination_relative,
      lna_checksum = state$lna_checksum %||% NA_character_,
      lna_bytes = state$lna_bytes %||% NA_real_,
      datatype = entry$bids$datatype,
      suffix = entry$bids$suffix,
      scope = entry$bids$scope,
      pipeline = entry$bids$pipeline %||% NA_character_,
      entity_sub = entry$entities$sub,
      codec_profile_id = entry$codec_profile_id,
      loss_mode = plan$codec_profile$loss_mode,
      codec_profile_json = plan$codec_profile,
      shared_basis_json = entry$shared_basis %||% list(status = "none"),
      associations_json = entry$expected_companions,
      metadata_json = entry$resolved_metadata,
      metadata_provenance_json = entry$metadata_provenance,
      n_volumes = entry$geometry$n_volumes %||% NA_real_,
      repetition_time_seconds = tr,
      reference_identifier = reference$identifier,
      reference_kind = reference$kind,
      reference_origin = reference$origin,
      reference_status = reference$status,
      reference_json = unclass(reference),
      grid_id = entry$geometry$grid_id,
      support_id = entry$geometry$support_id %||% NA_character_,
      status = state$status %||% "planned"
    )
    for (entity in entity_names) {
      row[[paste0("entity_", entity)]] <- entry$entities[[entity]] %||% NA_character_
    }
    row
  })
  data <- data.frame(matrix(nrow = n_rows, ncol = 0L))
  for (name in names(definitions)) {
    if (name %in% names(extras$values)) {
      data[[name]] <- extras$values[[name]]
    } else if (identical(definitions[[name]]$type, "json")) {
      data[[name]] <- I(lapply(rows, `[[`, name))
    } else {
      data[[name]] <- unlist(lapply(rows, `[[`, name), use.names = FALSE)
    }
  }
  validate_lna_manifest(structure(
    list(schema = .manifest_schema(definitions), data = data),
    class = c("lna_manifest", "list")
  ))
}

.manifest_format_scalar <- function(value, definition, field) {
  if (identical(definition$type, "json")) {
    return(.manifest_canonical_json(value, field))
  }
  if (identical(definition$type, "boolean")) return(if (value) "true" else "false")
  if (identical(definition$type, "integer")) return(format(value, scientific = FALSE, trim = TRUE))
  if (identical(definition$type, "number")) return(sprintf("%.17g", value))
  enc2utf8(value)
}

.manifest_to_tsv_data <- function(manifest) {
  data <- manifest$data
  schema <- manifest$schema
  out <- as.data.frame(
    stats::setNames(replicate(ncol(data), character(nrow(data)), simplify = FALSE), names(data)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (name in names(data)) {
    definition <- schema$columns[[name]]
    for (i in seq_len(nrow(data))) {
      missing <- .manifest_value_is_missing(data[[name]], i, definition$type)
      out[[name]][[i]] <- if (missing) {
        schema$null_encoding
      } else {
        value <- if (identical(definition$type, "json")) data[[name]][[i]] else data[[name]][[i]]
        .manifest_format_scalar(value, definition, paste0(name, "[", i, "]"))
      }
    }
  }
  out
}

.manifest_resolve_paths <- function(path) {
  path <- .manifest_scalar_string(path, "path")
  if (dir.exists(path) || !grepl("[.]tsv$|[.]json$", basename(path))) {
    root <- path
  } else {
    if (!basename(path) %in% c("lna_manifest.tsv", "lna_manifest.json")) {
      .manifest_abort(
        "manifest filename must be lna_manifest.tsv or lna_manifest.json", "path"
      )
    }
    root <- dirname(path)
  }
  list(
    root = root,
    tsv = file.path(root, "lna_manifest.tsv"),
    json = file.path(root, "lna_manifest.json")
  )
}

#' Write an LNA dataset manifest pair
#'
#' @param manifest A validated `lna_manifest`.
#' @param path Archive root, `lna_manifest.tsv`, or `lna_manifest.json` path.
#' @param overwrite Whether an existing manifest pair may be replaced.
#' @return Paths to the TSV and JSON files, invisibly.
#' @export
write_lna_manifest <- function(manifest, path, overwrite = FALSE) {
  manifest <- validate_lna_manifest(manifest)
  paths <- .manifest_resolve_paths(path)
  if (!dir.exists(paths$root)) {
    if (!dir.create(paths$root, recursive = TRUE, showWarnings = FALSE)) {
      .manifest_abort("could not create manifest directory", "path")
    }
  }
  existing <- c(paths$tsv, paths$json)[file.exists(c(paths$tsv, paths$json))]
  if (length(existing) && !isTRUE(overwrite)) {
    .manifest_abort(
      paste0("manifest file already exists: ", paste(basename(existing), collapse = ", ")),
      "path"
    )
  }
  tsv <- .manifest_to_tsv_data(manifest)
  utils::write.table(
    tsv, file = paths$tsv, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = manifest$schema$null_encoding, fileEncoding = "UTF-8"
  )
  jsonlite::write_json(
    manifest$schema, paths$json, auto_unbox = TRUE, null = "null",
    digits = NA, pretty = TRUE
  )
  invisible(paths[c("tsv", "json")])
}

.manifest_parse_scalar <- function(text, definition, field, null_encoding) {
  if (identical(text, null_encoding)) {
    if (!isTRUE(definition$nullable)) {
      .manifest_abort(paste0(field, " may not use the null encoding"), field)
    }
    return(NULL)
  }
  type <- definition$type
  if (type == "json") {
    value <- tryCatch(
      jsonlite::fromJSON(text, simplifyVector = FALSE),
      error = function(e) .manifest_abort(paste0(field, " is invalid JSON"), field)
    )
    return(.manifest_validate_scalar(value, definition, field))
  }
  if (type == "integer") {
    if (!grepl("^(0|[1-9][0-9]*)$", text)) {
      .manifest_abort(paste0(field, " is not a canonical non-negative integer"), field)
    }
    value <- suppressWarnings(as.numeric(text))
    return(.manifest_validate_scalar(value, definition, field))
  }
  if (type == "number") {
    if (!grepl("^[+-]?(([0-9]+([.][0-9]*)?)|([.][0-9]+))([eE][+-]?[0-9]+)?$", text)) {
      .manifest_abort(paste0(field, " is not a number"), field)
    }
    return(.manifest_validate_scalar(as.numeric(text), definition, field))
  }
  if (type == "boolean") {
    if (!text %in% c("true", "false")) {
      .manifest_abort(paste0(field, " must be true or false"), field)
    }
    return(identical(text, "true"))
  }
  .manifest_validate_scalar(text, definition, field)
}

.manifest_parse_table <- function(raw, schema) {
  n_rows <- nrow(raw)
  data <- raw
  for (name in names(schema$columns)) {
    definition <- schema$columns[[name]]
    parsed <- lapply(seq_len(n_rows), function(i) {
      .manifest_parse_scalar(
        raw[[name]][[i]], definition, paste0(name, "[", i, "]"),
        schema$null_encoding
      )
    })
    if (definition$type == "json") {
      data[[name]] <- I(parsed)
    } else if (definition$type %in% c("integer", "number")) {
      data[[name]] <- vapply(parsed, function(x) x %||% NA_real_, numeric(1))
    } else if (definition$type == "boolean") {
      data[[name]] <- vapply(parsed, function(x) x %||% NA, logical(1))
    } else {
      data[[name]] <- vapply(parsed, function(x) x %||% NA_character_, character(1))
    }
  }
  data
}

#' Read and validate an LNA dataset manifest pair
#'
#' The JSON sidecar is parsed and validated before any TSV values are typed.
#' Unsupported major versions fail closed. Declared future optional columns are
#' retained in both the schema and table.
#'
#' @param path Archive root, `lna_manifest.tsv`, or `lna_manifest.json` path.
#' @return A validated `lna_manifest`.
#' @export
read_lna_manifest <- function(path) {
  paths <- .manifest_resolve_paths(path)
  missing <- c(paths$json, paths$tsv)[!file.exists(c(paths$json, paths$tsv))]
  if (length(missing)) {
    .manifest_abort(
      paste0("manifest pair is incomplete; missing: ", paste(basename(missing), collapse = ", ")),
      "path"
    )
  }
  schema <- tryCatch(
    jsonlite::read_json(paths$json, simplifyVector = FALSE),
    error = function(e) .manifest_abort("lna_manifest.json is invalid JSON", "schema")
  )
  migration <- .manifest_migrate_legacy_schema(schema)
  schema <- .manifest_validate_schema(migration$schema)
  raw <- tryCatch(
    utils::read.delim(
      paths$tsv, header = TRUE, sep = "\t", quote = "", comment.char = "",
      colClasses = "character", check.names = FALSE, na.strings = character(),
      stringsAsFactors = FALSE, fileEncoding = "UTF-8"
    ),
    error = function(e) .manifest_abort(
      paste0("could not read lna_manifest.tsv: ", conditionMessage(e)), "data"
    )
  )
  raw <- .manifest_migrate_legacy_table(
    raw, schema, migration$added_columns
  )
  if (!identical(names(raw), schema$column_order)) {
    .manifest_abort(
      "TSV columns must exactly match schema.column_order", "column_order"
    )
  }
  validate_lna_manifest(structure(
    list(schema = schema, data = .manifest_parse_table(raw, schema)),
    class = c("lna_manifest", "list")
  ))
}

#' @export
print.lna_manifest <- function(x, ...) {
  cat("<lna_manifest>", x$schema$manifest_schema_version, "\n")
  cat("  rows:", nrow(x$data), "columns:", ncol(x$data), "\n")
  cat("  status:", paste(names(table(x$data$status)), table(x$data$status), collapse = ", "), "\n")
  invisible(x)
}
