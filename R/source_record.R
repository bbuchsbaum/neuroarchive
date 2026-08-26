# Typed, bidser-independent ingestion records for LNA dataset planning.

.lna_source_record_schema_version <- "1.0.0"
.lna_source_checksum_pattern <- "^sha256:[0-9a-f]{64}$"

.source_record_abort <- function(message, field) {
  abort_lna(
    message,
    field = field,
    .subclass = "lna_error_source_record",
    location = paste0("lna_source_record:", field)
  )
}

.source_record_scalar_string <- function(x, field, allow_null = FALSE) {
  if (allow_null && is.null(x)) return(NULL)
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    .source_record_abort(
      paste0(field, " must be one non-empty string"),
      field
    )
  }
  enc2utf8(x)
}

.source_record_id <- function(x, field = "record_id", allow_null = FALSE) {
  x <- .source_record_scalar_string(x, field, allow_null = allow_null)
  if (is.null(x)) return(NULL)
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9._:-]*$", x)) {
    .source_record_abort(
      paste0(
        field,
        " must start with an alphanumeric character and contain only ",
        "letters, numbers, '.', '_', ':', or '-'"
      ),
      field
    )
  }
  x
}

.source_record_checksum <- function(x, field) {
  x <- tolower(.source_record_scalar_string(x, field))
  if (!grepl(.lna_source_checksum_pattern, x)) {
    .source_record_abort(
      paste0(field, " must have the form sha256:<64 hexadecimal digits>"),
      field
    )
  }
  x
}

.source_record_nullable_checksum <- function(x, field) {
  if (is.null(x)) return(NULL)
  .source_record_checksum(x, field)
}

.source_record_named_list <- function(x, field, allow_empty = TRUE) {
  if (!is.list(x) || is.data.frame(x)) {
    .source_record_abort(paste0(field, " must be a named list"), field)
  }
  if (!length(x) && allow_empty) return(x)
  nms <- names(x)
  if (is.null(nms) || anyNA(nms) || any(!nzchar(nms)) || anyDuplicated(nms)) {
    .source_record_abort(
      paste0(field, " must have unique, non-empty names"),
      field
    )
  }
  x
}

.source_record_json_value <- function(x, field) {
  if (is.null(x)) return(NULL)
  if (is.environment(x) || is.function(x) || isS4(x) ||
      is.raw(x) || is.complex(x) || is.pairlist(x) || is.expression(x)) {
    .source_record_abort(
      paste0(field, " contains a value that cannot be represented in JSON"),
      field
    )
  }
  if (is.list(x)) {
    if (is.data.frame(x)) {
      .source_record_abort(
        paste0(field, " must use JSON objects/arrays, not data frames"),
        field
      )
    }
    nms <- names(x)
    if (!is.null(nms)) {
      named <- nzchar(nms)
      if (anyNA(nms) || (any(named) && !all(named)) || anyDuplicated(nms[named])) {
        .source_record_abort(
          paste0(field, " cannot mix named and unnamed JSON members"),
          field
        )
      }
    }
    for (i in seq_along(x)) {
      label <- if (!is.null(nms) && nzchar(nms[[i]])) nms[[i]] else i
      .source_record_json_value(x[[i]], paste0(field, ".", label))
    }
    return(x)
  }
  if (!is.atomic(x) || length(attributes(x)) || anyNA(x)) {
    .source_record_abort(
      paste0(field, " must contain JSON scalars/arrays with no NA values"),
      field
    )
  }
  if (is.numeric(x) && any(!is.finite(x))) {
    .source_record_abort(
      paste0(field, " must not contain NaN or infinite values"),
      field
    )
  }
  if (!is.character(x) && !is.logical(x) && !is.numeric(x)) {
    .source_record_abort(
      paste0(field, " contains an unsupported atomic type"),
      field
    )
  }
  x
}

.source_record_entities <- function(x) {
  if (is.atomic(x) && !is.null(names(x))) x <- as.list(x)
  x <- .source_record_named_list(x, "entities")
  if (!length(x)) return(x)

  nms <- names(x)
  if (any(!grepl("^[a-z][a-z0-9]*$", nms))) {
    .source_record_abort(
      "entity names must use lower-case ASCII letters and digits",
      "entities"
    )
  }
  for (nm in nms) {
    value <- .source_record_scalar_string(x[[nm]], paste0("entities.", nm))
    if (grepl("[/\\\\_]", value) || grepl("[[:cntrl:]]", value)) {
      .source_record_abort(
        paste0("entities.", nm, " contains a path or entity separator"),
        paste0("entities.", nm)
      )
    }
    x[[nm]] <- value
  }
  x[sort(names(x))]
}

.source_record_source <- function(x) {
  x <- .source_record_named_list(x, "source", allow_empty = FALSE)
  allowed <- c("path", "uri", "checksum")
  unknown <- setdiff(names(x), allowed)
  if (length(unknown)) {
    .source_record_abort(
      paste0("source has unknown field(s): ", paste(unknown, collapse = ", ")),
      "source"
    )
  }
  path <- .source_record_scalar_string(x$path, "source.path", allow_null = TRUE)
  uri <- .source_record_scalar_string(x$uri, "source.uri", allow_null = TRUE)
  if (is.null(path) && is.null(uri)) {
    .source_record_abort(
      "source requires at least one of path or uri",
      "source"
    )
  }
  if (!is.null(uri) && !grepl("^[A-Za-z][A-Za-z0-9+.-]*:", uri)) {
    .source_record_abort("source.uri must be an absolute URI", "source.uri")
  }
  list(
    path = path,
    uri = uri,
    checksum = .source_record_checksum(x$checksum, "source.checksum")
  )
}

.source_record_bids <- function(x) {
  x <- .source_record_named_list(x, "bids", allow_empty = FALSE)
  allowed <- c("datatype", "suffix", "extension", "scope", "pipeline")
  unknown <- setdiff(names(x), allowed)
  if (length(unknown)) {
    .source_record_abort(
      paste0("bids has unknown field(s): ", paste(unknown, collapse = ", ")),
      "bids"
    )
  }
  datatype <- .source_record_scalar_string(x$datatype, "bids.datatype")
  suffix <- .source_record_scalar_string(x$suffix, "bids.suffix")
  extension <- .source_record_scalar_string(x$extension, "bids.extension")
  scope <- .source_record_scalar_string(x$scope, "bids.scope")
  pipeline <- .source_record_scalar_string(
    x$pipeline, "bids.pipeline", allow_null = TRUE
  )
  if (!scope %in% c("raw", "derivative")) {
    .source_record_abort(
      "bids.scope must be either 'raw' or 'derivative'",
      "bids.scope"
    )
  }
  if (!startsWith(extension, ".") || grepl("[/\\\\]", extension)) {
    .source_record_abort(
      "bids.extension must start with '.' and contain no path separators",
      "bids.extension"
    )
  }
  list(
    datatype = datatype,
    suffix = suffix,
    extension = extension,
    scope = scope,
    pipeline = pipeline
  )
}

.source_record_affine <- function(x) {
  if (is.list(x)) {
    if (length(x) == 4L && all(vapply(x, length, integer(1)) == 4L)) {
      rows <- lapply(x, function(row) as.double(unlist(row, use.names = FALSE)))
      x <- matrix(unlist(rows, use.names = FALSE), nrow = 4L, byrow = TRUE)
    } else {
      x <- as.double(unlist(x, recursive = TRUE, use.names = FALSE))
    }
  }
  tryCatch(
    .lna_normalize_affine(x),
    lna_error_spatial_identity = function(e) {
      .source_record_abort(conditionMessage(e), "geometry.affine")
    }
  )
}

.source_record_dims <- function(x) {
  if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
  tryCatch(
    .lna_normalize_grid_dims(x),
    lna_error_spatial_identity = function(e) {
      .source_record_abort(conditionMessage(e), "geometry.dims")
    }
  )
}

.source_record_positive_integer <- function(x, field, allow_null = TRUE) {
  if (allow_null && is.null(x)) return(NULL)
  if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
  value <- suppressWarnings(as.integer(x))
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != value || value < 1L) {
    .source_record_abort(paste0(field, " must be a positive integer"), field)
  }
  value
}

.source_record_geometry <- function(x) {
  x <- .source_record_named_list(x, "geometry", allow_empty = FALSE)
  allowed <- c(
    "dims", "affine", "units", "voxel_ordering", "n_volumes",
    "grid_id", "support_id"
  )
  unknown <- setdiff(names(x), allowed)
  if (length(unknown)) {
    .source_record_abort(
      paste0("geometry has unknown field(s): ", paste(unknown, collapse = ", ")),
      "geometry"
    )
  }
  dims <- .source_record_dims(x$dims)
  affine <- .source_record_affine(x$affine)
  units <- .source_record_scalar_string(x$units, "geometry.units")
  voxel_ordering <- .source_record_scalar_string(
    x$voxel_ordering, "geometry.voxel_ordering"
  )
  grid_id <- .lna_grid_id(dims, affine, units, voxel_ordering)
  if (!is.null(x$grid_id)) {
    supplied <- .source_record_checksum(x$grid_id, "geometry.grid_id")
    if (!identical(supplied, grid_id)) {
      .source_record_abort(
        "geometry.grid_id does not match dims, affine, units, and voxel_ordering",
        "geometry.grid_id"
      )
    }
  }
  support_id <- .source_record_nullable_checksum(
    x$support_id, "geometry.support_id"
  )
  list(
    dims = dims,
    affine = affine,
    units = units,
    voxel_ordering = voxel_ordering,
    n_volumes = .source_record_positive_integer(
      x$n_volumes, "geometry.n_volumes"
    ),
    grid_id = grid_id,
    support_id = support_id
  )
}

.source_record_association <- function(x, index, record_id) {
  field <- paste0("associations[[", index, "]]")
  x <- .source_record_named_list(x, field, allow_empty = FALSE)
  allowed <- c(
    "association_id", "role", "target_record_id", "path", "uri",
    "checksum", "metadata"
  )
  unknown <- setdiff(names(x), allowed)
  if (length(unknown)) {
    .source_record_abort(
      paste0(field, " has unknown field(s): ", paste(unknown, collapse = ", ")),
      field
    )
  }
  association_id <- .source_record_id(
    x$association_id, paste0(field, ".association_id")
  )
  role <- .source_record_scalar_string(x$role, paste0(field, ".role"))
  target_record_id <- .source_record_id(
    x$target_record_id,
    paste0(field, ".target_record_id"),
    allow_null = TRUE
  )
  path <- .source_record_scalar_string(
    x$path, paste0(field, ".path"), allow_null = TRUE
  )
  uri <- .source_record_scalar_string(
    x$uri, paste0(field, ".uri"), allow_null = TRUE
  )
  if (is.null(target_record_id) && is.null(path) && is.null(uri)) {
    .source_record_abort(
      paste0(field, " requires target_record_id, path, or uri"),
      field
    )
  }
  if (!is.null(target_record_id) && identical(target_record_id, record_id)) {
    .source_record_abort(
      paste0(field, " cannot target its containing record"),
      paste0(field, ".target_record_id")
    )
  }
  if (!is.null(uri) && !grepl("^[A-Za-z][A-Za-z0-9+.-]*:", uri)) {
    .source_record_abort(
      paste0(field, ".uri must be an absolute URI"),
      paste0(field, ".uri")
    )
  }
  if ((!is.null(path) || !is.null(uri)) && is.null(x$checksum)) {
    .source_record_abort(
      paste0(field, " requires checksum for a path/URI association"),
      paste0(field, ".checksum")
    )
  }
  metadata <- x$metadata %||% list()
  .source_record_named_list(metadata, paste0(field, ".metadata"))
  .source_record_json_value(metadata, paste0(field, ".metadata"))
  list(
    association_id = association_id,
    role = role,
    target_record_id = target_record_id,
    path = path,
    uri = uri,
    checksum = .source_record_nullable_checksum(
      x$checksum, paste0(field, ".checksum")
    ),
    metadata = metadata
  )
}

.source_record_associations <- function(x, record_id) {
  if (!is.list(x) || is.data.frame(x)) {
    .source_record_abort("associations must be a list", "associations")
  }
  if (!length(x)) return(list())
  out <- lapply(
    seq_along(x),
    function(i) .source_record_association(x[[i]], i, record_id)
  )
  ids <- vapply(out, `[[`, character(1), "association_id")
  if (anyDuplicated(ids)) {
    duplicate <- unique(ids[duplicated(ids)])
    .source_record_abort(
      paste0("association_id is duplicated: ", paste(duplicate, collapse = ", ")),
      "associations.association_id"
    )
  }
  out[order(ids)]
}

.source_record_reference <- function(x, record_id) {
  if (is.null(x)) {
    return(unclass(lna_reference_assertion(record_id = record_id)))
  }
  unclass(tryCatch(
    .as_lna_reference_assertion(x, record_id = record_id),
    lna_error_reference_assertion = function(e) {
      abort_lna(
        conditionMessage(e),
        field = paste0("reference_assertion.", e$field %||% "unknown"),
        .subclass = "lna_error_source_record",
        location = "lna_source_record:reference_assertion",
        parent = e
      )
    }
  ))
}

.source_record_from_raw <- function(x) {
  x <- .source_record_named_list(x, "source_record", allow_empty = FALSE)
  required <- c(
    "schema_version", "record_id", "source", "bids", "entities",
    "resolved_metadata", "associations", "geometry", "reference_assertion"
  )
  unknown <- setdiff(names(x), required)
  missing <- setdiff(required, names(x))
  if (length(unknown)) {
    .source_record_abort(
      paste0("source_record has unknown field(s): ", paste(unknown, collapse = ", ")),
      "source_record"
    )
  }
  if (length(missing)) {
    .source_record_abort(
      paste0("source_record is missing field(s): ", paste(missing, collapse = ", ")),
      "source_record"
    )
  }
  schema_version <- .source_record_scalar_string(
    x$schema_version, "schema_version"
  )
  if (!identical(schema_version, .lna_source_record_schema_version)) {
    .source_record_abort(
      paste0(
        "unsupported source-record schema_version: ", schema_version,
        "; expected ", .lna_source_record_schema_version
      ),
      "schema_version"
    )
  }
  record_id <- .source_record_id(x$record_id)
  metadata <- .source_record_named_list(
    x$resolved_metadata, "resolved_metadata"
  )
  .source_record_json_value(metadata, "resolved_metadata")

  structure(
    list(
      schema_version = schema_version,
      record_id = record_id,
      source = .source_record_source(x$source),
      bids = .source_record_bids(x$bids),
      entities = .source_record_entities(x$entities),
      resolved_metadata = metadata,
      associations = .source_record_associations(x$associations, record_id),
      geometry = .source_record_geometry(x$geometry),
      reference_assertion = .source_record_reference(
        x$reference_assertion, record_id
      )
    ),
    class = "lna_source_record"
  )
}

#' Construct a typed LNA dataset source record
#'
#' A source record is the bidser-independent boundary between dataset discovery
#' and LNA planning. Entity values are strings by contract, so indices such as
#' `run = "01"` and `echo = "02"` retain their original spelling. The record is
#' normalized into a JSON-safe, copy-on-modify R value and has no mutating API.
#'
#' @param record_id Stable, dataset-local record identifier.
#' @param source_path Optional source path, usually relative to the source
#'   dataset root.
#' @param source_uri Optional absolute source URI. At least one of
#'   `source_path` and `source_uri` is required.
#' @param source_checksum SHA-256 checksum in `sha256:<64 hex>` form.
#' @param entities Named list of BIDS entity labels. Values must be strings.
#' @param datatype BIDS datatype directory, such as `"func"`.
#' @param suffix BIDS suffix, such as `"bold"`.
#' @param extension Original imaging extension, such as `".nii.gz"`.
#' @param scope Either `"raw"` or `"derivative"`.
#' @param pipeline Optional derivative pipeline identifier.
#' @param resolved_metadata Named, JSON-safe list after metadata inheritance has
#'   been resolved by an adapter.
#' @param associations List of typed association records. Each entry requires
#'   `association_id`, `role`, and either `target_record_id` or a checksummed
#'   `path`/`uri` locator.
#' @param geometry Named list containing `dims`, `affine`, `units`,
#'   `voxel_ordering`, and optional `n_volumes`, `grid_id`, or `support_id`.
#' @param reference_assertion Optional reference-space provenance assertion.
#'   Missing assertions become a record-unique unknown-native assertion.
#' @return An object of class `lna_source_record`.
#' @export
lna_source_record <- function(
    record_id,
    source_path = NULL,
    source_uri = NULL,
    source_checksum,
    entities = list(),
    datatype,
    suffix,
    extension,
    scope = "raw",
    pipeline = NULL,
    resolved_metadata = list(),
    associations = list(),
    geometry,
    reference_assertion = NULL) {
  raw <- list(
    schema_version = .lna_source_record_schema_version,
    record_id = record_id,
    source = list(
      path = source_path,
      uri = source_uri,
      checksum = source_checksum
    ),
    bids = list(
      datatype = datatype,
      suffix = suffix,
      extension = extension,
      scope = scope,
      pipeline = pipeline
    ),
    entities = entities,
    resolved_metadata = resolved_metadata,
    associations = associations,
    geometry = geometry,
    reference_assertion = reference_assertion
  )
  .source_record_from_raw(raw)
}

#' Coerce a decoded list to an LNA source record
#'
#' This is the inverse of `as.list()` for source records and is intended for
#' records decoded from JSON with `simplifyVector = FALSE`.
#'
#' @param x A named list following source-record schema version 1.0.0.
#' @return A validated `lna_source_record`.
#' @export
as_lna_source_record <- function(x) {
  if (inherits(x, "lna_source_record")) {
    x <- unclass(x)
  }
  .source_record_from_raw(x)
}

#' @export
as.list.lna_source_record <- function(x, ...) {
  unclass(x)
}

#' Validate a collection of LNA source records
#'
#' In addition to per-record validation, this checks record-ID uniqueness,
#' one-source-to-one-record identity, and closure of record-to-record
#' associations. It returns normalized copies and never modifies its inputs.
#'
#' @param records One `lna_source_record` or a list of source records/decoded
#'   record lists.
#' @param require_association_closure If `TRUE`, every `target_record_id` must
#'   occur in `records`.
#' @return Invisibly, a normalized source record or list of source records,
#'   matching the shape of `records`.
#' @export
validate_lna_source_records <- function(
    records,
    require_association_closure = TRUE) {
  single <- inherits(records, "lna_source_record") ||
    (is.list(records) && setequal(names(records), c(
      "schema_version", "record_id", "source", "bids", "entities",
      "resolved_metadata", "associations", "geometry", "reference_assertion"
    )))
  input <- if (single) list(records) else records
  if (!is.list(input) || !length(input)) {
    .source_record_abort(
      "records must contain at least one source record",
      "records"
    )
  }
  normalized <- lapply(input, as_lna_source_record)
  ids <- vapply(normalized, `[[`, character(1), "record_id")
  if (anyDuplicated(ids)) {
    duplicate <- unique(ids[duplicated(ids)])
    .source_record_abort(
      paste0("record_id is duplicated: ", paste(duplicate, collapse = ", ")),
      "record_id"
    )
  }

  source_keys <- vapply(normalized, function(record) {
    source <- record$source
    paste(source$path %||% "", source$uri %||% "", source$checksum, sep = "\r")
  }, character(1))
  if (anyDuplicated(source_keys)) {
    duplicated_ids <- ids[
      duplicated(source_keys) | duplicated(source_keys, fromLast = TRUE)
    ]
    .source_record_abort(
      paste0(
        "one source locator/checksum is assigned to multiple records: ",
        paste(duplicated_ids, collapse = ", ")
      ),
      "source"
    )
  }

  if (isTRUE(require_association_closure)) {
    targets <- unlist(lapply(normalized, function(record) {
      vapply(
        record$associations,
        function(association) association$target_record_id %||% NA_character_,
        character(1)
      )
    }), use.names = FALSE)
    targets <- unique(targets[!is.na(targets)])
    missing <- setdiff(targets, ids)
    if (length(missing)) {
      .source_record_abort(
        paste0(
          "record association target(s) are absent: ",
          paste(missing, collapse = ", ")
        ),
        "associations.target_record_id"
      )
    }
  }

  invisible(if (single) normalized[[1L]] else normalized)
}
