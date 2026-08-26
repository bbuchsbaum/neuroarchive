# Optional BIDS discovery adapter. The dataset engine itself remains bidser-free.

.lna_bidser_min_version <- "0.5.1"
.lna_bidser_required_exports <- c(
  "bids_project", "query_files", "bids_entities", "get_metadata",
  "event_files", "read_events", "derivative_pipelines", "n_volumes"
)

.bidser_adapter_abort <- function(message, field, ...) {
  abort_lna(
    message,
    ...,
    field = field,
    .subclass = "lna_error_bidser_adapter",
    location = paste0("lna_source_records_from_bidser:", field)
  )
}

.bidser_adapter_validate_contract <- function(version, exports) {
  version <- base::package_version(as.character(version))
  if (version < base::package_version(.lna_bidser_min_version)) {
    .bidser_adapter_abort(
      paste0(
        "installed bidser ", version, " is incompatible; bidser >= ",
        .lna_bidser_min_version,
        " is required for compound-extension and metadata-provenance support"
      ),
      "bidser_version",
      installed_version = as.character(version),
      required_version = .lna_bidser_min_version
    )
  }
  missing <- setdiff(.lna_bidser_required_exports, exports)
  if (length(missing)) {
    .bidser_adapter_abort(
      paste0(
        "bidser ", version, " is missing required public API(s): ",
        paste(missing, collapse = ", ")
      ),
      "bidser_api"
    )
  }
  as.character(version)
}

.bidser_adapter_api <- function() {
  if (!requireNamespace("bidser", quietly = TRUE)) {
    .bidser_adapter_abort(
      paste0(
        "bidser >= ", .lna_bidser_min_version,
        " is required only for BIDS discovery. Install the released package ",
        "or a reviewed 0.5.1+ source checkout, then retry."
      ),
      "bidser_version"
    )
  }
  version <- utils::packageVersion("bidser")
  exports <- getNamespaceExports("bidser")
  version <- .bidser_adapter_validate_contract(version, exports)
  api <- stats::setNames(lapply(.lna_bidser_required_exports, function(name) {
    getExportedValue("bidser", name)
  }), .lna_bidser_required_exports)
  api$version <- version
  api
}

.bidser_adapter_json <- function(x, field) {
  normalized <- tryCatch(
    jsonlite::fromJSON(
      jsonlite::toJSON(
        x, auto_unbox = TRUE, null = "null", na = "null",
        digits = NA, dataframe = "rows", matrix = "rowmajor"
      ),
      simplifyVector = FALSE
    ),
    error = function(e) {
      .bidser_adapter_abort(
        paste0(field, " is not JSON-safe: ", conditionMessage(e)), field,
        parent = e
      )
    }
  )
  normalized %||% list()
}

.bidser_adapter_exact_regex <- function(value) {
  if (is.null(value)) return(".*")
  escaped <- gsub("([][{}()+*^$|\\?.])", "\\\\\\1", value, perl = TRUE)
  paste0("^", escaped, "$")
}

.bidser_adapter_relative_path <- function(path, field) {
  path <- gsub("\\\\", "/", as.character(path))
  tryCatch(
    .manifest_validate_path(path, field),
    lna_error_manifest = function(e) {
      .bidser_adapter_abort(conditionMessage(e), field, parent = e)
    }
  )
}

.bidser_adapter_project_root <- function(absolute_path, relative_path) {
  root <- normalizePath(absolute_path, winslash = "/", mustWork = TRUE)
  pieces <- strsplit(relative_path, "/", fixed = TRUE)[[1L]]
  for (i in seq_along(pieces)) root <- dirname(root)
  root
}

.bidser_adapter_match_absolute <- function(relative_path, absolute_paths) {
  normalized <- gsub("\\\\", "/", as.character(absolute_paths))
  suffix <- paste0("/", relative_path)
  hits <- normalized[endsWith(normalized, suffix)]
  if (length(hits) != 1L) {
    .bidser_adapter_abort(
      paste0(
        "expected one absolute source for ", relative_path,
        "; found ", length(hits)
      ),
      "source_path",
      relative_path = relative_path
    )
  }
  normalizePath(hits, winslash = "/", mustWork = TRUE)
}

.bidser_adapter_entities <- function(api, relative_path) {
  row <- api$bids_entities(
    relative_path, include_path = FALSE, coerce = FALSE
  )
  if (!is.data.frame(row) || nrow(row) != 1L) {
    .bidser_adapter_abort(
      "bidser::bids_entities() did not return exactly one row",
      "entities"
    )
  }
  value <- function(name) {
    if (!name %in% names(row) || is.na(row[[name]][[1L]])) return(NULL)
    as.character(row[[name]][[1L]])
  }
  aliases <- list(
    sub = value("subid"), ses = value("session"),
    acq = value("acq") %||% value("acquisition"),
    ce = value("ce") %||% value("contrast"),
    rec = value("rec") %||% value("reconstruction")
  )
  keys <- c(
    "task", "dir", "run", "echo", "part", "chunk", "space", "res",
    "den", "desc", "hemi", "label", "from", "to", "mode",
    "recording", "proc", "split"
  )
  entities <- c(aliases, stats::setNames(lapply(keys, value), keys))
  entities[!vapply(entities, is.null, logical(1))]
}

.bidser_adapter_units <- function(header, relative_path) {
  code <- header$raw$xyzt_units %||% NA_integer_
  code <- suppressWarnings(as.integer(code[[1L]]))
  spatial_code <- if (is.na(code)) NA_integer_ else bitwAnd(code, 7L)
  units <- switch(
    as.character(spatial_code),
    "1" = "m", "2" = "mm", "3" = "um", NULL
  )
  if (is.null(units)) {
    .bidser_adapter_abort(
      paste0(
        "NIfTI spatial units are absent or unsupported for ", relative_path
      ),
      "geometry.units",
      relative_path = relative_path
    )
  }
  units
}

.bidser_adapter_geometry <- function(api, absolute_path, relative_path) {
  header <- tryCatch(
    neuroim2::header(absolute_path),
    error = function(e) {
      .bidser_adapter_abort(
        paste0("could not read NIfTI header: ", conditionMessage(e)),
        "geometry", relative_path = relative_path, parent = e
      )
    }
  )
  dims <- suppressWarnings(as.integer(header$dim))
  if (length(dims) < 3L || anyNA(dims[1:3]) || any(dims[1:3] < 1L)) {
    .bidser_adapter_abort(
      paste0("NIfTI header has invalid spatial dimensions: ", relative_path),
      "geometry.dims"
    )
  }
  affine <- header$trans
  if (is.null(affine)) {
    .bidser_adapter_abort(
      paste0("NIfTI header has no affine: ", relative_path),
      "geometry.affine"
    )
  }
  n_volumes <- tryCatch(
    unname(api$n_volumes(absolute_path))[[1L]],
    error = function(e) {
      .bidser_adapter_abort(
        paste0("could not determine volume count: ", conditionMessage(e)),
        "geometry.n_volumes", relative_path = relative_path, parent = e
      )
    }
  )
  n_volumes <- suppressWarnings(as.integer(n_volumes))
  if (length(dims) >= 4L && !identical(n_volumes, dims[[4L]])) {
    .bidser_adapter_abort(
      paste0("bidser volume count disagrees with NIfTI header: ", relative_path),
      "geometry.n_volumes",
      bidser_n_volumes = n_volumes,
      header_n_volumes = dims[[4L]]
    )
  }
  list(
    dims = dims[1:3],
    affine = affine,
    units = .bidser_adapter_units(header, relative_path),
    voxel_ordering = "first-axis-fastest",
    n_volumes = n_volumes
  )
}

.bidser_adapter_association_id <- function(role, path) {
  paste0(
    "assoc:", role, ":",
    substr(digest::digest(path, algo = "sha256", serialize = FALSE), 1L, 20L)
  )
}

.bidser_adapter_file_association <- function(
    role, relative_path, project_root, metadata = list()) {
  relative_path <- .bidser_adapter_relative_path(
    relative_path, paste0("associations.", role, ".path")
  )
  absolute_path <- file.path(
    project_root,
    do.call(file.path, as.list(strsplit(relative_path, "/", fixed = TRUE)[[1L]]))
  )
  if (!file.exists(absolute_path)) {
    .bidser_adapter_abort(
      paste0("associated BIDS file does not exist: ", relative_path),
      "associations.path"
    )
  }
  list(
    association_id = .bidser_adapter_association_id(role, relative_path),
    role = role,
    path = relative_path,
    checksum = .journal_file_checksum(absolute_path),
    metadata = metadata
  )
}

.bidser_adapter_metadata <- function(api, project, relative_path, project_root) {
  result <- tryCatch(
    api$get_metadata(
      project, relative_path, inherit = TRUE, scope = "auto",
      provenance = TRUE
    ),
    error = function(e) {
      .bidser_adapter_abort(
        paste0("bidser metadata resolution failed: ", conditionMessage(e)),
        "resolved_metadata", relative_path = relative_path, parent = e
      )
    }
  )
  if (!is.list(result) || !setequal(names(result), c("metadata", "sources")) ||
      !is.list(result$metadata) || !is.list(result$sources)) {
    .bidser_adapter_abort(
      "bidser metadata provenance response is incompatible with 0.5.1",
      "resolved_metadata"
    )
  }
  metadata <- .bidser_adapter_json(result$metadata, "resolved_metadata")
  associations <- lapply(result$sources, function(source) {
    required <- c("path", "inheritance_level", "precedence", "fields")
    if (!is.list(source) || length(setdiff(required, names(source)))) {
      .bidser_adapter_abort(
        "bidser returned a malformed metadata source", "metadata_provenance"
      )
    }
    .bidser_adapter_file_association(
      "metadata-source", source$path, project_root,
      metadata = list(
        inheritance_level = as.character(source$inheritance_level),
        precedence = as.integer(source$precedence),
        fields = sort(as.character(unlist(source$fields, use.names = FALSE))),
        media_type = "application/json"
      )
    )
  })
  list(metadata = metadata, associations = associations, sources = result$sources)
}

.bidser_adapter_event_associations <- function(
    api, project, entities, project_root) {
  event_paths <- tryCatch(
    api$event_files(
      project,
      subid = .bidser_adapter_exact_regex(entities$sub),
      task = .bidser_adapter_exact_regex(entities$task),
      run = .bidser_adapter_exact_regex(entities$run),
      session = .bidser_adapter_exact_regex(entities$ses),
      full_path = FALSE
    ),
    error = function(e) {
      .bidser_adapter_abort(
        paste0("bidser event discovery failed: ", conditionMessage(e)),
        "associations.events", parent = e
      )
    }
  )
  event_paths <- sort(unique(gsub("\\\\", "/", as.character(event_paths %||% character()))))
  if (!length(event_paths)) return(list())

  same_record <- vapply(event_paths, function(path) {
    candidate <- .bidser_adapter_entities(api, path)
    keys <- intersect(c("sub", "ses", "task", "run"), names(entities))
    all(vapply(keys, function(key) {
      identical(candidate[[key]], entities[[key]])
    }, logical(1)))
  }, logical(1))
  event_paths <- event_paths[same_record]
  out <- list()
  for (event_path in event_paths) {
    out[[length(out) + 1L]] <- .bidser_adapter_file_association(
      "events", event_path, project_root,
      metadata = list(media_type = "text/tab-separated-values")
    )
    json_path <- sub("[.]tsv$", ".json", event_path)
    if (!identical(json_path, event_path) &&
        file.exists(file.path(project_root, json_path))) {
      out[[length(out) + 1L]] <- .bidser_adapter_file_association(
        "events-json", json_path, project_root,
        metadata = list(media_type = "application/json")
      )
    }
  }
  out
}

.bidser_adapter_reference <- function(
    record_id, entities, metadata, metadata_sources) {
  space <- entities$space
  spatial_reference <- metadata$SpatialReference
  if (!is.null(space)) {
    return(lna_reference_assertion(
      record_id = record_id,
      identifier = space,
      kind = "template",
      origin = "bids-space-entity",
      evidence = list(space_entity = space)
    ))
  }
  if (!is.null(spatial_reference)) {
    if (is.list(spatial_reference)) {
      spatial_reference <- unlist(spatial_reference, recursive = TRUE, use.names = FALSE)
    }
    if (!is.character(spatial_reference) || length(spatial_reference) != 1L ||
        is.na(spatial_reference) || !nzchar(spatial_reference)) {
      .bidser_adapter_abort(
        "resolved SpatialReference is not one unambiguous identifier",
        "reference_assertion"
      )
    }
    candidates <- metadata_sources[vapply(metadata_sources, function(source) {
      "SpatialReference" %in% as.character(unlist(source$fields, use.names = FALSE))
    }, logical(1))]
    source_path <- if (length(candidates)) {
      candidates[[which.max(vapply(
        candidates, function(source) as.integer(source$precedence), integer(1)
      ))]]$path
    } else NULL
    evidence <- list(SpatialReference = spatial_reference)
    if (!is.null(source_path)) evidence$sidecar_path <- source_path
    return(lna_reference_assertion(
      record_id = record_id,
      identifier = spatial_reference,
      kind = "named-reference",
      origin = "spatial-reference-sidecar",
      evidence = evidence
    ))
  }
  NULL
}

.lna_source_records_from_bidser <- function(
    project,
    api,
    scope = c("all", "raw", "derivatives"),
    pipeline = NULL,
    include_events = TRUE) {
  scope <- match.arg(scope)
  if (!is.logical(include_events) || length(include_events) != 1L ||
      is.na(include_events)) {
    .bidser_adapter_abort("include_events must be TRUE or FALSE", "include_events")
  }
  query_args <- list(
    x = project,
    regex = ".*",
    full_path = FALSE,
    match_mode = "exact",
    scope = scope,
    pipeline = pipeline,
    return = "tibble",
    use_index = "auto",
    kind = "bold",
    datatype = "func",
    extension = c(".nii.gz", ".nii")
  )
  rows <- tryCatch(
    do.call(api$query_files, query_args),
    error = function(e) {
      .bidser_adapter_abort(
        paste0("bidser imaging query failed: ", conditionMessage(e)),
        "query", parent = e
      )
    }
  )
  if (!is.data.frame(rows)) {
    .bidser_adapter_abort("bidser query did not return a tibble", "query")
  }
  if (!nrow(rows)) return(list())
  required_columns <- c(
    "path", "scope", "pipeline", "extension", "datatype", "kind"
  )
  missing <- setdiff(required_columns, names(rows))
  if (length(missing)) {
    .bidser_adapter_abort(
      paste0("bidser query rows lack: ", paste(missing, collapse = ", ")),
      "query"
    )
  }
  rows <- rows[order(rows$path), , drop = FALSE]
  query_args$full_path <- TRUE
  query_args$return <- "paths"
  absolute_paths <- do.call(api$query_files, query_args)

  records <- lapply(seq_len(nrow(rows)), function(i) {
    relative_path <- .bidser_adapter_relative_path(
      rows$path[[i]], paste0("query.path[[", i, "]]" )
    )
    absolute_path <- .bidser_adapter_match_absolute(
      relative_path, absolute_paths
    )
    project_root <- .bidser_adapter_project_root(absolute_path, relative_path)
    entities <- .bidser_adapter_entities(api, relative_path)
    if (is.null(entities$sub)) {
      .bidser_adapter_abort(
        paste0("functional image lacks a subject entity: ", relative_path),
        "entities.sub"
      )
    }
    record_id <- paste0(
      "bids:",
      substr(digest::digest(relative_path, algo = "sha256", serialize = FALSE), 1L, 24L)
    )
    metadata_result <- .bidser_adapter_metadata(
      api, project, relative_path, project_root
    )
    event_associations <- if (isTRUE(include_events)) {
      .bidser_adapter_event_associations(
        api, project, entities, project_root
      )
    } else list()
    source_scope <- as.character(rows$scope[[i]])
    source_scope <- if (identical(source_scope, "derivatives")) {
      "derivative"
    } else source_scope
    pipeline_value <- rows$pipeline[[i]]
    if (length(pipeline_value) == 0L || is.na(pipeline_value) ||
        !nzchar(as.character(pipeline_value))) {
      pipeline_value <- NULL
    } else pipeline_value <- as.character(pipeline_value)
    suffix <- as.character(rows$kind[[i]])
    if (is.na(suffix) || !nzchar(suffix)) {
      .bidser_adapter_abort(
        paste0("functional image lacks a BIDS suffix: ", relative_path),
        "bids.suffix"
      )
    }

    lna_source_record(
      record_id = record_id,
      source_path = relative_path,
      source_checksum = .journal_file_checksum(absolute_path),
      entities = entities,
      datatype = as.character(rows$datatype[[i]]),
      suffix = suffix,
      extension = as.character(rows$extension[[i]]),
      scope = source_scope,
      pipeline = pipeline_value,
      resolved_metadata = metadata_result$metadata,
      associations = c(metadata_result$associations, event_associations),
      geometry = .bidser_adapter_geometry(api, absolute_path, relative_path),
      reference_assertion = .bidser_adapter_reference(
        record_id, entities, metadata_result$metadata,
        metadata_result$sources
      )
    )
  })
  validate_lna_source_records(records)
}

#' Convert a bidser project into typed LNA source records
#'
#' This optional adapter is the only neuroarchive layer that knows about
#' bidser. bidser owns BIDS discovery, exact entity parsing, inheritance,
#' derivative pipeline identity, events, and volume counts. neuroarchive maps
#' those public results into its own immutable source-record contract and reads
#' only NIfTI header geometry; it does not reimplement a BIDS layout.
#'
#' @param project A `bidser::bids_project()` object.
#' @param scope One of `"all"`, `"raw"`, or `"derivatives"`.
#' @param pipeline Optional derivative pipeline name.
#' @param include_events Whether matching event TSV/JSON files should be added
#'   as byte-preserved associations.
#' @return A validated list of `lna_source_record` objects, ordered by source
#'   path. No bidser object is retained in the records.
#' @export
lna_source_records_from_bidser <- function(
    project,
    scope = c("all", "raw", "derivatives"),
    pipeline = NULL,
    include_events = TRUE) {
  .lna_source_records_from_bidser(
    project = project,
    api = .bidser_adapter_api(),
    scope = scope,
    pipeline = pipeline,
    include_events = include_events
  )
}
