# Manifest-only dataset opening, namespaced views, and exact queries.

.query_abort <- function(message, field, ...) {
  abort_lna(
    message,
    ...,
    field = field,
    .subclass = "lna_error_dataset_query",
    location = paste0("query_lna:", field)
  )
}

.query_root <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    .query_abort("path must be one non-empty string", "path")
  }
  root <- if (identical(basename(path), "lna_dataset.json")) dirname(path) else path
  if (!dir.exists(root)) .query_abort("dataset root does not exist", "path")
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

.query_dataset <- function(x) {
  if (inherits(x, "lna_dataset_handle")) return(x)
  open_lna_dataset(x)
}

.query_view_columns <- function(data, namespace) {
  identity <- c("record_id", "lna_path")
  candidates <- switch(
    namespace,
    records = c(identity, "status"),
    entities = c(
      identity, "datatype", "suffix", "scope", "pipeline",
      grep("^entity_", names(data), value = TRUE)
    ),
    spatial = c(
      identity, "n_volumes", "reference_identifier", "reference_kind",
      "reference_origin", "reference_status", "reference_json", "grid_id",
      "support_id"
    ),
    codec = c(
      identity, "codec_profile_id", "loss_mode", "codec_profile_json"
    ),
    metadata = c(
      identity, "repetition_time_seconds", "metadata_json",
      "metadata_provenance_json", "associations_json"
    ),
    integrity = c(
      identity, "source_path", "source_uri", "source_checksum", "source_bytes",
      "source_extension", "lna_checksum", "lna_bytes", "status"
    ),
    all = names(data),
    .query_abort("unsupported manifest namespace", "namespace")
  )
  unique(candidates[candidates %in% names(data)])
}

#' Open an LNA exchange dataset without opening its HDF5 artifacts
#'
#' This validates `lna_dataset.json` and the manifest schema/table pair before
#' exposing any records. Artifact HDF5 files are not opened or inspected.
#'
#' @param path Archive root or `lna_dataset.json` path.
#' @return An immutable `lna_dataset_handle` backed by the validated manifest.
#' @export
open_lna_dataset <- function(path) {
  root <- .query_root(path)
  declaration <- tryCatch(
    read_lna_dataset_declaration(root),
    error = function(e) .query_abort(
      paste0("dataset declaration is invalid: ", conditionMessage(e)),
      "declaration", parent = e
    )
  )
  manifest <- tryCatch(
    read_lna_manifest(root),
    error = function(e) .query_abort(
      paste0("dataset manifest is invalid: ", conditionMessage(e)),
      "manifest", parent = e
    )
  )
  if (!identical(
    declaration$manifest$schema_version,
    manifest$schema$manifest_schema_version
  )) {
    .query_abort(
      "declaration and manifest schema versions disagree", "manifest"
    )
  }
  structure(list(
    root = root,
    declaration = declaration,
    manifest = manifest
  ), class = c("lna_dataset_handle", "list"))
}

#' Return a namespaced view of an LNA dataset manifest
#'
#' @param dataset Dataset path or object returned by [open_lna_dataset()].
#' @param namespace One of `"records"`, `"entities"`, `"spatial"`, `"codec"`,
#'   `"metadata"`, `"integrity"`, or `"all"`.
#' @return A deterministic data-frame view containing `record_id` and
#'   `lna_path`.
#' @export
lna_dataset_view <- function(
    dataset,
    namespace = c(
      "records", "entities", "spatial", "codec", "metadata", "integrity", "all"
    )) {
  namespace <- match.arg(namespace)
  dataset <- .query_dataset(dataset)
  data <- dataset$manifest$data
  columns <- .query_view_columns(data, namespace)
  order <- order(data$lna_path, data$record_id)
  out <- data[order, columns, drop = FALSE]
  rownames(out) <- NULL
  structure(out, class = c("lna_manifest_view", class(out)), namespace = namespace)
}

.query_filter_column <- function(name, columns) {
  aliases <- c(
    codec = "codec_profile_id",
    reference = "reference_identifier",
    grid = "grid_id",
    support = "support_id"
  )
  if (name %in% names(aliases)) name <- unname(aliases[[name]])
  if (name %in% columns) return(name)
  if (identical(name, "sub")) return("entity_sub")
  if (startsWith(name, "entity_")) return(name)
  if (grepl("^[a-z][a-z0-9]*$", name)) return(paste0("entity_", name))
  .query_abort(paste0("unknown query field: ", name), paste0("filters.", name))
}

.query_patterns <- function(x, field) {
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x) || !length(x) || anyNA(x) || any(!nzchar(x))) {
    .query_abort(
      paste0(field, " must contain one or more non-empty strings"), field
    )
  }
  enc2utf8(x)
}

.query_matches <- function(values, patterns, mode, field) {
  present <- !is.na(values)
  out <- rep(FALSE, length(values))
  if (!any(present)) return(out)
  if (identical(mode, "exact")) {
    out[present] <- values[present] %in% patterns
    return(out)
  }
  regex <- if (identical(mode, "glob")) {
    vapply(patterns, utils::glob2rx, character(1), trim.head = FALSE,
      trim.tail = FALSE)
  } else patterns
  matched <- tryCatch(
    suppressWarnings(vapply(values[present], function(value) {
      any(vapply(regex, grepl, logical(1), x = value, perl = TRUE))
    }, logical(1))),
    error = function(e) .query_abort(
      paste0("invalid ", mode, " pattern for ", field, ": ", conditionMessage(e)),
      field, parent = e
    )
  )
  out[present] <- matched
  out
}

#' Query an LNA dataset using validated manifest fields
#'
#' @param dataset Dataset path or open handle.
#' @param ... Named filters. BIDS entities use their unprefixed names (`sub`,
#'   `ses`, `task`, `run`, `echo`, and so on); core manifest columns and aliases
#'   `codec`, `reference`, `grid`, and `support` are also accepted.
#' @param match Matching mode. `"exact"` is the safe default. Regex and glob
#'   matching must be requested explicitly.
#' @param view Manifest namespace returned after filtering.
#' @return A deterministic `lna_query_result` containing at least `record_id`
#'   and archive-relative `lna_path`.
#' @export
query_lna <- function(
    dataset,
    ...,
    match = c("exact", "regex", "glob"),
    view = c(
      "records", "entities", "spatial", "codec", "metadata", "integrity", "all"
    )) {
  match <- match.arg(match)
  view <- match.arg(view)
  dataset <- .query_dataset(dataset)
  filters <- list(...)
  if (length(filters) && (is.null(names(filters)) || any(!nzchar(names(filters))) ||
      anyDuplicated(names(filters)))) {
    .query_abort("filters must have unique non-empty names", "filters")
  }
  data <- dataset$manifest$data
  keep <- rep(TRUE, nrow(data))
  for (name in names(filters)) {
    column <- .query_filter_column(name, names(data))
    if (!column %in% names(data)) {
      keep[] <- FALSE
      next
    }
    patterns <- .query_patterns(filters[[name]], paste0("filters.", name))
    keep <- keep & .query_matches(
      as.character(data[[column]]), patterns, match, paste0("filters.", name)
    )
  }
  columns <- .query_view_columns(data, view)
  out <- data[keep, columns, drop = FALSE]
  out <- out[order(out$lna_path, out$record_id), , drop = FALSE]
  rownames(out) <- NULL
  structure(
    out,
    class = c("lna_query_result", class(out)),
    match = match,
    namespace = view,
    filters = filters
  )
}

#' @export
print.lna_dataset_handle <- function(x, ...) {
  cat("<lna_dataset_handle>", x$declaration$profile_version, "\n")
  cat("  root:", x$root, "\n")
  cat("  records:", nrow(x$manifest$data), "\n")
  invisible(x)
}
