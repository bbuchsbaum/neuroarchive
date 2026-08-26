# Explicit preflight-fit-write-validate-finalize execution graph.

.lna_execution_schema_version <- "1.0.0"
.lna_execution_phases <- c("preflight", "fit", "write", "validate", "finalize")

.execution_abort <- function(message, field, ...) {
  abort_lna(
    message,
    ...,
    field = field,
    .subclass = "lna_error_execution",
    location = paste0("lna_execution:", field)
  )
}

.execution_scalar_string <- function(x, field) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    .execution_abort(paste0(field, " must be one non-empty string"), field)
  }
  enc2utf8(x)
}

.execution_string_vector <- function(x, field, allow_empty = TRUE) {
  if (allow_empty && (is.null(x) || (is.list(x) && !length(x)))) {
    return(character())
  }
  if (is.list(x)) x <- unlist(x, use.names = FALSE)
  if (!is.character(x) || anyNA(x) || any(!nzchar(x)) || anyDuplicated(x) ||
      (!allow_empty && !length(x))) {
    .execution_abort(
      paste0(field, " must contain unique non-empty strings"), field
    )
  }
  enc2utf8(x)
}

.execution_id <- function(x, field) {
  x <- .execution_scalar_string(x, field)
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9._:-]*$", x)) {
    .execution_abort(
      paste0(field, " contains an invalid execution identifier"), field
    )
  }
  x
}

.execution_fit_resource <- function(resource, index, record_ids) {
  field <- paste0("fit_resources[[", index, "]]")
  if (!is.list(resource) || is.data.frame(resource) || is.null(names(resource))) {
    .execution_abort(paste0(field, " must be a named object"), field)
  }
  allowed <- c(
    "resource_id", "resource_type", "input_record_ids",
    "required_by_record_ids", "parameters"
  )
  unknown <- setdiff(names(resource), allowed)
  required <- c("resource_id", "resource_type", "input_record_ids")
  missing <- setdiff(required, names(resource))
  if (length(unknown) || length(missing)) {
    .execution_abort(
      paste0(
        field, " must use fields ", paste(allowed, collapse = ", "),
        if (length(missing)) paste0("; missing ", paste(missing, collapse = ", ")) else "",
        if (length(unknown)) paste0("; unknown ", paste(unknown, collapse = ", ")) else ""
      ),
      field
    )
  }
  resource_id <- .execution_id(resource$resource_id, paste0(field, ".resource_id"))
  resource_type <- .execution_scalar_string(
    resource$resource_type, paste0(field, ".resource_type")
  )
  input_record_ids <- sort(.execution_string_vector(
    resource$input_record_ids, paste0(field, ".input_record_ids"),
    allow_empty = FALSE
  ))
  required_by <- resource$required_by_record_ids %||% record_ids
  required_by <- sort(.execution_string_vector(
    required_by, paste0(field, ".required_by_record_ids"), allow_empty = FALSE
  ))
  absent <- setdiff(c(input_record_ids, required_by), record_ids)
  if (length(absent)) {
    .execution_abort(
      paste0(field, " refers to unknown record(s): ", paste(unique(absent), collapse = ", ")),
      field
    )
  }
  parameters <- resource$parameters %||% list()
  tryCatch(
    .source_record_json_value(parameters, paste0(field, ".parameters")),
    lna_error_source_record = function(e) {
      .execution_abort(conditionMessage(e), paste0(field, ".parameters"))
    }
  )
  list(
    resource_id = resource_id,
    resource_type = resource_type,
    input_record_ids = input_record_ids,
    required_by_record_ids = required_by,
    parameters = parameters
  )
}

.execution_fit_resources <- function(resources, record_ids) {
  if (!is.list(resources) || is.data.frame(resources)) {
    .execution_abort("fit_resources must be a list", "fit_resources")
  }
  if (!length(resources)) return(list())
  normalized <- lapply(seq_along(resources), function(i) {
    .execution_fit_resource(resources[[i]], i, record_ids)
  })
  ids <- vapply(normalized, `[[`, character(1), "resource_id")
  if (anyDuplicated(ids)) {
    .execution_abort("fit_resources have duplicate resource_id values", "fit_resources")
  }
  normalized[order(ids)]
}

.execution_node <- function(
    node_id,
    phase,
    kind,
    depends_on = character(),
    record_id = NULL,
    resource_id = NULL,
    source_checksum = NULL,
    artifact_path = NULL,
    required_fit_resources = character()) {
  list(
    node_id = node_id,
    phase = phase,
    kind = kind,
    depends_on = sort(depends_on),
    record_id = record_id,
    resource_id = resource_id,
    source_checksum = source_checksum,
    artifact_path = artifact_path,
    required_fit_resources = sort(required_fit_resources)
  )
}

.execution_topological_order <- function(nodes) {
  ids <- vapply(nodes, `[[`, character(1), "node_id")
  dependencies <- stats::setNames(lapply(nodes, `[[`, "depends_on"), ids)
  phase_by_id <- stats::setNames(vapply(nodes, `[[`, character(1), "phase"), ids)
  completed <- character()
  remaining <- ids
  while (length(remaining)) {
    ready <- remaining[vapply(remaining, function(id) {
      all(dependencies[[id]] %in% completed)
    }, logical(1))]
    if (!length(ready)) {
      .execution_abort(
        paste0("execution graph contains a dependency cycle involving: ",
          paste(sort(remaining), collapse = ", ")),
        "nodes"
      )
    }
    rank <- match(phase_by_id[ready], .lna_execution_phases)
    ready <- ready[order(rank, ready)]
    completed <- c(completed, ready)
    remaining <- setdiff(remaining, ready)
  }
  completed
}

.build_lna_execution <- function(entries, fit_resources = list()) {
  record_ids <- vapply(entries, `[[`, character(1), "record_id")
  resources <- .execution_fit_resources(fit_resources, record_ids)
  nodes <- list(.execution_node(
    "phase:preflight", "preflight", "inventory-and-preflight"
  ))
  if (length(resources)) {
    nodes <- c(nodes, lapply(resources, function(resource) {
      .execution_node(
        paste0("fit:", resource$resource_id),
        "fit",
        "fit-resource",
        depends_on = "phase:preflight",
        resource_id = resource$resource_id
      )
    }))
  }
  fit_nodes <- if (length(resources)) {
    paste0("fit:", vapply(resources, `[[`, character(1), "resource_id"))
  } else {
    character()
  }
  nodes[[length(nodes) + 1L]] <- .execution_node(
    "phase:fit-complete",
    "fit",
    "phase-barrier",
    depends_on = c("phase:preflight", fit_nodes)
  )
  write_nodes <- lapply(entries, function(entry) {
    required <- vapply(resources[vapply(resources, function(resource) {
      entry$record_id %in% resource$required_by_record_ids
    }, logical(1))], `[[`, character(1), "resource_id")
    .execution_node(
      paste0("write:", entry$record_id),
      "write",
      "artifact-write",
      depends_on = "phase:fit-complete",
      record_id = entry$record_id,
      source_checksum = entry$source$checksum,
      artifact_path = entry$destination_relative,
      required_fit_resources = required
    )
  })
  nodes <- c(nodes, write_nodes)
  write_ids <- vapply(write_nodes, `[[`, character(1), "node_id")
  nodes[[length(nodes) + 1L]] <- .execution_node(
    "phase:write-complete",
    "write",
    "phase-barrier",
    depends_on = write_ids
  )
  validation_nodes <- lapply(entries, function(entry) {
    .execution_node(
      paste0("validate:", entry$record_id),
      "validate",
      "artifact-validation",
      depends_on = c("phase:write-complete", paste0("write:", entry$record_id)),
      record_id = entry$record_id,
      source_checksum = entry$source$checksum,
      artifact_path = entry$destination_relative
    )
  })
  nodes <- c(nodes, validation_nodes)
  nodes[[length(nodes) + 1L]] <- .execution_node(
    "phase:finalize",
    "finalize",
    "manifest-finalization",
    depends_on = vapply(validation_nodes, `[[`, character(1), "node_id")
  )
  ids <- vapply(nodes, `[[`, character(1), "node_id")
  nodes <- nodes[order(ids)]
  execution <- list(
    schema_version = .lna_execution_schema_version,
    phase_order = .lna_execution_phases,
    fit_resources = resources,
    nodes = nodes
  )
  execution$topological_order <- .execution_topological_order(nodes)
  validate_lna_execution(execution, record_ids)
}

.normalize_execution_node <- function(node, index) {
  field <- paste0("nodes[[", index, "]]")
  if (!is.list(node) || is.data.frame(node)) {
    .execution_abort(paste0(field, " must be an object"), field)
  }
  required <- c(
    "node_id", "phase", "kind", "depends_on", "record_id", "resource_id",
    "source_checksum", "artifact_path", "required_fit_resources"
  )
  missing <- setdiff(required, names(node))
  unknown <- setdiff(names(node), required)
  if (length(missing) || length(unknown)) {
    .execution_abort(
      paste0(field, " must contain exactly: ", paste(required, collapse = ", ")),
      field
    )
  }
  node_id <- .execution_id(node$node_id, paste0(field, ".node_id"))
  phase <- .execution_scalar_string(node$phase, paste0(field, ".phase"))
  if (!phase %in% .lna_execution_phases) {
    .execution_abort(paste0(field, ".phase is unsupported"), paste0(field, ".phase"))
  }
  kind <- .execution_scalar_string(node$kind, paste0(field, ".kind"))
  nullable_id <- function(x, name) {
    if (is.null(x)) return(NULL)
    .execution_id(x, paste0(field, ".", name))
  }
  source_checksum <- node$source_checksum
  if (!is.null(source_checksum)) {
    source_checksum <- tryCatch(
      .source_record_checksum(source_checksum, paste0(field, ".source_checksum")),
      lna_error_source_record = function(e) {
        .execution_abort(conditionMessage(e), paste0(field, ".source_checksum"))
      }
    )
  }
  artifact_path <- node$artifact_path
  if (!is.null(artifact_path)) {
    artifact_path <- .manifest_validate_path(
      artifact_path, paste0(field, ".artifact_path")
    )
  }
  .execution_node(
    node_id = node_id,
    phase = phase,
    kind = kind,
    depends_on = .execution_string_vector(
      node$depends_on, paste0(field, ".depends_on")
    ),
    record_id = nullable_id(node$record_id, "record_id"),
    resource_id = nullable_id(node$resource_id, "resource_id"),
    source_checksum = source_checksum,
    artifact_path = artifact_path,
    required_fit_resources = .execution_string_vector(
      node$required_fit_resources, paste0(field, ".required_fit_resources")
    )
  )
}

.execution_expect_node <- function(
    node,
    phase,
    kind,
    depends_on,
    field,
    record_id = NULL,
    resource_id = NULL,
    required_fit_resources = character()) {
  mismatches <- character()
  if (!identical(node$phase, phase)) mismatches <- c(mismatches, "phase")
  if (!identical(node$kind, kind)) mismatches <- c(mismatches, "kind")
  if (!setequal(node$depends_on, depends_on)) mismatches <- c(mismatches, "depends_on")
  if (!identical(node$record_id, record_id)) mismatches <- c(mismatches, "record_id")
  if (!identical(node$resource_id, resource_id)) mismatches <- c(mismatches, "resource_id")
  if (!setequal(node$required_fit_resources, required_fit_resources)) {
    mismatches <- c(mismatches, "required_fit_resources")
  }
  if (length(mismatches)) {
    .execution_abort(
      paste0(field, " has incompatible ", paste(mismatches, collapse = ", ")),
      field
    )
  }
  invisible(TRUE)
}

#' Validate an LNA dataset execution graph
#'
#' @param execution Execution graph from a dataset plan.
#' @param record_ids Record identities represented by the plan.
#' @return The normalized execution graph.
#' @export
validate_lna_execution <- function(execution, record_ids) {
  if (!is.list(execution) || is.data.frame(execution)) {
    .execution_abort("execution must be an object", "execution")
  }
  required <- c(
    "schema_version", "phase_order", "fit_resources", "nodes",
    "topological_order"
  )
  missing <- setdiff(required, names(execution))
  if (length(missing)) {
    .execution_abort(
      paste0("execution is missing: ", paste(missing, collapse = ", ")),
      "execution"
    )
  }
  version <- .execution_scalar_string(execution$schema_version, "schema_version")
  major <- suppressWarnings(as.integer(strsplit(version, ".", fixed = TRUE)[[1L]][[1L]]))
  if (!grepl("^[0-9]+[.][0-9]+[.][0-9]+$", version) || !identical(major, 1L)) {
    .execution_abort(
      paste0("unsupported execution schema version: ", version), "schema_version"
    )
  }
  phase_order <- .execution_string_vector(execution$phase_order, "phase_order")
  if (!identical(phase_order, .lna_execution_phases)) {
    .execution_abort("phase_order is incompatible", "phase_order")
  }
  record_ids <- sort(.execution_string_vector(record_ids, "record_ids", allow_empty = FALSE))
  resources <- .execution_fit_resources(execution$fit_resources, record_ids)
  if (!is.list(execution$nodes) || is.data.frame(execution$nodes) || !length(execution$nodes)) {
    .execution_abort("nodes must be a non-empty list", "nodes")
  }
  nodes <- lapply(seq_along(execution$nodes), function(i) {
    .normalize_execution_node(execution$nodes[[i]], i)
  })
  ids <- vapply(nodes, `[[`, character(1), "node_id")
  if (anyDuplicated(ids)) .execution_abort("node_id values must be unique", "nodes")
  missing_dependencies <- setdiff(
    unique(unlist(lapply(nodes, `[[`, "depends_on"), use.names = FALSE)), ids
  )
  if (length(missing_dependencies)) {
    .execution_abort(
      paste0("missing dependency node(s): ", paste(missing_dependencies, collapse = ", ")),
      "nodes"
    )
  }
  for (node in nodes) {
    if (node$node_id %in% node$depends_on) {
      .execution_abort(
        paste0("node depends on itself: ", node$node_id), "nodes"
      )
    }
    if (!is.null(node$record_id) && !node$record_id %in% record_ids) {
      .execution_abort(
        paste0("node refers to unknown record: ", node$record_id), "nodes"
      )
    }
  }
  computed <- .execution_topological_order(nodes)
  resource_ids <- vapply(resources, `[[`, character(1), "resource_id")
  expected_fit <- if (length(resource_ids)) paste0("fit:", resource_ids) else character()
  expected_write <- paste0("write:", record_ids)
  expected_validate <- paste0("validate:", record_ids)
  special <- c(
    "phase:preflight", "phase:fit-complete", "phase:write-complete",
    "phase:finalize"
  )
  expected_ids <- c(special, expected_fit, expected_write, expected_validate)
  if (!setequal(ids, expected_ids)) {
    .execution_abort(
      paste0(
        "execution node set is incompatible; missing: ",
        paste(setdiff(expected_ids, ids), collapse = ", "),
        "; unexpected: ", paste(setdiff(ids, expected_ids), collapse = ", ")
      ),
      "nodes"
    )
  }
  node_map <- stats::setNames(nodes, ids)
  .execution_expect_node(
    node_map[["phase:preflight"]], "preflight", "inventory-and-preflight",
    character(), "nodes.phase:preflight"
  )
  for (resource in resources) {
    id <- paste0("fit:", resource$resource_id)
    .execution_expect_node(
      node_map[[id]], "fit", "fit-resource", "phase:preflight",
      paste0("nodes.", id), resource_id = resource$resource_id
    )
  }
  .execution_expect_node(
    node_map[["phase:fit-complete"]], "fit", "phase-barrier",
    c("phase:preflight", expected_fit), "nodes.phase:fit-complete"
  )
  for (record_id in record_ids) {
    id <- paste0("write:", record_id)
    required_resources <- resource_ids[vapply(resources, function(resource) {
      record_id %in% resource$required_by_record_ids
    }, logical(1))]
    .execution_expect_node(
      node_map[[id]], "write", "artifact-write", "phase:fit-complete",
      paste0("nodes.", id), record_id = record_id,
      required_fit_resources = required_resources
    )
  }
  .execution_expect_node(
    node_map[["phase:write-complete"]], "write", "phase-barrier",
    expected_write, "nodes.phase:write-complete"
  )
  for (record_id in record_ids) {
    id <- paste0("validate:", record_id)
    .execution_expect_node(
      node_map[[id]], "validate", "artifact-validation",
      c("phase:write-complete", paste0("write:", record_id)),
      paste0("nodes.", id), record_id = record_id
    )
  }
  .execution_expect_node(
    node_map[["phase:finalize"]], "finalize", "manifest-finalization",
    expected_validate, "nodes.phase:finalize"
  )
  supplied <- .execution_string_vector(
    execution$topological_order, "topological_order", allow_empty = FALSE
  )
  if (!identical(supplied, computed)) {
    .execution_abort(
      "topological_order is not the deterministic order for this graph",
      "topological_order"
    )
  }
  nodes <- nodes[order(ids)]
  list(
    schema_version = version,
    phase_order = phase_order,
    fit_resources = resources,
    nodes = nodes,
    topological_order = computed
  )
}

.dataset_plan_payload <- function(plan) {
  payload <- unclass(plan)
  payload$plan_hash <- NULL
  payload
}

#' Validate a serializable LNA dataset plan
#'
#' @param plan A dataset plan or decoded plan object.
#' @return The normalized `lna_dataset_plan`.
#' @export
validate_lna_dataset_plan <- function(plan) {
  if (!is.list(plan) || is.data.frame(plan)) {
    .execution_abort("plan must be an object", "plan")
  }
  required <- c(
    "schema_version", "archive_mode", "output_root", "codec_profile",
    "metadata_policy", "artifacts", "execution", "summary", "plan_hash"
  )
  missing <- setdiff(required, names(plan))
  if (length(missing)) {
    .execution_abort(
      paste0("plan is missing: ", paste(missing, collapse = ", ")), "plan"
    )
  }
  major <- suppressWarnings(as.integer(
    strsplit(plan$schema_version, ".", fixed = TRUE)[[1L]][[1L]]
  ))
  if (!is.character(plan$schema_version) || length(plan$schema_version) != 1L ||
      !grepl("^[0-9]+[.][0-9]+[.][0-9]+$", plan$schema_version) ||
      !identical(major, 1L)) {
    .execution_abort(
      paste0("unsupported dataset plan schema version: ", plan$schema_version),
      "schema_version"
    )
  }
  if (!identical(plan$archive_mode, "one-source-one-lna")) {
    .execution_abort("archive_mode is unsupported", "archive_mode")
  }
  plan$codec_profile <- .dataset_plan_codec(plan$codec_profile)
  if (!is.list(plan$artifacts) || !length(plan$artifacts)) {
    .execution_abort("plan.artifacts must be a non-empty list", "artifacts")
  }
  record_ids <- vapply(plan$artifacts, function(entry) {
    .execution_id(entry$record_id, "artifacts.record_id")
  }, character(1))
  if (anyDuplicated(record_ids)) {
    .execution_abort("artifact record_id values must be unique", "artifacts")
  }
  paths <- vapply(plan$artifacts, function(entry) {
    .manifest_validate_path(entry$destination_relative, "artifacts.destination_relative")
  }, character(1))
  if (anyDuplicated(paths)) {
    .execution_abort("artifact destination paths must be unique", "artifacts")
  }
  plan$metadata_policy <- unclass(validate_lna_metadata_policy(
    plan$metadata_policy
  ))
  validate_lna_metadata_closure(plan)
  normalized_metadata <- .metadata_apply_plan(
    plan$artifacts, plan$metadata_policy
  )
  plan$artifacts <- normalized_metadata$entries
  plan$metadata_policy <- normalized_metadata$policy
  plan$execution <- validate_lna_execution(plan$execution, record_ids)
  node_map <- .execution_nodes_by_id(plan)
  for (entry in plan$artifacts) {
    source_checksum <- tryCatch(
      .source_record_checksum(entry$source$checksum, "artifacts.source.checksum"),
      lna_error_source_record = function(e) {
        .execution_abort(conditionMessage(e), "artifacts.source.checksum")
      }
    )
    for (prefix in c("write", "validate")) {
      node <- node_map[[paste0(prefix, ":", entry$record_id)]]
      if (!identical(node$source_checksum, source_checksum) ||
          !identical(node$artifact_path, entry$destination_relative)) {
        .execution_abort(
          paste0(prefix, " node disagrees with artifact identity for ", entry$record_id),
          "execution.nodes"
        )
      }
    }
  }
  if (!is.list(plan$summary) ||
      !identical(as.integer(plan$summary$n_records), length(plan$artifacts)) ||
      !identical(as.integer(plan$summary$n_artifacts), length(plan$artifacts))) {
    .execution_abort("plan.summary disagrees with artifacts", "summary")
  }
  expected_hash <- .dataset_plan_hash(.dataset_plan_payload(plan))
  if (!identical(plan$plan_hash, expected_hash)) {
    .execution_abort(
      "plan_hash does not match the serialized plan payload", "plan_hash",
      expected_hash = expected_hash, actual_hash = plan$plan_hash
    )
  }
  structure(plan, class = c("lna_dataset_plan", "list"))
}

#' Serialize an LNA dataset plan
#'
#' @param plan A validated `lna_dataset_plan`.
#' @param path Destination JSON path.
#' @param overwrite Whether an existing file may be replaced.
#' @return `path`, invisibly.
#' @export
write_lna_dataset_plan <- function(plan, path, overwrite = FALSE) {
  plan <- validate_lna_dataset_plan(plan)
  path <- .execution_scalar_string(path, "path")
  if (file.exists(path) && !isTRUE(overwrite)) {
    .execution_abort("dataset plan file already exists", "path")
  }
  parent <- dirname(path)
  if (!dir.exists(parent) &&
      !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
    .execution_abort("could not create dataset plan directory", "path")
  }
  jsonlite::write_json(
    unclass(plan), path, auto_unbox = TRUE, null = "null", digits = NA,
    pretty = TRUE
  )
  invisible(path)
}

#' @rdname write_lna_dataset_plan
#' @export
read_lna_dataset_plan <- function(path) {
  path <- .execution_scalar_string(path, "path")
  if (!file.exists(path)) .execution_abort("dataset plan file does not exist", "path")
  decoded <- tryCatch(
    jsonlite::read_json(path, simplifyVector = FALSE),
    error = function(e) .execution_abort("dataset plan is invalid JSON", "path")
  )
  validate_lna_dataset_plan(decoded)
}

.execution_nodes_by_id <- function(plan) {
  stats::setNames(
    plan$execution$nodes,
    vapply(plan$execution$nodes, `[[`, character(1), "node_id")
  )
}

.validate_execution_state <- function(state, plan) {
  plan <- validate_lna_dataset_plan(plan)
  if (!is.list(state) || !identical(state$plan_hash, plan$plan_hash)) {
    .execution_abort("execution state belongs to a different plan", "state")
  }
  completed <- .execution_string_vector(
    state$completed_nodes %||% character(), "state.completed_nodes"
  )
  nodes <- .execution_nodes_by_id(plan)
  absent <- setdiff(completed, names(nodes))
  if (length(absent)) {
    .execution_abort(
      paste0("execution state contains unknown node(s): ", paste(absent, collapse = ", ")),
      "state.completed_nodes"
    )
  }
  for (node_id in completed) {
    missing <- setdiff(nodes[[node_id]]$depends_on, completed)
    if (length(missing)) {
      .execution_abort(
        paste0("completed node ", node_id, " lacks completed dependencies: ",
          paste(missing, collapse = ", ")),
        "state.completed_nodes"
      )
    }
  }
  structure(
    list(
      schema_version = .lna_execution_schema_version,
      plan_hash = plan$plan_hash,
      completed_nodes = completed
    ),
    class = c("lna_execution_state", "list")
  )
}

#' Traverse a dataset execution plan safely
#'
#' @param plan A validated `lna_dataset_plan`.
#' @param completed_nodes Dependency-closed node IDs already completed.
#' @param state An `lna_execution_state`.
#' @param node_id Node to mark complete.
#' @param phase Optional phase filter for ready nodes.
#' @return `lna_execution_state()` and `advance_lna_execution()` return state;
#'   `ready_lna_execution_nodes()` returns ready node records.
#' @export
lna_execution_state <- function(plan, completed_nodes = character()) {
  .validate_execution_state(
    list(plan_hash = plan$plan_hash, completed_nodes = completed_nodes), plan
  )
}

#' @rdname lna_execution_state
#' @export
ready_lna_execution_nodes <- function(plan, state, phase = NULL) {
  plan <- validate_lna_dataset_plan(plan)
  state <- .validate_execution_state(state, plan)
  if (!is.null(phase)) {
    phase <- .execution_scalar_string(phase, "phase")
    if (!phase %in% .lna_execution_phases) {
      .execution_abort("phase is unsupported", "phase")
    }
  }
  nodes <- .execution_nodes_by_id(plan)
  remaining <- setdiff(plan$execution$topological_order, state$completed_nodes)
  ready <- remaining[vapply(remaining, function(node_id) {
    all(nodes[[node_id]]$depends_on %in% state$completed_nodes)
  }, logical(1))]
  if (!is.null(phase)) {
    ready <- ready[vapply(ready, function(id) identical(nodes[[id]]$phase, phase), logical(1))]
  }
  unname(nodes[ready])
}

#' @rdname lna_execution_state
#' @export
advance_lna_execution <- function(plan, state, node_id) {
  plan <- validate_lna_dataset_plan(plan)
  state <- .validate_execution_state(state, plan)
  node_id <- .execution_id(node_id, "node_id")
  nodes <- .execution_nodes_by_id(plan)
  if (!node_id %in% names(nodes)) {
    .execution_abort(paste0("unknown execution node: ", node_id), "node_id")
  }
  if (node_id %in% state$completed_nodes) return(state)
  missing <- setdiff(nodes[[node_id]]$depends_on, state$completed_nodes)
  if (length(missing)) {
    .execution_abort(
      paste0("cannot complete ", node_id, "; incomplete dependencies: ",
        paste(missing, collapse = ", ")),
      "node_id", incomplete_dependencies = missing
    )
  }
  state$completed_nodes <- c(state$completed_nodes, node_id)
  .validate_execution_state(state, plan)
}
