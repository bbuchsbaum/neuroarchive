# Dataset-level grouping, fitting, and provenance for shared spatial bases.

.shared_dataset_abort <- function(message, field, ...) {
  abort_lna(
    message,
    ...,
    field = field,
    .subclass = "lna_error_dataset_shared_basis",
    location = paste0("lna_shared_basis:", field)
  )
}

.shared_effective_support <- function(record, mask = NULL) {
  geometry <- record$geometry
  if (is.null(geometry$support_id)) {
    if (!is.null(mask) && !all(as.logical(mask))) {
      .shared_dataset_abort(
        "a masked support must be asserted by support_id in the source record",
        "support_id", record_id = record$record_id
      )
    }
    full <- array(TRUE, dim = geometry$dims)
    return(list(
      support_id = .lna_support_id(geometry$grid_id, full),
      support_type = "full-grid", mask = full
    ))
  }
  if (is.null(mask)) {
    return(list(
      support_id = geometry$support_id,
      support_type = "masked", mask = NULL
    ))
  }
  mask <- array(as.logical(mask), dim = geometry$dims)
  actual <- .lna_support_id(geometry$grid_id, mask)
  if (!identical(actual, geometry$support_id)) {
    .shared_dataset_abort(
      "support mask does not match the source-record support_id",
      "support_id", record_id = record$record_id,
      expected = geometry$support_id, actual = actual
    )
  }
  list(support_id = actual, support_type = "masked", mask = mask)
}

.shared_reference_key <- function(assertion) {
  paste(assertion$kind, assertion$identifier, assertion$origin, sep = "|")
}

#' Group source records eligible to share a spatial basis
#'
#' Group identity is the conjunction of exact `grid_id`, exact `support_id`,
#' and provenance-asserted reference identity. Record-unique unknown-native
#' references therefore never pool across source images. A missing support is
#' interpreted only as an explicit full-grid support; masked support never
#' falls back to an unknown or implicit mask.
#'
#' @param records Source records.
#' @param support_masks Optional named list of logical masks keyed by record ID.
#'   A mask is required for every record that asserts a non-null support ID.
#' @return A list of `lna_shared_basis_group` objects.
#' @export
lna_shared_basis_groups <- function(records, support_masks = list()) {
  normalized <- validate_lna_source_records(records)
  if (inherits(normalized, "lna_source_record")) normalized <- list(normalized)
  if (!is.list(support_masks) || is.data.frame(support_masks) ||
      (length(support_masks) && (is.null(names(support_masks)) ||
        any(!nzchar(names(support_masks))) || anyDuplicated(names(support_masks))))) {
    .shared_dataset_abort("support_masks must be uniquely named by record ID", "support_masks")
  }
  ids <- vapply(normalized, `[[`, character(1), "record_id")
  unknown <- setdiff(names(support_masks), ids)
  if (length(unknown)) {
    .shared_dataset_abort(
      paste0("support_masks contains unknown record IDs: ", paste(unknown, collapse = ", ")),
      "support_masks"
    )
  }
  supports <- stats::setNames(lapply(normalized, function(record) {
    .shared_effective_support(record, support_masks[[record$record_id]])
  }), ids)
  keys <- vapply(normalized, function(record) {
    support <- supports[[record$record_id]]
    paste(
      record$geometry$grid_id,
      support$support_id,
      .shared_reference_key(record$reference_assertion),
      sep = "|"
    )
  }, character(1))
  partitions <- split(seq_along(normalized), keys)
  groups <- lapply(partitions, function(index) {
    members <- normalized[index]
    member_ids <- vapply(members, `[[`, character(1), "record_id")
    first <- members[[1L]]
    support <- supports[[member_ids[[1L]]]]
    member_masks <- lapply(member_ids, function(id) supports[[id]]$mask)
    if (identical(support$support_type, "masked")) {
      missing <- member_ids[vapply(member_masks, is.null, logical(1))]
      if (length(missing)) {
        .shared_dataset_abort(
          "masked shared-basis groups require an explicit mask for every member",
          "support_masks", record_ids = missing
        )
      }
      digests <- vapply(member_masks, function(mask) {
        .lna_support_id(first$geometry$grid_id, mask)
      }, character(1))
      if (length(unique(digests)) != 1L) {
        .shared_dataset_abort("group members do not have identical support", "support_id")
      }
    }
    group_payload <- list(
      grid_id = first$geometry$grid_id,
      support_id = support$support_id,
      support_type = support$support_type,
      reference_identifier = first$reference_assertion$identifier,
      reference_kind = first$reference_assertion$kind,
      reference_origin = first$reference_assertion$origin,
      record_ids = sort(member_ids)
    )
    structure(list(
      group_id = paste0(
        "shared-", substr(sub("^sha256:", "", .dataset_plan_hash(group_payload)), 1L, 20L)
      ),
      records = members,
      record_ids = sort(member_ids),
      grid_id = group_payload$grid_id,
      support_id = group_payload$support_id,
      support_type = group_payload$support_type,
      mask = support$mask %||% member_masks[[1L]],
      dims = first$geometry$dims,
      affine = first$geometry$affine,
      spatial_units = first$geometry$units,
      voxel_ordering = first$geometry$voxel_ordering,
      reference_assertion = first$reference_assertion,
      group_hash = .dataset_plan_hash(group_payload)
    ), class = c("lna_shared_basis_group", "list"))
  })
  groups[order(vapply(groups, `[[`, character(1), "group_id"))]
}

#' @export
print.lna_shared_basis_group <- function(x, ...) {
  cat("<lna_shared_basis_group>", x$group_id, "\n")
  cat("  records:", length(x$record_ids), "grid:", x$grid_id, "\n")
  cat("  support:", x$support_type, x$support_id, "\n")
  cat("  reference:", x$reference_assertion$identifier, "\n")
  invisible(x)
}

.shared_source_matrix <- function(record, source_root, mask) {
  path <- .writer_resolve_source(list(source = record$source, record_id = record$record_id), source_root)
  if (is.null(path) || !file.exists(path) ||
      !identical(.journal_file_checksum(path), record$source$checksum)) {
    .shared_dataset_abort(
      "shared-basis input is missing or has changed", "source",
      record_id = record$record_id
    )
  }
  image <- neuroim2::read_vec(path)
  values <- as.array(image)
  spatial_dims <- as.integer(dim(values)[seq_len(3L)])
  if (!identical(spatial_dims, as.integer(record$geometry$dims))) {
    .shared_dataset_abort("shared-basis input dimensions changed", "geometry.dims")
  }
  n_volumes <- if (length(dim(values)) >= 4L) dim(values)[[4L]] else 1L
  matrix_values <- matrix(values, nrow = prod(spatial_dims), ncol = n_volumes)
  active <- as.logical(mask)
  out <- t(matrix_values[active, , drop = FALSE])
  if (any(!is.finite(out))) {
    .shared_dataset_abort("shared-basis inputs must contain finite values", "source")
  }
  out
}

.shared_basis_id <- function(group, n_atoms) {
  reference <- gsub("[^A-Za-z0-9+.-]", "-", group$reference_assertion$identifier)
  paste0(
    group$group_id, "-", substr(reference, 1L, 48L), "-k", as.integer(n_atoms)
  )
}

#' Fit one content-addressed basis per compatible source-record group
#'
#' This is the explicit fit pass. Each group is loaded once, pooled along the
#' observation axis, and fitted with a deterministic PCA/SVD basis. The return
#' value is directly consumable as `fit_resources` by [plan_lna_dataset()].
#'
#' @param groups Output from [lna_shared_basis_groups()].
#' @param source_root Source dataset root.
#' @param basis_dir Target `bases` directory, normally `<archive>/bases`.
#' @param n_atoms Maximum number of basis atoms.
#' @param center Whether to center pooled voxel series during fitting and
#'   per-artifact projection.
#' @return A list of shared-basis fit resource receipts.
#' @export
fit_lna_shared_bases <- function(
    groups,
    source_root,
    basis_dir,
    n_atoms = 32L,
    center = TRUE) {
  if (!is.list(groups) || !length(groups) ||
      any(!vapply(groups, inherits, logical(1), "lna_shared_basis_group"))) {
    .shared_dataset_abort("groups must come from lna_shared_basis_groups()", "groups")
  }
  n_atoms <- suppressWarnings(as.integer(n_atoms))
  if (length(n_atoms) != 1L || is.na(n_atoms) || n_atoms < 1L) {
    .shared_dataset_abort("n_atoms must be one positive integer", "n_atoms")
  }
  if (!is.logical(center) || length(center) != 1L || is.na(center)) {
    .shared_dataset_abort("center must be TRUE or FALSE", "center")
  }
  source_root <- .update_normalize_root(source_root, "source_root")
  basis_dir <- .execution_scalar_string(basis_dir, "basis_dir")
  if (!dir.exists(basis_dir) && !dir.create(basis_dir, recursive = TRUE, showWarnings = FALSE)) {
    .shared_dataset_abort("could not create basis_dir", "basis_dir")
  }
  lapply(groups, function(group) {
    mask <- group$mask
    if (is.null(mask)) {
      .shared_dataset_abort("shared-basis fit has no explicit support mask", "mask")
    }
    matrices <- lapply(group$records, .shared_source_matrix, source_root = source_root, mask = mask)
    pooled <- do.call(rbind, matrices)
    fitted <- if (isTRUE(center)) scale(pooled, center = TRUE, scale = FALSE) else pooled
    rank <- min(n_atoms, nrow(fitted), ncol(fitted))
    if (rank < 1L) .shared_dataset_abort("shared-basis group has zero fitted rank", "n_atoms")
    loadings <- svd(fitted, nu = 0L, nv = rank)$v
    source_checksums <- stats::setNames(
      vapply(group$records, function(record) record$source$checksum, character(1)),
      vapply(group$records, `[[`, character(1), "record_id")
    )
    fit_parameters <- list(
      algorithm = "deterministic-base-svd",
      requested_atoms = n_atoms,
      fitted_atoms = rank,
      center = isTRUE(center)
    )
    provenance <- list(
      group_hash = group$group_hash,
      fit_record_ids = sort(group$record_ids),
      source_checksums = as.list(source_checksums),
      reference_identifier = group$reference_assertion$identifier,
      reference_kind = group$reference_assertion$kind,
      reference_origin = group$reference_assertion$origin,
      fit_parameters = fit_parameters
    )
    basis_id <- .shared_basis_id(group, rank)
    artifact <- write_shared_basis(
      list(
        loadings = loadings,
        mask = mask,
        kind = "dataset-pca",
        params = provenance,
        grid = list(
          dims = group$dims, transform = group$affine,
          spatial_units = group$spatial_units,
          voxel_ordering = group$voxel_ordering
        ),
        support_type = group$support_type
      ),
      basis_dir, basis_id = basis_id, basis_kind = "dataset-pca",
      params = provenance, update_registry = TRUE
    )
    list(
      resource_id = group$group_id,
      resource_type = "shared-basis",
      input_record_ids = sort(group$record_ids),
      required_by_record_ids = sort(group$record_ids),
      parameters = list(
        basis_id = basis_id,
        basis_checksum = artifact$checksum,
        basis_filename = artifact$filename,
        basis_kind = "dataset-pca",
        basis_locator = NULL,
        grid_id = artifact$grid_id,
        support_id = artifact$support_id,
        support_type = artifact$support_type,
        reference_identifier = group$reference_assertion$identifier,
        reference_kind = group$reference_assertion$kind,
        reference_origin = group$reference_assertion$origin,
        fit_record_ids = sort(group$record_ids),
        source_checksums = as.list(source_checksums),
        fit_parameters = fit_parameters,
        update_policy = list(
          append = "reuse-if-membership-and-source-checksums-unchanged",
          refresh = "refit-if-membership-or-source-checksum-changed",
          repack = "reuse-only-by-exact-basis-checksum"
        )
      )
    )
  })
}

.dataset_shared_basis_parameters <- function(resource) {
  p <- resource$parameters
  required <- c(
    "basis_id", "basis_checksum", "basis_filename", "basis_kind", "basis_locator",
    "grid_id", "support_id", "support_type", "reference_identifier",
    "reference_kind", "reference_origin", "fit_record_ids", "source_checksums",
    "fit_parameters", "update_policy"
  )
  if (!is.list(p) || !setequal(names(p), required)) {
    .shared_dataset_abort(
      "shared-basis resource parameters are incomplete or contain unknown fields",
      "fit_resources.parameters", resource_id = resource$resource_id
    )
  }
  p$basis_checksum <- .source_record_checksum(p$basis_checksum, "basis_checksum")
  p$grid_id <- .source_record_checksum(p$grid_id, "grid_id")
  p$support_id <- .source_record_checksum(p$support_id, "support_id")
  if (!p$support_type %in% c("masked", "full-grid")) {
    .shared_dataset_abort("shared-basis support_type is invalid", "support_type")
  }
  .validate_basis_id(p$basis_id)
  .validate_registry_filename(p$basis_filename)
  p$fit_record_ids <- sort(.execution_string_vector(
    p$fit_record_ids, "fit_record_ids", allow_empty = FALSE
  ))
  if (!is.list(p$source_checksums) || is.null(names(p$source_checksums)) ||
      any(!nzchar(names(p$source_checksums))) || anyDuplicated(names(p$source_checksums))) {
    .shared_dataset_abort("source_checksums must be named by record ID", "source_checksums")
  }
  p$source_checksums <- stats::setNames(lapply(names(p$source_checksums), function(id) {
    .source_record_checksum(p$source_checksums[[id]], paste0("source_checksums.", id))
  }), names(p$source_checksums))
  if (!is.list(p$fit_parameters) || !is.list(p$update_policy)) {
    .shared_dataset_abort("fit_parameters and update_policy must be objects", "fit_resources")
  }
  p
}

.dataset_apply_shared_basis_resources <- function(entries, fit_resources, codec_profile) {
  uses_shared <- "embed.shared_basis" %in% codec_profile$transforms
  if (!length(fit_resources)) {
    if (uses_shared) {
      .shared_dataset_abort(
        "embed.shared_basis requires a fit resource for every artifact",
        "fit_resources"
      )
    }
    return(entries)
  }
  record_ids <- vapply(entries, `[[`, character(1), "record_id")
  resources <- .execution_fit_resources(fit_resources, record_ids)
  shared <- resources[vapply(resources, function(resource) {
    identical(resource$resource_type, "shared-basis")
  }, logical(1))]
  if (!length(shared)) {
    if (uses_shared) {
      .shared_dataset_abort(
        "embed.shared_basis requires shared-basis fit resources",
        "fit_resources"
      )
    }
    return(entries)
  }
  if (!uses_shared) {
    .shared_dataset_abort(
      "shared-basis fit resources require the embed.shared_basis transform",
      "codec_profile.transforms"
    )
  }
  assigned <- character()
  for (resource in shared) {
    p <- .dataset_shared_basis_parameters(resource)
    if (!setequal(resource$input_record_ids, p$fit_record_ids)) {
      .shared_dataset_abort("fit membership disagrees with basis provenance", "fit_record_ids")
    }
    overlap <- intersect(assigned, resource$required_by_record_ids)
    if (length(overlap)) {
      .shared_dataset_abort("a record cannot use multiple shared bases", "fit_resources")
    }
    assigned <- c(assigned, resource$required_by_record_ids)
    for (record_id in resource$required_by_record_ids) {
      index <- match(record_id, record_ids)
      entry <- entries[[index]]
      effective_support <- entry$geometry$support_id
      if (is.null(effective_support)) {
        full <- array(TRUE, dim = as.integer(entry$geometry$dims))
        effective_support <- .lna_support_id(entry$geometry$grid_id, full)
        if (!identical(p$support_type, "full-grid")) {
          .shared_dataset_abort(
            "masked shared basis requires an asserted source-record support_id",
            "support_id", record_id = record_id
          )
        }
        entry$geometry$support_id <- effective_support
      }
      if (!identical(entry$geometry$grid_id, p$grid_id) ||
          !identical(effective_support, p$support_id)) {
        .shared_dataset_abort(
          "shared-basis grid or support does not match its consuming record",
          "spatial_identity", record_id = record_id,
          expected_grid_id = p$grid_id, actual_grid_id = entry$geometry$grid_id,
          expected_support_id = p$support_id, actual_support_id = effective_support
        )
      }
      if (!identical(entry$reference_assertion$identifier, p$reference_identifier) ||
          !identical(entry$reference_assertion$kind, p$reference_kind)) {
        .shared_dataset_abort(
          "shared-basis reference assertion does not match its consuming record",
          "reference_identifier", record_id = record_id
        )
      }
      expected_source <- p$source_checksums[[record_id]]
      if (is.null(expected_source) || !identical(expected_source, entry$source$checksum)) {
        .shared_dataset_abort(
          "shared-basis fit source checksum disagrees with its consuming record",
          "source_checksums", record_id = record_id
        )
      }
      entry$shared_basis <- c(
        list(status = "shared", resource_id = resource$resource_id), p
      )
      params <- entry$codec_transform_params %||% list()
      existing <- params[["embed.shared_basis"]] %||% list()
      params[["embed.shared_basis"]] <- utils::modifyList(existing, list(
        basis_id = p$basis_id,
        basis_checksum = p$basis_checksum,
        basis_kind = p$basis_kind,
        basis_locator = p$basis_locator,
        center = isTRUE(p$fit_parameters$center)
      ))
      entry$codec_transform_params <- params
      entries[[index]] <- entry
    }
  }
  missing <- setdiff(record_ids, assigned)
  if (length(missing)) {
    .shared_dataset_abort(
      "every artifact in a shared-basis codec profile requires exactly one basis resource",
      "fit_resources", record_ids = missing
    )
  }
  entries
}

.dataset_verify_shared_basis_resource <- function(resource, plan, archive_root) {
  if (!identical(resource$resource_type, "shared-basis")) return(invisible(TRUE))
  p <- .dataset_shared_basis_parameters(resource)
  file <- file.path(archive_root, "bases", p$basis_filename)
  if (!file.exists(file)) {
    .shared_dataset_abort("planned shared-basis artifact is missing", "basis_file")
  }
  basis <- read_shared_basis(file, verify_checksum = TRUE)
  mismatches <- c(
    if (!identical(basis$checksum, p$basis_checksum)) "checksum" else character(),
    if (!identical(basis$grid_id, p$grid_id)) "grid_id" else character(),
    if (!identical(basis$support_id, p$support_id)) "support_id" else character(),
    if (!identical(basis$support_type, p$support_type)) "support_type" else character()
  )
  if (length(mismatches)) {
    .shared_dataset_abort(
      paste0("shared-basis artifact disagrees with fit receipt: ", paste(mismatches, collapse = ", ")),
      "basis_file", resource_id = resource$resource_id
    )
  }
  invisible(TRUE)
}
