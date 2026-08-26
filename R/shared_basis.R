# Shared basis artifact storage, registry, and resolution
#
# Shared bases are stored as standalone HDF5 files (.lna_basis.h5) with a
# content-addressed filename (sha256-<hash>.lna_basis.h5). A JSON registry
# (bases/registry.json) provides human-readable lookup by basis_id aliases.
#
# Resolution order for locating a shared basis:
#   1. basis_locator path relative to the consuming archive
#   2. bases/registry.json in the dataset directory
#   3. user-configured search paths (lna_options("basis_search_paths"))
#   4. fail with actionable error

#' Write a shared basis artifact to disk
#'
#' Serializes a fmrilatent template (or raw loadings matrix + mask) to a
#' standalone `.lna_basis.h5` file with content-addressed filename. Optionally
#' updates a `registry.json` in the target directory.
#'
#' @param template A fmrilatent template object (HierarchicalBasisTemplate or
#'   ParcelBasisTemplate), or a list with components `loadings` (Matrix) and
#'   `mask` (LogicalNeuroVol or logical array). Raw lists may also provide a
#'   `grid` list with `dims`, `transform`, `spatial_units`, and
#'   `voxel_ordering`, plus `support_type = "masked"` or the explicit
#'   `"full-grid"`. A `LogicalNeuroVol` supplies its grid automatically.
#' @param dir Directory to write into (typically `derivatives/lna/bases/`).
#' @param basis_id Human-readable identifier (e.g., "schaefer400_slepian-k8").
#' @param basis_kind Character string describing the basis family. Inferred
#'   from template metadata when possible.
#' @param params Optional list of parameters stored as JSON in the artifact.
#' @param update_registry Logical; if TRUE (default), update or create
#'   `registry.json` in `dir`.
#' @return Invisibly, a list with `file`, `checksum`, `basis_id`, `filename`,
#'   `grid_id`, `support_id`, and `support_type`.
#' @export
write_shared_basis <- function(template, dir, basis_id,
                               basis_kind = NULL, params = NULL,
                               update_registry = TRUE) {
  .validate_basis_id(basis_id)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  # Capture timestamp once for consistent use in HDF5 and registry

  created_str <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")

  # Extract loadings and mask from template or raw list
  extracted <- .extract_basis_components(template)
  loadings <- extracted$loadings
  mask_arr <- extracted$mask_arr
  if (is.null(basis_kind)) basis_kind <- extracted$kind
  if (is.null(params)) params <- extracted$params
  spatial_identity <- .shared_basis_spatial_identity(extracted, mask_arr)

  # Compute artifact checksum from the immutable basis contents.
  checksum <- .compute_artifact_checksum(
    loadings,
    mask_arr,
    basis_kind,
    params,
    spatial_identity = spatial_identity
  )
  checksum_str <- paste0("sha256:", checksum)
  mask_checksum <- .compute_mask_checksum(mask_arr)

  # Content-addressed filename
  filename <- paste0("sha256-", checksum, ".lna_basis.h5")
  filepath <- file.path(dir, filename)

  if (file.exists(filepath)) {
    if (!.verify_basis_checksum(filepath, checksum_str)) {
      abort_lna(
        sprintf("Existing basis artifact does not match expected checksum: %s", filepath),
        .subclass = "lna_error_checksum",
        location = "write_shared_basis"
      )
    }
  } else {
    h5 <- hdf5r::H5File$new(filepath, mode = "w")
    on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)

    storage_kind <- .shared_basis_storage_kind(loadings)
    params_json <- .encode_basis_params(params)

    # /meta group with immutable artifact metadata
    meta_grp <- h5$create_group("meta")
    meta_grp$create_attr("basis_kind", basis_kind %||% "")
    meta_grp$create_attr("checksum", checksum_str)
    meta_grp$create_attr("mask_checksum", mask_checksum)
    meta_grp$create_attr("n_atoms", as.integer(ncol(loadings)))
    meta_grp$create_attr("n_voxels", as.integer(nrow(loadings)))
    meta_grp$create_attr("created", created_str)
    meta_grp$create_attr("loadings_storage", storage_kind)
    if (nzchar(spatial_identity$grid_id)) {
      meta_grp$create_attr("grid_id", spatial_identity$grid_id)
      meta_grp$create_attr("support_id", spatial_identity$support_id)
      meta_grp$create_attr("support_type", spatial_identity$support_type)
      meta_grp$create_attr("spatial_units", spatial_identity$spatial_units)
      meta_grp$create_attr("voxel_ordering", spatial_identity$voxel_ordering)
      meta_grp$create_attr("grid_dims", spatial_identity$dims)
    }
    safe_h5_close(meta_grp)

    root <- h5[["/"]]
    .write_shared_basis_loadings(root, loadings, storage_kind)

    # /mask dataset (store as integer array preserving dimensions)
    mask_int <- array(as.integer(mask_arr), dim = dim(mask_arr))
    h5_write_dataset(root, "mask", mask_int)

    if (!is.null(params_json)) {
      h5_write_dataset(root, "params", matrix(params_json, nrow = 1L, ncol = 1L))
    }
    safe_h5_close(root)
  }

  # Update registry
  if (update_registry) {
    .update_basis_registry(dir, basis_id, list(
      checksum = checksum_str,
      filename = filename,
      kind = basis_kind %||% "",
      mask_checksum = mask_checksum,
      grid_id = spatial_identity$grid_id,
      support_id = spatial_identity$support_id,
      support_type = spatial_identity$support_type,
      n_atoms = as.integer(ncol(loadings)),
      n_voxels = as.integer(nrow(loadings)),
      created = created_str
    ))
  }

  invisible(list(
    file = filepath,
    checksum = checksum_str,
    basis_id = basis_id,
    filename = filename,
    grid_id = spatial_identity$grid_id,
    support_id = spatial_identity$support_id,
    support_type = spatial_identity$support_type
  ))
}

#' Read a shared basis artifact from disk
#'
#' Opens an `.lna_basis.h5` file and returns the loadings matrix, mask, and
#' metadata.
#'
#' @param file Path to the `.lna_basis.h5` file.
#' @param verify_checksum Logical; if TRUE (default), verify the content
#'   checksum matches the stored value.
#' @return A list with `loadings` (matrix or sparse Matrix), `mask` (logical array),
#'   `basis_id`, `basis_ids`, `basis_kind`, `checksum`, `mask_checksum`,
#'   `grid_id`, `support_id`, `support_type`, and `params`.
#' @export
read_shared_basis <- function(file, verify_checksum = TRUE) {
  if (!file.exists(file)) {
    abort_lna(
      sprintf("Shared basis file not found: %s", file),
      .subclass = "lna_error_basis_not_found",
      location = "read_shared_basis"
    )
  }

  h5 <- hdf5r::H5File$new(file, mode = "r")
  on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)

  # Read metadata
  meta_grp <- h5[["meta"]]
  basis_kind <- as.character(.h5_read_attr(meta_grp, "basis_kind") %||% "")[[1]]
  checksum <- as.character(.h5_read_attr(meta_grp, "checksum") %||% "")[[1]]
  mask_checksum <- as.character(.h5_read_attr(meta_grp, "mask_checksum") %||% "")[[1]]
  storage_kind <- as.character(.h5_read_attr(meta_grp, "loadings_storage") %||% "dense")[[1]]
  grid_id <- as.character(.h5_read_attr(meta_grp, "grid_id") %||% "")[[1]]
  support_id <- as.character(.h5_read_attr(meta_grp, "support_id") %||% "")[[1]]
  support_type <- as.character(.h5_read_attr(meta_grp, "support_type") %||% "")[[1]]
  spatial_units <- as.character(.h5_read_attr(meta_grp, "spatial_units") %||% "")[[1]]
  voxel_ordering <- as.character(.h5_read_attr(meta_grp, "voxel_ordering") %||% "")[[1]]
  grid_dims <- as.integer(.h5_read_attr(meta_grp, "grid_dims") %||% integer())
  safe_h5_close(meta_grp)

  # Read loadings
  root <- h5[["/"]]
  loadings <- .read_shared_basis_loadings(root, storage_kind)

  # Read mask
  mask_raw <- h5_read(root, "mask")
  mask_arr <- array(as.logical(mask_raw), dim = dim(mask_raw))

  # Read params if present
  params <- NULL
  if (root$exists("params")) {
    params_json <- h5_read(root, "params")
    if (is.matrix(params_json) || is.array(params_json)) {
      params_json <- as.vector(params_json)[[1]]
    }
    params <- tryCatch(
      jsonlite::fromJSON(params_json, simplifyVector = FALSE),
      error = function(e) NULL
    )
  }

  basis_ids <- .find_basis_ids_for_file(file)
  basis_id <- if (length(basis_ids) == 1L) basis_ids[[1]] else NULL

  # Verify checksum
  if (verify_checksum) {
    spatial_identity <- list(
      grid_id = grid_id,
      support_id = support_id,
      support_type = support_type,
      spatial_units = spatial_units,
      voxel_ordering = voxel_ordering,
      dims = grid_dims
    )
    actual <- paste0(
      "sha256:",
      .compute_artifact_checksum(
        loadings,
        mask_arr,
        basis_kind,
        params,
        spatial_identity = spatial_identity
      )
    )
    if (!identical(actual, checksum)) {
      basis_label <- basis_id %||% basename(file)
      abort_lna(
        sprintf("Checksum mismatch for basis '%s': expected %s, got %s",
                basis_label, checksum, actual),
        .subclass = "lna_error_checksum",
        location = "read_shared_basis"
      )
    }
  }

  list(
    loadings = loadings,
    mask = mask_arr,
    basis_id = basis_id,
    basis_ids = basis_ids,
    basis_kind = basis_kind,
    checksum = checksum,
    mask_checksum = mask_checksum,
    grid_id = grid_id,
    support_id = support_id,
    support_type = support_type,
    spatial_units = spatial_units,
    voxel_ordering = voxel_ordering,
    grid_dims = grid_dims,
    params = params
  )
}

#' Resolve a shared basis by identity
#'
#' Locates a shared basis artifact using a multi-step resolution chain:
#' locator path, registry lookup, configured search paths.
#'
#' @param basis_id Human-readable basis identifier.
#' @param checksum Expected checksum string (e.g., "sha256:a1b2c3...").
#' @param locator Optional relative path hint (from the archive's directory).
#' @param archive_path Optional path to the consuming `.lna.h5` file, used
#'   to resolve relative locators and find nearby registries.
#' @return Path to the verified `.lna_basis.h5` file.
#' @export
resolve_basis <- function(basis_id, checksum, locator = NULL,
                          archive_path = NULL) {
  archive_dir <- if (!is.null(archive_path)) dirname(archive_path) else NULL

  # 1. Try locator relative to archive
  if (!is.null(locator) && !is.null(archive_dir)) {
    candidate <- .resolve_locator_safe(locator, archive_dir)
    if (!is.null(candidate) && file.exists(candidate) &&
        .verify_basis_checksum(candidate, checksum)) {
      return(candidate)
    }
  }

  # 2. Walk up to find bases/registry.json
  if (!is.null(archive_dir)) {
    registry <- .find_basis_registry(archive_dir)
    if (!is.null(registry)) {
      entry <- registry$bases[[basis_id]]
      if (!is.null(entry) && identical(entry$checksum, checksum)) {
        .validate_registry_filename(entry$filename)
        registry_dir <- dirname(registry$.path)
        candidate <- normalizePath(
          file.path(registry_dir, entry$filename),
          winslash = "/", mustWork = FALSE
        )
        if (file.exists(candidate) && .verify_basis_checksum(candidate, checksum)) {
          return(candidate)
        }
      }
    }
  }

  # 3. User-configured search paths
  search_paths <- tryCatch(
    lna_options("basis_search_paths")[[1]],
    error = function(e) character(0)
  )
  if (length(search_paths) > 0) {
    hash_part <- sub("^sha256:", "", checksum)
    target_filename <- paste0("sha256-", hash_part, ".lna_basis.h5")
    for (dir in search_paths) {
      candidate <- file.path(dir, target_filename)
      if (file.exists(candidate) && .verify_basis_checksum(candidate, checksum)) {
        return(candidate)
      }
    }
  }

  # 4. Fail
  abort_lna(
    sprintf(
      "Cannot resolve basis '%s' (checksum: %s). Searched: locator=%s, registry, search_paths.",
      basis_id, checksum, locator %||% "<none>"
    ),
    .subclass = "lna_error_basis_not_found",
    location = "resolve_basis"
  )
}

#' List registered shared bases
#'
#' @param dir Directory containing `registry.json`.
#' @return A data frame with columns: basis_id, checksum, filename, kind,
#'   n_atoms, n_voxels, created. Empty data frame if no registry exists.
#' @export
list_shared_bases <- function(dir) {
  registry_path <- file.path(dir, "registry.json")
  if (!file.exists(registry_path)) {
    return(data.frame(
      basis_id = character(), checksum = character(),
      filename = character(), kind = character(),
      n_atoms = integer(), n_voxels = integer(),
      created = character(), stringsAsFactors = FALSE
    ))
  }

  registry <- jsonlite::fromJSON(registry_path, simplifyVector = FALSE)
  bases <- registry$bases %||% list()
  if (length(bases) == 0) {
    return(data.frame(
      basis_id = character(), checksum = character(),
      filename = character(), kind = character(),
      n_atoms = integer(), n_voxels = integer(),
      created = character(), stringsAsFactors = FALSE
    ))
  }

  rows <- lapply(names(bases), function(id) {
    b <- bases[[id]]
    data.frame(
      basis_id = id,
      checksum = b$checksum %||% "",
      filename = b$filename %||% "",
      kind = b$kind %||% "",
      n_atoms = as.integer(b$n_atoms %||% NA),
      n_voxels = as.integer(b$n_voxels %||% NA),
      created = b$created %||% "",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Extract loadings and mask from a template or raw list
#' @keywords internal
.extract_basis_components <- function(template) {
  # A plain raw-list template is already fully specified. Resolve it before
  # loading the optional fmrilatent namespace; this keeps the core artifact
  # path genuinely optional-dependency-free.
  if (identical(class(template), "list") &&
      !is.null(template$loadings) && !is.null(template$mask)) {
    loadings <- template$loadings
    mask_arr <- as.array(template$mask)
    kind <- template$kind %||% "unknown"
    params <- template$params %||% list()
    return(list(
      loadings = loadings,
      mask_arr = mask_arr,
      mask_obj = template$mask,
      kind = kind,
      params = params,
      grid = template$grid %||% NULL,
      support_type = template$support_type %||% NULL
    ))
  }

  if (requireNamespace("fmrilatent", quietly = TRUE)) {
    fmri_ns <- asNamespace("fmrilatent")
    is_template_fun <- if (exists("is_template", envir = fmri_ns, mode = "function")) {
      get("is_template", envir = fmri_ns, mode = "function")
    } else {
      NULL
    }
    template_loadings_fun <- if (exists("template_loadings", envir = fmri_ns, mode = "function")) {
      get("template_loadings", envir = fmri_ns, mode = "function")
    } else {
      NULL
    }
    template_mask_fun <- if (exists("template_mask", envir = fmri_ns, mode = "function")) {
      get("template_mask", envir = fmri_ns, mode = "function")
    } else {
      NULL
    }
    template_meta_fun <- if (exists("template_meta", envir = fmri_ns, mode = "function")) {
      get("template_meta", envir = fmri_ns, mode = "function")
    } else {
      NULL
    }

    # Try the template protocol
    if (is.function(is_template_fun) &&
        is.function(template_loadings_fun) &&
        is.function(template_mask_fun) &&
        isTRUE(is_template_fun(template))) {
      loadings <- template_loadings_fun(template)
      mask_obj <- template_mask_fun(template)
      mask_arr <- as.array(mask_obj)
      meta <- if (is.function(template_meta_fun)) {
        tryCatch(template_meta_fun(template), error = function(e) list())
      } else {
        list()
      }
      kind <- meta$family %||% class(template)[1]
      params <- meta
      return(list(
        loadings = loadings,
        mask_arr = mask_arr,
        mask_obj = mask_obj,
        kind = kind,
        params = params,
        grid = meta$grid %||% NULL,
        support_type = meta$support_type %||% NULL
      ))
    }
  }

  # Raw list fallback

  if (is.list(template) && !is.null(template$loadings) && !is.null(template$mask)) {
    loadings <- template$loadings
    mask_arr <- as.array(template$mask)
    kind <- template$kind %||% "unknown"
    params <- template$params %||% list()
    return(list(
      loadings = loadings,
      mask_arr = mask_arr,
      mask_obj = template$mask,
      kind = kind,
      params = params,
      grid = template$grid %||% NULL,
      support_type = template$support_type %||% NULL
    ))
  }

  abort_lna(
    "template must be a fmrilatent template object or list with $loadings and $mask",
    .subclass = "lna_error_validation",
    location = "write_shared_basis"
  )
}

.shared_basis_spatial_identity <- function(extracted, mask_arr) {
  grid <- extracted$grid
  mask_obj <- extracted$mask_obj
  if (is.null(grid) && methods::is(mask_obj, "NeuroObj")) {
    space <- neuroim2::space(mask_obj)
    grid <- list(
      dims = dim(mask_arr),
      transform = neuroim2::trans(space),
      spatial_units = "mm",
      voxel_ordering = "first-axis-fastest"
    )
  }

  if (is.null(grid)) {
    return(list(
      grid_id = "",
      support_id = "",
      support_type = "",
      spatial_units = "",
      voxel_ordering = "",
      dims = integer()
    ))
  }

  header <- list(
    dims = grid$dims %||% dim(mask_arr),
    transform = grid$transform %||% grid$affine,
    spatial_units = grid$spatial_units %||% grid$units,
    voxel_ordering = grid$voxel_ordering
  )
  identity <- .lna_grid_identity_from_header(
    header,
    units = header$spatial_units %||% "mm",
    voxel_ordering = header$voxel_ordering %||% "first-axis-fastest",
    location = "write_shared_basis:grid"
  )
  if (!identical(identity$dims, as.integer(dim(mask_arr)))) {
    abort_lna(
      "shared basis mask dimensions do not match its declared grid",
      field = "grid_dims",
      expected = identity$dims,
      actual = as.integer(dim(mask_arr)),
      .subclass = "lna_error_spatial_identity",
      location = "write_shared_basis:grid"
    )
  }

  support_type <- extracted$support_type %||% "masked"
  if (!support_type %in% c("masked", "full-grid")) {
    abort_lna(
      "support_type must be 'masked' or the explicit 'full-grid'",
      field = "support_type",
      .subclass = "lna_error_spatial_identity",
      location = "write_shared_basis:support"
    )
  }
  if (identical(support_type, "full-grid") && !all(as.logical(mask_arr))) {
    abort_lna(
      "full-grid shared bases require support at every grid voxel",
      field = "support_type",
      .subclass = "lna_error_spatial_identity",
      location = "write_shared_basis:support"
    )
  }

  list(
    grid_id = identity$grid_id,
    support_id = .lna_support_id(identity$grid_id, mask_arr),
    support_type = support_type,
    spatial_units = identity$units,
    voxel_ordering = identity$voxel_ordering,
    dims = identity$dims
  )
}

#' Compute SHA256 checksum of the immutable artifact contents
#' @keywords internal
.compute_artifact_checksum <- function(loadings, mask_arr, basis_kind, params,
                                       spatial_identity = NULL) {
  con <- rawConnection(raw(), open = "wb")
  on.exit(close(con), add = TRUE)

  has_spatial_identity <- !is.null(spatial_identity$grid_id) &&
    nzchar(spatial_identity$grid_id)
  .write_raw_string(
    con,
    if (has_spatial_identity) "lna_shared_basis_v2" else "lna_shared_basis_v1"
  )
  .write_raw_string(con, basis_kind %||% "")
  .write_raw_string(con, .encode_basis_params(params) %||% "")
  if (has_spatial_identity) {
    .write_raw_string(con, spatial_identity$grid_id)
    .write_raw_string(con, spatial_identity$support_id)
    .write_raw_string(con, spatial_identity$support_type)
    .write_raw_string(con, spatial_identity$spatial_units)
    .write_raw_string(con, spatial_identity$voxel_ordering)
    .write_raw_int(con, spatial_identity$dims)
  }
  .write_raw_mask(con, mask_arr)
  .write_raw_loadings(con, loadings)

  digest::digest(rawConnectionValue(con), algo = "sha256", serialize = FALSE)
}

.compute_mask_checksum <- function(mask_arr) {
  con <- rawConnection(raw(), open = "wb")
  on.exit(close(con), add = TRUE)

  .write_raw_string(con, "lna_shared_basis_mask_v1")
  .write_raw_mask(con, mask_arr)

  paste0(
    "sha256:",
    digest::digest(rawConnectionValue(con), algo = "sha256", serialize = FALSE)
  )
}

#' Verify checksum of a basis file without fully reading it
#' @keywords internal
.verify_basis_checksum <- function(filepath, expected_checksum) {
  tryCatch({
    basis_data <- read_shared_basis(filepath, verify_checksum = TRUE)
    identical(basis_data$checksum, expected_checksum)
  }, error = function(e) FALSE)
}

#' Read a single HDF5 attribute value
#' @keywords internal
.h5_read_attr <- function(obj, name) {
  if (!obj$attr_exists(name)) return(NULL)
  attr <- obj$attr_open(name)
  on.exit(safe_h5_close(attr), add = TRUE)
  attr$read()
}

#' Find registry.json by walking up from a directory
#' @keywords internal
.find_basis_registry <- function(start_dir, root_dir = NULL) {
  dir <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  if (!is.null(root_dir)) {
    root_dir <- normalizePath(root_dir, winslash = "/", mustWork = FALSE)
  }
  # Walk up at most 5 levels

  for (i in seq_len(5)) {
    # Check bases/registry.json
    candidate <- file.path(dir, "bases", "registry.json")
    if (file.exists(candidate)) {
      reg <- jsonlite::fromJSON(candidate, simplifyVector = FALSE)
      reg$.path <- candidate
      return(reg)
    }
    # Check registry.json directly (if we're already in bases/)
    candidate2 <- file.path(dir, "registry.json")
    if (file.exists(candidate2)) {
      reg <- jsonlite::fromJSON(candidate2, simplifyVector = FALSE)
      reg$.path <- candidate2
      return(reg)
    }
    parent <- dirname(dir)
    if (identical(parent, dir)) break
    # Do not walk above root_dir when provided
    if (!is.null(root_dir) && !startsWith(parent, root_dir)) break
    dir <- parent
  }
  NULL
}

#' Update or create registry.json
#' @keywords internal
.update_basis_registry <- function(dir, basis_id, entry) {
  registry_path <- file.path(dir, "registry.json")
  if (file.exists(registry_path)) {
    registry <- jsonlite::fromJSON(registry_path, simplifyVector = FALSE)
  } else {
    registry <- list(version = "1.0", bases = list())
  }
  registry$bases[[basis_id]] <- entry
  writeLines(
    jsonlite::toJSON(registry, auto_unbox = TRUE, pretty = TRUE),
    registry_path
  )
  invisible(registry_path)
}

.shared_basis_storage_kind <- function(loadings) {
  if (inherits(loadings, "sparseMatrix")) "dgCMatrix" else "dense"
}

.write_shared_basis_loadings <- function(root, loadings, storage_kind) {
  if (identical(storage_kind, "dgCMatrix")) {
    sparse <- methods::as(loadings, "dgCMatrix")
    h5_write_dataset(root, "loadings_i", as.integer(sparse@i))
    h5_write_dataset(root, "loadings_p", as.integer(sparse@p))
    h5_write_dataset(root, "loadings_x", as.numeric(sparse@x))
    h5_write_dataset(root, "loadings_dim", as.integer(sparse@Dim))
    return(invisible(NULL))
  }

  h5_write_dataset(root, "loadings", unname(as.matrix(loadings)))
  invisible(NULL)
}

.read_shared_basis_loadings <- function(root, storage_kind) {
  if (identical(storage_kind, "dgCMatrix")) {
    i <- as.integer(h5_read(root, "loadings_i"))
    p <- as.integer(h5_read(root, "loadings_p"))
    x <- as.numeric(h5_read(root, "loadings_x"))
    dims <- as.integer(h5_read(root, "loadings_dim"))
    return(methods::new("dgCMatrix", i = i, p = p, x = x, Dim = dims))
  }

  h5_read(root, "loadings")
}

.encode_basis_params <- function(params) {
  if (is.null(params) || length(params) == 0) return(NULL)
  jsonlite::toJSON(
    .canonicalize_hash_value(params),
    auto_unbox = TRUE,
    null = "null",
    pretty = FALSE
  )
}

.canonicalize_hash_value <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    if (!is.null(names(x))) {
      x <- x[order(names(x))]
    }
    return(lapply(x, .canonicalize_hash_value))
  }

  if (is.atomic(x) && !is.null(names(x))) {
    return(x[order(names(x))])
  }

  x
}

.write_raw_string <- function(con, value) {
  bytes <- charToRaw(enc2utf8(as.character(value %||% "")))
  writeBin(as.integer(length(bytes)), con, size = 4L, endian = "little")
  if (length(bytes) > 0) writeBin(bytes, con)
}

.write_raw_int <- function(con, value) {
  writeBin(as.integer(value), con, size = 4L, endian = "little")
}

.write_raw_double <- function(con, value) {
  writeBin(as.double(value), con, size = 8L, endian = "little")
}

.write_raw_mask <- function(con, mask_arr) {
  .write_raw_string(con, "mask")
  .write_raw_int(con, dim(mask_arr))
  writeBin(as.raw(as.integer(as.logical(mask_arr))), con)
}

.write_raw_loadings <- function(con, loadings) {
  if (inherits(loadings, "sparseMatrix")) {
    sparse <- methods::as(loadings, "dgCMatrix")
    .write_raw_string(con, "dgCMatrix")
    .write_raw_int(con, sparse@Dim)
    .write_raw_int(con, sparse@i)
    .write_raw_int(con, sparse@p)
    .write_raw_double(con, sparse@x)
    return(invisible(NULL))
  }

  dense <- unname(as.matrix(loadings))
  .write_raw_string(con, "dense")
  .write_raw_int(con, dim(dense))
  .write_raw_double(con, as.vector(dense))
  invisible(NULL)
}

#' Validate that a basis_id is a safe, well-formed identifier
#' @keywords internal
.validate_basis_id <- function(basis_id) {
  if (!is.character(basis_id) || length(basis_id) != 1L || !nzchar(basis_id)) {
    abort_lna(
      "basis_id must be a non-empty character string",
      .subclass = "lna_error_validation",
      location = "write_shared_basis:basis_id"
    )
  }
  if (nchar(basis_id) > 256L) {
    abort_lna(
      "basis_id must be at most 256 characters",
      .subclass = "lna_error_validation",
      location = "write_shared_basis:basis_id"
    )
  }
  if (grepl("[[:cntrl:]]", basis_id)) {
    abort_lna(
      "basis_id must not contain control characters",
      .subclass = "lna_error_validation",
      location = "write_shared_basis:basis_id"
    )
  }
  if (grepl("[/\\\\]", basis_id)) {
    abort_lna(
      "basis_id must not contain path separators",
      .subclass = "lna_error_validation",
      location = "write_shared_basis:basis_id"
    )
  }
  invisible(TRUE)
}

#' Resolve a locator path safely, rejecting path traversal
#' @keywords internal
.resolve_locator_safe <- function(locator, archive_dir) {
  # Reject absolute paths
  if (grepl("^[/\\\\]", locator) || grepl("^[A-Za-z]:", locator)) {
    return(NULL)
  }
  candidate <- normalizePath(
    file.path(archive_dir, locator),
    winslash = "/", mustWork = FALSE
  )
  # Verify the resolved path stays within an ancestor of archive_dir.
  # Count how many ".." components the locator has and ensure the resolved

  # path starts with a corresponding ancestor of archive_dir.
  norm_archive <- normalizePath(archive_dir, winslash = "/", mustWork = FALSE)
  # Walk up from archive_dir by the number of ".." segments in locator,
  # plus one extra level of tolerance (the archive is typically nested
  # inside a dataset root).
  dotdot_count <- sum(unlist(strsplit(locator, "[/\\\\]")) == "..")
  anchor <- norm_archive
  for (i in seq_len(min(dotdot_count + 1L, 10L))) {
    anchor <- dirname(anchor)
  }
  if (!startsWith(candidate, anchor)) {
    return(NULL)
  }
  candidate
}

#' Validate that a registry filename matches the content-addressed pattern
#' @keywords internal
.validate_registry_filename <- function(filename) {
  if (grepl("\\.\\.", filename) || grepl("[/\\\\]", filename)) {
    abort_lna(
      sprintf("Registry filename contains path traversal characters: %s", filename),
      .subclass = "lna_error_validation",
      location = "resolve_basis:registry_filename"
    )
  }
  if (!grepl("^sha256-[0-9a-f]{64}\\.lna_basis\\.h5$", filename)) {
    abort_lna(
      sprintf("Registry filename does not match expected content-addressed pattern: %s", filename),
      .subclass = "lna_error_validation",
      location = "resolve_basis:registry_filename"
    )
  }
  invisible(TRUE)
}

.find_basis_ids_for_file <- function(file) {
  registry <- .find_basis_registry(dirname(file))
  if (is.null(registry) || is.null(registry$bases)) {
    return(character())
  }

  filename <- basename(file)
  ids <- names(Filter(function(entry) identical(entry$filename %||% "", filename), registry$bases))
  as.character(ids)
}
