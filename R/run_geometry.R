# Run-scoped spatial geometry for low-level multi-run LNA archives.

.lna_run_geometry_schema_version <- "1.0.0"
.lna_spatial_header_fields <- c(
  "geometry_schema_version", "run_id", "dims", "spacing", "origin",
  "transform", "affine", "spatial_units", "units", "voxel_ordering",
  "n_volumes", "grid_id"
)

.run_geometry_abort <- function(message, field = "run_geometry", runs = NULL) {
  abort_lna(
    message,
    field = field,
    runs = runs,
    .subclass = "lna_error_run_geometry",
    location = paste0("run_geometry:", field)
  )
}

.lna_input_runs <- function(x, run_ids) {
  values <- if (is.list(x) && !methods::is(x, "NeuroObj")) x else list(x)
  if (length(values) != length(run_ids)) {
    .run_geometry_abort(
      "number of input runs does not match run_ids",
      "run_ids",
      run_ids
    )
  }
  names(values) <- run_ids
  values
}

.lna_header_has_geometry <- function(x) {
  is.list(x) && length(intersect(
    names(x), c("dims", "transform", "affine", "grid_id")
  )) > 0L
}

.lna_extract_geometry_fields <- function(x) {
  x[intersect(names(x), .lna_spatial_header_fields)]
}

.run_geometry_vector <- function(x, length_required, field, positive = FALSE) {
  if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
  if (!is.numeric(x) || length(x) != length_required || any(!is.finite(x))) {
    .run_geometry_abort(
      paste0(field, " must contain ", length_required, " finite numeric values"),
      field
    )
  }
  x <- as.double(x)
  if (positive && any(x <= 0)) {
    .run_geometry_abort(paste0(field, " must be strictly positive"), field)
  }
  x
}

.run_geometry_scalar_string <- function(x, field) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    .run_geometry_abort(paste0(field, " must be one non-empty string"), field)
  }
  enc2utf8(x)
}

.run_geometry_nvolumes <- function(x, input, field) {
  if (is.null(x)) {
    input_dims <- dim(input)
    return(if (length(input_dims) >= 4L) as.integer(input_dims[[4L]]) else 1L)
  }
  if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
  value <- suppressWarnings(as.integer(x))
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != value || value < 1L) {
    .run_geometry_abort(paste0(field, " must be a positive integer"), field)
  }
  value
}

.lna_normalize_run_geometry <- function(geometry, run_id, input = NULL) {
  if (!is.list(geometry) || is.data.frame(geometry) || !length(geometry) ||
      is.null(names(geometry)) || any(!nzchar(names(geometry)))) {
    .run_geometry_abort(
      paste0("geometry for ", run_id, " must be a non-empty named list"),
      paste0(run_id, ".geometry"),
      run_id
    )
  }
  dims <- tryCatch(
    .lna_normalize_grid_dims(geometry$dims),
    lna_error_spatial_identity = function(e) {
      .run_geometry_abort(
        conditionMessage(e), paste0(run_id, ".dims"), run_id
      )
    }
  )
  affine_value <- geometry$transform %||% geometry$affine
  transform <- tryCatch(
    .lna_normalize_affine(affine_value),
    lna_error_spatial_identity = function(e) {
      .run_geometry_abort(
        conditionMessage(e), paste0(run_id, ".transform"), run_id
      )
    }
  )

  input_dims <- dim(input)
  if (!is.null(input) && (is.null(input_dims) || length(input_dims) < 3L)) {
    .run_geometry_abort(
      paste0("input ", run_id, " has no three-dimensional spatial domain"),
      paste0(run_id, ".input_dims"),
      run_id
    )
  }
  if (!is.null(input) && !identical(as.integer(input_dims[1:3]), dims)) {
    .run_geometry_abort(
      paste0(
        "geometry dims for ", run_id, " do not match input dims: expected ",
        paste(input_dims[1:3], collapse = "x"), " but received ",
        paste(dims, collapse = "x")
      ),
      paste0(run_id, ".dims"),
      run_id
    )
  }

  derived_spacing <- sqrt(colSums(transform[1:3, 1:3, drop = FALSE]^2))
  spacing <- if (is.null(geometry$spacing)) {
    derived_spacing
  } else {
    .run_geometry_vector(
      geometry$spacing, 3L, paste0(run_id, ".spacing"), positive = TRUE
    )
  }
  origin <- if (is.null(geometry$origin)) {
    as.double(transform[1:3, 4L])
  } else {
    .run_geometry_vector(geometry$origin, 3L, paste0(run_id, ".origin"))
  }
  spatial_units <- .run_geometry_scalar_string(
    geometry$spatial_units %||% geometry$units %||% "mm",
    paste0(run_id, ".spatial_units")
  )
  voxel_ordering <- .run_geometry_scalar_string(
    geometry$voxel_ordering %||% "first-axis-fastest",
    paste0(run_id, ".voxel_ordering")
  )
  grid_id <- .lna_grid_id(dims, transform, spatial_units, voxel_ordering)
  if (!is.null(geometry$grid_id)) {
    supplied <- tolower(.run_geometry_scalar_string(
      geometry$grid_id, paste0(run_id, ".grid_id")
    ))
    if (!grepl("^sha256:[0-9a-f]{64}$", supplied) ||
        !identical(supplied, grid_id)) {
      .run_geometry_abort(
        paste0("stored grid_id for ", run_id, " does not match its geometry"),
        paste0(run_id, ".grid_id"),
        run_id
      )
    }
  }

  list(
    geometry_schema_version = .lna_run_geometry_schema_version,
    run_id = run_id,
    dims = dims,
    spacing = spacing,
    origin = origin,
    transform = transform,
    spatial_units = spatial_units,
    voxel_ordering = voxel_ordering,
    n_volumes = .run_geometry_nvolumes(
      geometry$n_volumes, input, paste0(run_id, ".n_volumes")
    ),
    grid_id = grid_id
  )
}

.lna_derive_one_run_geometry <- function(input, run_id) {
  if (!methods::is(input, "NeuroObj")) return(NULL)
  header <- derive_header_from_input(input)
  if (is.null(header)) {
    .run_geometry_abort(
      paste0("could not derive NeuroSpace geometry for ", run_id),
      paste0(run_id, ".geometry"),
      run_id
    )
  }
  header
}

.lna_explicit_geometry_map <- function(run_geometry, run_ids) {
  if (length(run_ids) == 1L && .lna_header_has_geometry(run_geometry)) {
    return(stats::setNames(list(run_geometry), run_ids))
  }
  if (!is.list(run_geometry) || is.data.frame(run_geometry) ||
      is.null(names(run_geometry)) || any(!nzchar(names(run_geometry))) ||
      anyDuplicated(names(run_geometry))) {
    .run_geometry_abort(
      "multi-run geometry must be a uniquely named list keyed by run_id",
      "run_geometry",
      run_ids
    )
  }
  missing <- setdiff(run_ids, names(run_geometry))
  extra <- setdiff(names(run_geometry), run_ids)
  if (length(missing) || length(extra)) {
    detail <- c(
      if (length(missing)) paste0("missing: ", paste(missing, collapse = ", ")),
      if (length(extra)) paste0("unknown: ", paste(extra, collapse = ", "))
    )
    .run_geometry_abort(
      paste0("run_geometry keys do not match run_ids (", paste(detail, collapse = "; "), ")"),
      "run_geometry",
      unique(c(missing, extra))
    )
  }
  run_geometry[run_ids]
}

.lna_resolve_run_geometry <- function(
    x,
    run_ids,
    run_geometry = NULL,
    header = list(),
    geometry_policy = c("optional", "fmri")) {
  geometry_policy <- match.arg(geometry_policy)
  inputs <- .lna_input_runs(x, run_ids)

  if (is.null(run_geometry)) {
    derived <- lapply(run_ids, function(run_id) {
      .lna_derive_one_run_geometry(inputs[[run_id]], run_id)
    })
    names(derived) <- run_ids
    have <- !vapply(derived, is.null, logical(1))
    if (all(have)) {
      raw_map <- derived
    } else if (any(have)) {
      missing <- run_ids[!have]
      .run_geometry_abort(
        paste0(
          "geometry was derivable for only some runs; supply run_geometry for: ",
          paste(missing, collapse = ", ")
        ),
        "run_geometry",
        missing
      )
    } else if (.lna_header_has_geometry(header)) {
      if (length(run_ids) != 1L) {
        .run_geometry_abort(
          paste0(
            "one global spatial header is ambiguous for multiple runs: ",
            paste(run_ids, collapse = ", "),
            "; supply run_geometry keyed by run_id"
          ),
          "header",
          run_ids
        )
      }
      raw_map <- stats::setNames(list(.lna_extract_geometry_fields(header)), run_ids)
    } else if (identical(geometry_policy, "fmri")) {
      .run_geometry_abort(
        paste0(
          "fMRI geometry is required for every run; missing: ",
          paste(run_ids, collapse = ", ")
        ),
        "run_geometry",
        run_ids
      )
    } else {
      return(list())
    }
  } else {
    raw_map <- .lna_explicit_geometry_map(run_geometry, run_ids)
  }

  normalized <- lapply(run_ids, function(run_id) {
    .lna_normalize_run_geometry(
      raw_map[[run_id]], run_id, input = inputs[[run_id]]
    )
  })
  names(normalized) <- run_ids

  if (length(run_ids) > 1L && .lna_header_has_geometry(header)) {
    global <- .lna_normalize_run_geometry(
      .lna_extract_geometry_fields(header), run_ids[[1L]], input = inputs[[1L]]
    )
    grids <- vapply(normalized, `[[`, character(1), "grid_id")
    if (length(unique(grids)) != 1L || !identical(global$grid_id, grids[[1L]])) {
      .run_geometry_abort(
        paste0(
          "global spatial header conflicts with run-scoped geometry for: ",
          paste(run_ids, collapse = ", ")
        ),
        "header",
        run_ids
      )
    }
  }
  normalized
}

.lna_global_header_for_write <- function(header, run_geometry) {
  non_spatial <- header[setdiff(names(header), .lna_spatial_header_fields)]
  if (!length(run_geometry)) return(header)
  grids <- vapply(run_geometry, `[[`, character(1), "grid_id")
  if (length(unique(grids)) == 1L) {
    compatible <- run_geometry[[1L]]
    if (length(run_geometry) > 1L) {
      compatible$run_id <- NULL
      compatible$n_volumes <- NULL
    }
    return(utils::modifyList(non_spatial, compatible))
  }
  non_spatial
}

.lna_assert_aggregate_geometry <- function(x, run_ids, run_geometry, transforms) {
  if (!length(transforms) || !identical(transforms[[1L]], "aggregate_runs") ||
      length(run_ids) < 2L) {
    return(invisible(TRUE))
  }
  inputs <- .lna_input_runs(x, run_ids)
  dims <- lapply(inputs, function(value) as.integer(dim(value)[1:3]))
  reference_dims <- dims[[1L]]
  incompatible_dims <- run_ids[!vapply(
    dims, identical, logical(1), reference_dims
  )]
  if (length(incompatible_dims)) {
    .run_geometry_abort(
      paste0(
        "aggregate_runs cannot combine different spatial dimensions; ",
        "incompatible runs: ", paste(incompatible_dims, collapse = ", ")
      ),
      "aggregate_runs.dims",
      incompatible_dims
    )
  }
  if (length(run_geometry)) {
    grids <- vapply(run_geometry, `[[`, character(1), "grid_id")
    incompatible_grids <- run_ids[grids != grids[[1L]]]
    if (length(incompatible_grids)) {
      .run_geometry_abort(
        paste0(
          "aggregate_runs cannot combine different spatial grids; ",
          "incompatible runs: ", paste(incompatible_grids, collapse = ", ")
        ),
        "aggregate_runs.grid_id",
        incompatible_grids
      )
    }
  }
  invisible(TRUE)
}
