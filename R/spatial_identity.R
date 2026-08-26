# Canonical spatial identity primitives
#
# These identifiers are pure functions of normalized values. Their byte
# encodings are explicitly versioned and little-endian so they do not depend
# on R serialization, object attributes, or host byte order.

.lna_identity_sha256 <- function(bytes) {
  paste0(
    "sha256:",
    digest::digest(bytes, algo = "sha256", serialize = FALSE)
  )
}

.lna_is_canonical_sha256 <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) &&
    grepl("^sha256:[0-9a-f]{64}$", value)
}

.lna_canonical_sha256 <- function(value, field, allow_null = FALSE,
                                  canonicalize = FALSE,
                                  location = ".lna_canonical_sha256",
                                  subclass = "lna_error_checksum") {
  if (is.null(value) && isTRUE(allow_null)) return(NULL)
  accepted <- is.character(value) && length(value) == 1L && !is.na(value) &&
    grepl("^sha256:[0-9a-fA-F]{64}$", value)
  if (!accepted) {
    abort_lna(
      paste0(field, " must be sha256:<64 lowercase hex digits>"),
      field = field,
      code = "noncanonical_sha256",
      .subclass = subclass,
      location = location
    )
  }
  normalized <- tolower(value)
  if (!isTRUE(canonicalize) && !identical(value, normalized)) {
    abort_lna(
      paste0(field, " must use canonical lowercase hexadecimal"),
      field = field,
      code = "noncanonical_sha256",
      .subclass = subclass,
      location = location
    )
  }
  normalized
}

.lna_identity_write_string <- function(connection, value) {
  value <- enc2utf8(value)
  bytes <- charToRaw(value)
  writeBin(as.integer(length(bytes)), connection, size = 4L, endian = "little")
  if (length(bytes) > 0L) writeBin(bytes, connection)
  invisible(NULL)
}

.lna_identity_scalar_string <- function(value, field) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value)) {
    abort_lna(
      sprintf("%s must be one non-empty string", field),
      field = field,
      .subclass = "lna_error_spatial_identity",
      location = ".lna_spatial_identity"
    )
  }
  enc2utf8(value)
}

.lna_numeric_id <- function(value, domain) {
  domain <- .lna_identity_scalar_string(domain, "domain")
  if (!is.numeric(value) || any(!is.finite(value))) {
    abort_lna(
      "canonical numeric identity requires finite numeric data",
      field = "value",
      .subclass = "lna_error_spatial_identity",
      location = ".lna_numeric_id"
    )
  }
  connection <- rawConnection(raw(), "wb")
  on.exit(close(connection), add = TRUE)
  writeBin(as.double(value), connection, size = 8L, endian = "little")
  .lna_identity_sha256(
    c(charToRaw(domain), rawConnectionValue(connection))
  )
}

.lna_normalize_affine <- function(affine) {
  if (is.data.frame(affine)) affine <- as.matrix(affine)
  if (!is.numeric(affine) || length(affine) != 16L ||
      any(!is.finite(affine))) {
    abort_lna(
      "affine must contain exactly 16 finite numeric values",
      field = "affine",
      .subclass = "lna_error_spatial_identity",
      location = ".lna_affine_id"
    )
  }
  matrix(as.double(affine), nrow = 4L, ncol = 4L)
}

.lna_normalize_grid_dims <- function(dims) {
  if (!is.numeric(dims) || length(dims) != 3L || any(!is.finite(dims))) {
    abort_lna(
      "dims must contain exactly three finite numeric values",
      field = "dims",
      .subclass = "lna_error_spatial_identity",
      location = ".lna_grid_id"
    )
  }
  normalized <- suppressWarnings(as.integer(dims))
  if (any(dims != normalized) || any(normalized < 1L)) {
    abort_lna(
      "dims must contain exactly three positive integers",
      field = "dims",
      .subclass = "lna_error_spatial_identity",
      location = ".lna_grid_id"
    )
  }
  normalized
}

.lna_normalize_support <- function(support) {
  if (methods::is(support, "NeuroVol")) support <- as.array(support)
  if (!is.atomic(support) || is.character(support) || is.complex(support) ||
      anyNA(support)) {
    abort_lna(
      "support must be finite logical or binary numeric data",
      field = "support",
      .subclass = "lna_error_spatial_identity",
      location = ".lna_support_id"
    )
  }
  if (is.numeric(support) && any(!is.finite(support))) {
    abort_lna(
      "support must be finite logical or binary numeric data",
      field = "support",
      .subclass = "lna_error_spatial_identity",
      location = ".lna_support_id"
    )
  }
  if (any(!support %in% c(FALSE, TRUE, 0, 1))) {
    abort_lna(
      "numeric support values must be exactly zero or one",
      field = "support",
      .subclass = "lna_error_spatial_identity",
      location = ".lna_support_id"
    )
  }
  as.logical(as.vector(support))
}

#' Compute the canonical affine identity
#'
#' The `lna_affine_v1` byte domain intentionally preserves the checksum used
#' by the first parcel-map guard, while centralizing its implementation.
#'
#' @param affine A finite numeric 4 by 4 affine, or an equivalent length-16
#'   numeric container in column-major order.
#' @return A `sha256:<hex>` identity string.
#' @keywords internal
.lna_affine_id <- function(affine) {
  normalized <- .lna_normalize_affine(affine)
  .lna_numeric_id(normalized, "lna_affine_v1")
}

#' Compute the canonical spatial-grid identity
#'
#' Grid identity schema `lna_grid_v1` covers the three spatial dimensions,
#' full affine identity, spatial units, and voxel-ordering convention.
#'
#' @param dims Three positive integer spatial dimensions.
#' @param affine A finite 4 by 4 affine.
#' @param units One non-empty spatial-units string.
#' @param voxel_ordering One non-empty voxel-ordering convention string.
#' @return A `sha256:<hex>` identity string.
#' @keywords internal
.lna_grid_id <- function(dims, affine, units, voxel_ordering) {
  dims <- .lna_normalize_grid_dims(dims)
  affine_id <- .lna_affine_id(affine)
  units <- .lna_identity_scalar_string(units, "units")
  voxel_ordering <- .lna_identity_scalar_string(
    voxel_ordering,
    "voxel_ordering"
  )

  connection <- rawConnection(raw(), "wb")
  on.exit(close(connection), add = TRUE)
  .lna_identity_write_string(connection, "lna_grid_v1")
  .lna_identity_write_string(connection, "dims")
  writeBin(dims, connection, size = 4L, endian = "little")
  .lna_identity_write_string(connection, "affine_id")
  .lna_identity_write_string(connection, affine_id)
  .lna_identity_write_string(connection, "units")
  .lna_identity_write_string(connection, units)
  .lna_identity_write_string(connection, "voxel_ordering")
  .lna_identity_write_string(connection, voxel_ordering)
  .lna_identity_sha256(rawConnectionValue(connection))
}

#' Compute the canonical support identity for a grid
#'
#' Support identity schema `lna_support_v1` hashes the canonical grid identity
#' and the support bits in first-to-last R vector order. Array dimensions and
#' other container attributes are deliberately ignored because shape and voxel
#' order are properties of `grid_id`.
#'
#' @param grid_id A canonical `sha256:<hex>` grid identity.
#' @param support Logical or binary numeric support values.
#' @return A `sha256:<hex>` identity string.
#' @keywords internal
.lna_support_id <- function(grid_id, support) {
  grid_id <- .lna_canonical_sha256(
    grid_id, "grid_id", location = ".lna_support_id",
    subclass = "lna_error_spatial_identity"
  )
  support <- .lna_normalize_support(support)

  connection <- rawConnection(raw(), "wb")
  on.exit(close(connection), add = TRUE)
  .lna_identity_write_string(connection, "lna_support_v1")
  .lna_identity_write_string(connection, "grid_id")
  .lna_identity_write_string(connection, grid_id)
  .lna_identity_write_string(connection, "support_length")
  writeBin(as.integer(length(support)), connection, size = 4L, endian = "little")
  .lna_identity_write_string(connection, "support_bits")
  if (length(support) > 0L) {
    writeBin(as.raw(as.integer(support)), connection)
  }
  .lna_identity_sha256(rawConnectionValue(connection))
}

#' Compute canonical LNA spatial identities
#'
#' These public constructors expose the versioned byte identities used by LNA
#' geometry guards. A grid identity covers spatial dimensions, the full affine,
#' spatial units, and voxel ordering. A support identity additionally covers a
#' binary voxel support on that grid.
#'
#' @param affine A finite numeric 4 by 4 affine.
#' @param dims Three positive integer spatial dimensions.
#' @param units One non-empty spatial-units string.
#' @param voxel_ordering One non-empty voxel-ordering convention string.
#' @param grid_id A value returned by `lna_grid_identity()`.
#' @param support Logical or binary numeric support values.
#' @return A `sha256:<hex>` identity string.
#' @export
lna_affine_identity <- function(affine) {
  .lna_affine_id(affine)
}

#' @rdname lna_affine_identity
#' @export
lna_grid_identity <- function(
    dims,
    affine,
    units = "mm",
    voxel_ordering = "first-axis-fastest") {
  .lna_grid_id(dims, affine, units, voxel_ordering)
}

#' @rdname lna_affine_identity
#' @export
lna_support_identity <- function(grid_id, support) {
  .lna_support_id(grid_id, support)
}

.lna_grid_identity_from_header <- function(
    header,
    units = NULL,
    voxel_ordering = NULL,
    location = ".lna_grid_identity_from_header") {
  if (is.null(header$dims) || length(header$dims) < 3L) {
    abort_lna(
      "spatial header is missing grid dimensions",
      field = "dims",
      .subclass = "lna_error_spatial_identity",
      location = location
    )
  }
  if (is.null(header$transform)) {
    abort_lna(
      "spatial header is missing its affine transform",
      field = "transform",
      .subclass = "lna_error_spatial_identity",
      location = location
    )
  }
  units <- units %||% header$spatial_units %||% header$units %||% "mm"
  voxel_ordering <- voxel_ordering %||%
    header$voxel_ordering %||% "first-axis-fastest"
  dims <- .lna_normalize_grid_dims(header$dims[1:3])
  affine <- .lna_normalize_affine(header$transform)
  list(
    grid_id = .lna_grid_id(dims, affine, units, voxel_ordering),
    dims = dims,
    affine_id = .lna_affine_id(affine),
    units = .lna_identity_scalar_string(units, "units"),
    voxel_ordering = .lna_identity_scalar_string(
      voxel_ordering,
      "voxel_ordering"
    )
  )
}
