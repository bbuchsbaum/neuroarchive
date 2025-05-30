#' Discover and Validate Transform Descriptors in HDF5 Group
#'
#' @description Lists objects within the `/transforms` group, parses names
#'   following the `NN_type.json` pattern, validates the sequence, and returns
#'   metadata about the transforms.
#'
#' @param h5_group An `H5Group` object representing the `/transforms` group.
#'
#' @return A `tibble::tibble` with columns:
#'   * `name` (character): The full object name (e.g., "00_mask.json").
#'   * `type` (character): The transform type extracted from the name (e.g., "mask").
#'   * `index` (integer): The zero-based index extracted from the name (e.g., 0).
#'   Returns an empty tibble if the group is empty.
#'
#' @details
#'   Invalid descriptor names or indices trigger errors. Sequence validation
#'   uses `abort_lna()` and signals the subclass `lna_error_sequence` when the
#'   numeric indices are not contiguous starting from zero. Other malformed
#'   names currently raise standard errors via `stop()`.
#'
#' @import hdf5r
#' @importFrom tibble tibble
#' @keywords internal
discover_transforms <- function(h5_group) {
  assert_h5group(h5_group)

  obj_names <- tryCatch({
    names(h5_group)
  }, error = function(e) {

    print("Error occurred during h5_group$names():")
    print(conditionMessage(e))
    abort_lna(
      "Failed to list names in HDF5 group.",
      .subclass = "lna_error_io",
      location = "discover_transforms",
      parent = e
    )
  })

  # Filter out report files (which end with _report.json) to avoid pattern match conflicts
  obj_names <- obj_names[!grepl("_report\\.json$", obj_names)]

  # Handle empty group (after filtering)
  if (length(obj_names) == 0) {
    return(tibble::tibble(name = character(), type = character(), index = integer()))
  }

  # Regex to capture NN, type from NN_type.json
  # NN must be digits, type must be characters except '_'
  pattern <- "^(\\d+)_([A-Za-z0-9._]+)\\.json$"
  matches <- regexec(pattern, obj_names)

  extracted_data <- lapply(seq_along(matches), function(i) {
    match_info <- matches[[i]]
    if (match_info[1] == -1) { # No match for this name
      # Consider warning or error? Spec says ensure NN_ is contiguous,
      # implies non-matching names might be an error or ignored.
      # Let's error for now if non-matching names are found.
      abort_lna(
        paste0(
          "Invalid object name found in /transforms: ",
          obj_names[i],
          ". Expected format NN_type.json."
        ),
        .subclass = "lna_error_descriptor",
        location = "discover_transforms"
      )
      # return(NULL) # Alternative: ignore non-matching files
    }
    full_name <- obj_names[i]
    index_str <- substr(full_name, match_info[2], match_info[2] + attr(match_info, "match.length")[2] - 1)
    type_str  <- substr(full_name, match_info[3], match_info[3] + attr(match_info, "match.length")[3] - 1)

    # Convert index, handle potential non-integer strings caught by regex (though unlikely)
    index_int <- suppressWarnings(as.integer(index_str))
    if (is.na(index_int)) {
      abort_lna(
        paste0(
          "Invalid numeric index found in transform name: ",
          full_name
        ),
        .subclass = "lna_error_descriptor",
        location = "discover_transforms"
      )
    }

    list(name = full_name, type = type_str, index = index_int)
  })

  # Filter out NULLs if we chose to ignore non-matching names previously
  # extracted_data <- extracted_data[!sapply(extracted_data, is.null)]

  if (length(extracted_data) == 0) {
      # This case occurs if obj_names was not empty but nothing matched the pattern
      # and we chose to ignore non-matching names. Should probably error if we expect
      # all names to match.
       abort_lna(
         "No valid transform descriptors (NN_type.json) found in non-empty /transforms group.",
         .subclass = "lna_error_descriptor",
         location = "discover_transforms"
       )
      # return(tibble::tibble(name = character(), type = character(), index = integer()))
  }


  # Combine into a temporary data frame or tibble for sorting
  temp_df <- do.call(rbind.data.frame, c(extracted_data, stringsAsFactors = FALSE))
  # Ensure correct types after rbind
  temp_df$index <- as.integer(temp_df$index)

  # Sort by index
  sorted_df <- temp_df[order(temp_df$index), , drop = FALSE]

  # Validate sequence: indices must be 0, 1, 2, ..., n-1
  expected_indices <- seq(0, nrow(sorted_df) - 1)
  if (!identical(sorted_df$index, expected_indices)) {
    abort_lna(
      paste0(
        "Transform descriptor indices are not contiguous starting from 0. Found indices: ",
        paste(sorted_df$index, collapse = ", ")
      ),
      .subclass = "lna_error_sequence",
      location = "discover_transforms"
    )
  }

  # Convert to final tibble
  result_tibble <- tibble::tibble(
    name = as.character(sorted_df$name),
    type = as.character(sorted_df$type),
    index = as.integer(sorted_df$index)
  )

  return(result_tibble)
} 