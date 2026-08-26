# Installed-package acceptance court for the Neuroarchive BIDS exchange profile.
#
# Usage:
#   Rscript bids-exchange-installed-court.R /absolute/output/directory
#
# The caller must set NEUROARCHIVE_TARBALL_SHA256 to the exact installed source
# tarball identity. This script uses exported package APIs only.

suppressPackageStartupMessages({
  library(neuroarchive)
  library(neuroim2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !nzchar(args[[1L]])) {
  stop("one output directory argument is required", call. = FALSE)
}
court_root <- normalizePath(args[[1L]], winslash = "/", mustWork = FALSE)
if (file.exists(court_root) || dir.exists(court_root)) {
  stop("output directory must not already exist", call. = FALSE)
}
if (!dir.create(court_root, recursive = TRUE, showWarnings = FALSE)) {
  stop("could not create court output directory", call. = FALSE)
}
tarball_hash <- Sys.getenv("NEUROARCHIVE_TARBALL_SHA256", unset = "")
if (!grepl("^sha256:[0-9a-f]{64}$", tarball_hash)) {
  stop("NEUROARCHIVE_TARBALL_SHA256 must identify the installed tarball", call. = FALSE)
}

source_root <- file.path(court_root, "source-bids")
archive_root <- file.path(court_root, "archive")
receiver_root <- file.path(court_root, "receiver")
restored_root <- file.path(court_root, "restored-bids")
dir.create(source_root, recursive = TRUE)

write_source_record <- function(
    record_id,
    relative,
    entities,
    values,
    scope = "raw",
    pipeline = NULL,
    reference_assertion = NULL,
    with_events = FALSE) {
  source <- file.path(source_root, relative)
  dir.create(dirname(source), recursive = TRUE, showWarnings = FALSE)
  image <- neuroim2::NeuroVec(
    values,
    neuroim2::NeuroSpace(dim(values), trans = diag(c(2, 2, 2, 1)))
  )
  neuroim2::write_vec(image, source)
  metadata_relative <- sub("[.]nii[.]gz$", ".json", relative)
  metadata_source <- file.path(source_root, metadata_relative)
  metadata <- list(
    RepetitionTime = 0.8,
    TaskName = entities$task,
    SliceTiming = seq(0, 0.7, length.out = dim(values)[[3L]])
  )
  jsonlite::write_json(
    metadata, metadata_source, auto_unbox = TRUE, null = "null", pretty = TRUE
  )
  associations <- list(list(
    association_id = paste0("metadata:", record_id),
    role = "metadata-source",
    path = metadata_relative,
    checksum = lna_file_checksum(metadata_source),
    metadata = list(
      inheritance_level = "run", precedence = 100L,
      fields = c("RepetitionTime", "TaskName", "SliceTiming"),
      media_type = "application/json"
    )
  ))
  event_checksum <- NULL
  event_relative <- NULL
  if (isTRUE(with_events)) {
    event_entities <- entities[
      intersect(names(entities), c("sub", "ses", "task", "run"))
    ]
    event_stem <- paste(vapply(names(event_entities), function(name) {
      paste0(name, "-", event_entities[[name]])
    }, character(1)), collapse = "_")
    event_relative <- file.path(
      dirname(relative), paste0(event_stem, "_events.tsv")
    )
    event_source <- file.path(source_root, event_relative)
    writeLines(
      c("onset\tduration\ttrial_type", "0\t1\tgo", "2\t1\tstop"),
      event_source, useBytes = TRUE
    )
    event_checksum <- lna_file_checksum(event_source)
    associations[[length(associations) + 1L]] <- list(
      association_id = paste0("events:", record_id),
      role = "events",
      path = gsub("\\\\", "/", event_relative),
      checksum = event_checksum,
      metadata = list(media_type = "text/tab-separated-values")
    )
  }
  record <- lna_source_record(
    record_id = record_id,
    source_path = gsub("\\\\", "/", relative),
    source_checksum = lna_file_checksum(source),
    entities = entities,
    datatype = "func",
    suffix = "bold",
    extension = ".nii.gz",
    scope = scope,
    pipeline = pipeline,
    resolved_metadata = metadata,
    associations = associations,
    geometry = list(
      dims = dim(values)[1:3],
      affine = diag(c(2, 2, 2, 1)),
      units = "mm",
      voxel_ordering = "first-axis-fastest",
      n_volumes = dim(values)[[4L]]
    ),
    reference_assertion = reference_assertion
  )
  list(
    record = record,
    values = values,
    event_relative = event_relative,
    event_checksum = event_checksum
  )
}

base_values <- array(
  as.numeric((seq_len(3072) * 7919L) %% 10007L),
  dim = c(8L, 8L, 8L, 6L)
)
raw_one <- write_source_record(
  "raw-rest",
  "sub-01/ses-A/func/sub-01_ses-A_task-rest_acq-mb_run-01_bold.nii.gz",
  list(sub = "01", ses = "A", task = "rest", acq = "mb", run = "01"),
  base_values,
  with_events = TRUE
)
raw_two <- write_source_record(
  "raw-nback",
  "sub-02/func/sub-02_task-nback_run-01_bold.nii.gz",
  list(sub = "02", task = "nback", run = "01"),
  base_values + 100,
  with_events = TRUE
)
derivative_id <- "derivative-rest"
derivative <- write_source_record(
  derivative_id,
  paste0(
    "derivatives/fmriprep/sub-01/ses-A/func/",
    "sub-01_ses-A_task-rest_run-01_space-MNI152NLin2009cAsym_",
    "res-2_desc-preproc_bold.nii.gz"
  ),
  list(
    sub = "01", ses = "A", task = "rest", run = "01",
    space = "MNI152NLin2009cAsym", res = "2", desc = "preproc"
  ),
  base_values + 200,
  scope = "derivative",
  pipeline = "fmriprep",
  reference_assertion = lna_reference_assertion(
    record_id = derivative_id,
    identifier = "MNI152NLin2009cAsym",
    kind = "template",
    origin = "bids-space-entity",
    evidence = list(space_entity = "MNI152NLin2009cAsym")
  )
)
initial <- list(raw_one, raw_two, derivative)
initial_records <- lapply(initial, `[[`, "record")

packed <- write_lna_dataset(
  initial_records,
  archive_root,
  source_root,
  metadata_policy = lna_metadata_policy(portable = TRUE),
  owner_id = "installed-bids-court-pack"
)
stopifnot(
  length(packed$written_record_ids) == length(initial_records),
  all(packed$manifest$data$status == "complete"),
  identical(
    query_lna(
      open_lna_dataset(archive_root),
      scope = "derivative", pipeline = "fmriprep",
      space = "MNI152NLin2009cAsym", desc = "preproc"
    )$record_id,
    derivative_id
  )
)
initial_transfer <- transfer_lna_dataset(
  archive_root, receiver_root, direction = "download", chunk_bytes = 4096L
)
stopifnot(initial_transfer$verification$passed)

added <- write_source_record(
  "raw-motor",
  "sub-03/func/sub-03_task-motor_run-01_bold.nii.gz",
  list(sub = "03", task = "motor", run = "01"),
  base_values + 300,
  with_events = TRUE
)
reviewed_update <- plan_lna_dataset_update(
  archive_root,
  c(initial_records, list(added$record)),
  source_root,
  mode = "append"
)
updated <- update_lna_dataset(
  reviewed_update,
  owner_id = "installed-bids-court-update"
)
stopifnot(
  identical(updated$receipt$update_plan_hash, reviewed_update$update_plan_hash),
  dir.exists(updated$previous_generation)
)

expected <- c(initial, list(added))
names(expected) <- vapply(expected, function(item) item$record$record_id, character(1))
unlink(source_root, recursive = TRUE, force = TRUE)
stopifnot(!dir.exists(source_root))

transferred <- transfer_lna_dataset(
  archive_root, receiver_root, direction = "download", chunk_bytes = 4096L
)
verified <- verify_lna_transfer(
  receiver_root,
  expected_inventory_hash = transferred$inventory$inventory_hash
)
stopifnot(verified$passed, verified$authenticated)
received <- open_lna_dataset(receiver_root)
final_paths <- stats::setNames(
  received$manifest$data$lna_path, received$manifest$data$record_id
)
stopifnot(
  all(vapply(
    transferred$outcomes[unname(final_paths[names(expected) != "raw-motor"])],
    identical, logical(1), "reused"
  )),
  identical(transferred$outcomes[[final_paths[["raw-motor"]]]], "copied")
)
stopifnot(
  nrow(lna_dataset_view(received, "records")) == length(expected),
  identical(
    query_lna(received, sub = "03", task = "motor", run = "01")$record_id,
    "raw-motor"
  )
)

restoration <- restore_bids(received, restored_root, closure = "closed")
stopifnot(
  identical(restoration$status, "closed"),
  !length(restoration$missing_dependencies)
)
manifest <- lna_dataset_view(received, "all")
numerical_expected <- list()
for (record_id in names(expected)) {
  row <- manifest[manifest$record_id == record_id, , drop = FALSE]
  restored_relative <- sub("[.]lna[.]h5$", ".nii.gz", row$lna_path[[1L]])
  restored_path <- file.path(restored_root, restored_relative)
  stopifnot(file.exists(restored_path))
  numerical_expected[[record_id]] <- list(
    restored_relative = restored_relative,
    values = expected[[record_id]]$values
  )
  event_relative <- expected[[record_id]]$event_relative
  if (!is.null(event_relative)) {
    stopifnot(identical(
      lna_file_checksum(file.path(restored_root, event_relative)),
      expected[[record_id]]$event_checksum
    ))
  }
}
numerical_expected_path <- file.path(court_root, "numerical-expected.rds")
numerical_receipt_path <- file.path(court_root, "numerical-court.json")
numerical_log <- file.path(court_root, "numerical-court.log")
saveRDS(numerical_expected, numerical_expected_path, version = 3L)
numerical_script <- system.file(
  "consumer", "bids-exchange-numerical-court.R", package = "neuroarchive"
)
if (!nzchar(numerical_script)) stop("installed numerical court is missing", call. = FALSE)
numerical_status <- suppressWarnings(system2(
  file.path(R.home("bin"), "Rscript"),
  c(
    "--vanilla", shQuote(numerical_script), shQuote(numerical_expected_path),
    shQuote(restored_root), shQuote(numerical_receipt_path)
  ),
  stdout = numerical_log,
  stderr = numerical_log
))
if (!identical(as.integer(numerical_status), 0L)) {
  stop(
    "fresh-process numerical court failed: ",
    paste(readLines(numerical_log, warn = FALSE), collapse = "\n"),
    call. = FALSE
  )
}
numerical <- jsonlite::read_json(
  numerical_receipt_path, simplifyVector = FALSE
)
stopifnot(identical(numerical$status, "pass"))
message("BIDS exchange court: numerical subprocess PASS")

package_versions <- vapply(c(
  "neuroarchive", "hdf5r", "neuroim2", "jsonlite", "digest", "Rcpp",
  "bidser", "fmrilatent"
), function(package) {
  package_path <- find.package(package, quiet = TRUE)
  if (!nzchar(package_path)) {
    "not-installed"
  } else {
    as.character(utils::packageVersion(package, lib.loc = dirname(package_path)))
  }
}, character(1))
message("BIDS exchange court: dependency inventory PASS")

receipt <- list(
  schema_version = "1.0.0",
  status = "pass",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
  exact_candidate = list(
    tarball_sha256 = tarball_hash,
    package_version = as.character(utils::packageVersion("neuroarchive")),
    consumer_script_sha256 = lna_file_checksum(system.file(
      "consumer", "bids-exchange-installed-court.R", package = "neuroarchive"
    )),
    package_versions = as.list(package_versions),
    r_version = R.version.string,
    platform = R.version$platform
  ),
  workflows = list(
    pack = "pass",
    exact_query = "pass",
    append_update = "pass",
    source_deleted_before_transfer = TRUE,
    transfer = "pass",
    changed_only_retransfer = "pass",
    authenticated_verify = "pass",
    source_independent_closed_restore = "pass",
    lossless_numerical_restore = "pass",
    association_byte_identity = "pass"
  ),
  archive = list(
    plan_hash = received$declaration$plan_hash,
    manifest_schema_version = received$manifest$schema$manifest_schema_version,
    manifest_tsv_sha256 = lna_file_checksum(
      file.path(receiver_root, "lna_manifest.tsv")
    ),
    manifest_json_sha256 = lna_file_checksum(
      file.path(receiver_root, "lna_manifest.json")
    ),
    declaration_sha256 = lna_file_checksum(
      file.path(receiver_root, "lna_dataset.json")
    ),
    inventory_hash = transferred$inventory$inventory_hash,
    verified_file_count = verified$n_files,
    record_ids = as.list(sort(manifest$record_id))
  ),
  update = read_lna_dataset_update(archive_root),
  restoration = read_lna_restoration(restored_root),
  numerical = numerical
)
message("BIDS exchange court: durable receipt assembled")
jsonlite::write_json(
  receipt,
  file.path(court_root, "installed-court.json"),
  auto_unbox = TRUE,
  null = "null",
  digits = NA,
  pretty = TRUE
)
message("BIDS exchange installed-candidate court: PASS")
