# Assemble and gate the durable exact-candidate BIDS exchange receipt.
#
# Usage:
#   Rscript finalize-bids-exchange-certification.R \
#     installed-court.json validator.json neuroarchive.tar.gz \
#     00check.log bidser.tar.gz certification.json

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L || any(!nzchar(args))) {
  stop(
    paste0(
      "expected installed court, validator JSON, neuroarchive tarball, ",
      "R CMD check log, bidser tarball, and output paths"
    ),
    call. = FALSE
  )
}
names(args) <- c(
  "court", "validator", "tarball", "check_log", "bidser_tarball", "output"
)
required_inputs <- args[setdiff(names(args), "output")]
missing <- required_inputs[!file.exists(required_inputs)]
if (length(missing)) {
  stop("missing certification input: ", paste(missing, collapse = ", "), call. = FALSE)
}

`%||%` <- function(left, right) if (is.null(left)) right else left

sha256 <- function(path) {
  paste0("sha256:", digest::digest(
    path, algo = "sha256", serialize = FALSE, file = TRUE
  ))
}
read_json <- function(path) {
  jsonlite::read_json(path, simplifyVector = FALSE)
}

court <- read_json(args[["court"]])
validator <- read_json(args[["validator"]])
tarball_hash <- sha256(args[["tarball"]])
bidser_tarball_hash <- sha256(args[["bidser_tarball"]])
if (!identical(court$status, "pass") ||
    !identical(court$exact_candidate$tarball_sha256, tarball_hash)) {
  stop("installed court does not belong to the exact candidate tarball", call. = FALSE)
}
required_workflows <- c(
  "pack", "exact_query", "append_update", "transfer",
  "changed_only_retransfer", "authenticated_verify",
  "source_independent_closed_restore", "lossless_numerical_restore",
  "association_byte_identity"
)
for (name in required_workflows) {
  if (!identical(court$workflows[[name]], "pass")) {
    stop("installed workflow did not pass: ", name, call. = FALSE)
  }
}
if (!identical(court$workflows$source_deleted_before_transfer, TRUE)) {
  stop("installed court did not delete source data before transfer", call. = FALSE)
}

errors <- validator$issues$errors %||% list()
warnings <- validator$issues$warnings %||% list()
if (length(errors)) {
  stop("BIDS Validator reported errors", call. = FALSE)
}
warning_keys <- sort(unique(vapply(warnings, `[[`, character(1), "key")))
accepted_warnings <- list(
  INCONSISTENT_SUBJECTS = paste0(
    "intentional heterogeneous task coverage in the multi-subject fixture; ",
    "the warning is advisory and every represented file is structurally valid"
  ),
  MISSING_SESSION = paste0(
    "intentional sparse-session fixture: sub-01 has ses-A while sub-02 and ",
    "sub-03 exercise legal no-session layouts"
  ),
  TOO_FEW_AUTHORS = paste0(
    "synthetic generated restoration has one declared generator author; ",
    "this is bibliographic advice, not a structural defect"
  )
)
if (!identical(warning_keys, sort(names(accepted_warnings)))) {
  stop(
    "validator warning keys differ from the reviewed disposition set: ",
    paste(warning_keys, collapse = ", "), call. = FALSE
  )
}
warning_dispositions <- lapply(warning_keys, function(key) list(
  key = key,
  disposition = "accepted-non-structural",
  rationale = accepted_warnings[[key]]
))

check_lines <- readLines(args[["check_log"]], warn = FALSE, encoding = "UTF-8")
status_line <- grep("^Status:", check_lines, value = TRUE)
if (length(status_line) != 1L ||
    grepl("[1-9][0-9]* ERROR|[1-9][0-9]* WARNING", status_line)) {
  stop("R CMD check log is not a zero-error, zero-warning result", call. = FALSE)
}

receipt <- list(
  schema_version = "1.0.0",
  status = "pass",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
  exact_candidate = list(
    neuroarchive_tarball_sha256 = tarball_hash,
    bidser_0_5_1_tarball_sha256 = bidser_tarball_hash,
    installed_consumer_sha256 = court$exact_candidate$consumer_script_sha256,
    installed_package_version = court$exact_candidate$package_version,
    package_versions = court$exact_candidate$package_versions,
    r_version = court$exact_candidate$r_version,
    platform = court$exact_candidate$platform
  ),
  gates = list(
    installed_workflows = "pass",
    source_deleted_before_transfer = TRUE,
    bids_validator_version = "1.15.0",
    bids_validator_errors = 0L,
    bids_validator_warning_count = length(warning_keys),
    r_cmd_check_status = sub("^Status:[[:space:]]*", "", status_line)
  ),
  warning_dispositions = warning_dispositions,
  archive = court$archive,
  update = court$update,
  restoration = court$restoration,
  installed_court = court,
  bids_validator = validator,
  evidence_hashes = list(
    installed_court_json_sha256 = sha256(args[["court"]]),
    bids_validator_json_sha256 = sha256(args[["validator"]]),
    r_cmd_check_log_sha256 = sha256(args[["check_log"]])
  )
)
dir.create(dirname(args[["output"]]), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(
  receipt, args[["output"]], auto_unbox = TRUE, null = "null",
  digits = NA, pretty = TRUE
)
message("Exact installed-candidate BIDS certification: PASS")
