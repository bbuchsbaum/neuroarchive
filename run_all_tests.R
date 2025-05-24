#!/usr/bin/env Rscript

# Script to run all testthat tests and generate a report

run_tests <- function() {
  test_dir <- file.path("tests", "testthat")
  report_file <- "test_run_report.md"

  # Check if test directory exists
  if (!dir.exists(test_dir)) {
    message(sprintf("Test directory '%s' not found.", test_dir))
    return(invisible(NULL))
  }

  test_files <- list.files(test_dir, pattern = "^test-.*\.R$", full.names = FALSE)

  if (length(test_files) == 0) {
    message(sprintf("No test files found in '%s'.", test_dir))
    return(invisible(NULL))
  }

  results <- data.frame(
    File = character(),
    Filter = character(),
    Status = character(),
    Summary = character(),
    Output = character(), # To store full output for debugging
    stringsAsFactors = FALSE
  )

  message(sprintf("Found %d test files. Running tests...
", length(test_files)))

  for (test_file in test_files) {
    filter_name <- gsub("^test-(.*)\.R$", "\\1", test_file)
    message(sprintf("Running test: %s (filter: '%s')", test_file, filter_name))

    r_command <- sprintf("library(devtools); devtools::test(filter='%s', reporter='summary')", filter_name)
    
    # Using tryCatch to handle potential errors during test execution itself
    # stdout and stderr are captured
    cmd_output_list <- tryCatch({
      system2("Rscript", args = c("-e", r_command), stdout = TRUE, stderr = TRUE)
    }, error = function(e) {
      # This catches errors in executing system2 itself, not typical test failures
      message(sprintf("Error executing test script for %s: %s", test_file, e$message))
      return(list(stdout = "", stderr = sprintf("Execution Error: %s", e$message), status = 999)) # status 999 for execution error
    })

    # system2 output might be a character vector with an attribute "status"
    # stdout is the character vector itself if stdout=TRUE and stderr=FALSE
    # If both are TRUE, it's a character vector, possibly with an attribute "stderr" (logical vector)
    # and an attribute "status" for the exit code of Rscript.
    
    full_output <- ""
    if (is.character(cmd_output_list)) {
        full_output <- paste(cmd_output_list, collapse = "
")
    } else if (is.list(cmd_output_list) && !is.null(cmd_output_list$stdout) && !is.null(cmd_output_list$stderr)) {
        # Our custom error structure
        full_output <- paste(c(cmd_output_list$stdout, cmd_output_list$stderr), collapse = "
")
    }


    exit_status_attr <- attr(cmd_output_list, "status")
    execution_had_error <- FALSE
    if (is.list(cmd_output_list) && !is.null(cmd_output_list$status) && cmd_output_list$status == 999) {
        execution_had_error <- TRUE
    } else if (!is.null(exit_status_attr) && exit_status_attr != 0) {
        # Rscript itself exited with an error (e.g., package not found, syntax error in test script before summary)
        message(sprintf("Rscript exited with status %d for %s.", exit_status_attr, test_file))
        # We'll still try to parse the output for a summary, but mark as potential issue
    }

    summary_line <- ""
    status <- "ERROR (No Summary)" # Default status

    # Regex to find the summary line, e.g., [ FAIL 0 | WARN 0 | SKIP 0 | PASS 16 ]
    # Making it more robust to variations in spacing
    summary_match <- regmatches(full_output, regexpr("\[\s*FAIL\s*(\d+)\s*\|\s*WARN\s*(\d+)\s*\|\s*SKIP\s*(\d+)\s*\|\s*PASS\s*(\d+)\s*\]", full_output))

    if (length(summary_match) > 0) {
      summary_line <- summary_match[[1]]
      # Extract the number of failures
      fail_count_match <- regmatches(summary_line, regexpr("FAIL\s*(\d+)", summary_line))
      if (length(fail_count_match) > 0) {
          fail_count_str <- gsub("FAIL\s*", "", fail_count_match[[1]])
          fail_count <- as.integer(fail_count_str)
          if (is.na(fail_count)) {
            status <- "ERROR (Parse Fail)"
          } else if (fail_count == 0) {
            status <- "PASS"
          } else {
            status <- "FAIL"
          }
      } else {
          status <- "ERROR (No Fail Count)"
      }
    } else if (execution_had_error) {
        status <- "ERROR (Execution)"
        summary_line <- "Test script execution failed"
    } else if (!is.null(exit_status_attr) && exit_status_attr != 0) {
        status <- sprintf("ERROR (Rscript Exit %d)", exit_status_attr)
        summary_line <- "Rscript execution indicated an error."
    }


    results <- rbind(results, data.frame(
      File = test_file,
      Filter = filter_name,
      Status = status,
      Summary = summary_line,
      Output = substr(full_output, 1, 200), # Store a snippet of output
      stringsAsFactors = FALSE
    ))
    message(sprintf("  Status: %s, Summary: %s
", status, summary_line))
  }

  # Generate Markdown report
  report_content <- c(
    paste("# Test Run Report -", Sys.time()),
    "",
    "| Test File          | Filter             | Status | Summary                                            | Output Snippet (First 200 chars) |",
    "|--------------------|--------------------|--------|----------------------------------------------------|----------------------------------|",
    apply(results, 1, function(row) {
      # Truncate summary and output if they are too long for the table
      summary_display <- substr(row["Summary"], 1, 50)
      if (nchar(row["Summary"]) > 50) summary_display <- paste0(summary_display, "...")
      
      output_display <- gsub("
", " ", row["Output"]) # Replace newlines for table
      output_display <- substr(output_display, 1, 50)
      if (nchar(row["Output"]) > 50 || grepl("
", row["Output"])) output_display <- paste0(output_display, "...")

      sprintf("| %s | %s | %s | %s | %s |",
              row["File"], row["Filter"], row["Status"], summary_display, output_display)
    })
  )

  # Write report to file
  writeLines(report_content, report_file)
  message(sprintf("
Report generated: %s", report_file))

  # Print report to console
  message("
--- Test Run Report ---")
  for (line in report_content) {
    message(line)
  }
  message("--- End of Report ---")

  # In a non-interactive session (like Rscript), print results explicitly if not done by message
  if (!interactive()) {
    print(results[, c("File", "Filter", "Status", "Summary")])
  }
  
  return(invisible(results))
}

# Run the function
if (sys.nframe() == 0) {
  run_tests()
} 