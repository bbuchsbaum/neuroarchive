#!/bin/sh
# Simple test runner for the neuroarchive package
if ! command -v R >/dev/null 2>&1; then
  echo "R is not installed. Skipping R unit tests." >&2
  exit 0
fi

Rscript -e 'testthat::test_local("tests/testthat")'
