# Neuroarchive Test Errors

This document tracks test failures and errors encountered during `devtools::test()`, organized by error type.

## Syntax and Function Errors

### TICKET-001: Unrecognized escape sequence in regular expression
- **Location**: `R/plan.R:136:19`
- **Error Message**: `'\.' is an unrecognized escape in character string`
- **Root Cause**: Insufficient escaping in regular expression pattern
- **Resolution**: ✅ Fixed by using `fixed = TRUE` in `grepl()` calls

### TICKET-002: Unclosed function definition
- **Location**: `R/utils_hdf5.R:488:0`
- **Error Message**: `unexpected end of input`
- **Root Cause**: Missing closing brace `}` at the end of the `close_h5_safely` function
- **Resolution**: ✅ Fixed by adding the missing closing brace

### TICKET-003: R6Class finalize() method should be private
- **Location**: `R/reader.R`
- **Warning**: `R6Class lna_reader: finalize() method is public, but it should be private as of R6 2.4.0`
- **Root Cause**: Finalizer method was defined in the `public` section instead of `private`
- **Resolution**: ✅ Fixed by moving method to `private` section

### TICKET-007: Non-function error in close_h5_safely
- **Location**: `neuroarchive/R/api.R:101:3` and `test-api.R:57:7`
- **Error**: `Error in 'close_h5_safely(h5)': attempt to apply non-function`
- **Root Cause**: Function reference issue, likely related to loading problems
- **Status**: ⚠️ NEEDS INVESTIGATION

## API and Function Export Errors

### TICKET-008: Missing progressr export
- **Location**: `test-api.R:117:3` and `test-api.R:130:3`
- **Error**: `'handlers_is_empty' is not an exported object from 'namespace:progressr'`
- **Root Cause**: Using a non-exported function from progressr package
- **Solution**: Either: 
  1. Use `progressr:::handlers_is_empty()` with proper namespace access
  2. Implement our own check for empty handlers
  3. Update code to use officially exported progressr functions
- **Status**: 🔴 TO FIX

## Expected Test Errors (Likely Part of Test Design)

### TICKET-009: Missing HDF5 files in tests
- **Location**: Multiple tests including `test-aliases.R:23:7` and `test-api.R:87:7`
- **Error**: `Failed to open HDF5 file 'bar.h5'` and `Failed to open HDF5 file 'somefile.h5'`
- **Context**: These errors occur in tests, likely when mocking is attempted but not fully implemented
- **Analysis**: May be expected in test scenarios that check error handling
- **Status**: ℹ️ VERIFY TEST INTENT

### TICKET-010: Input validation in tests
- **Location**: Multiple tests including `test-api.R:8:3`, `test-api.R:16:3`
- **Error**: `Error in abort_lna("input data must have at least 3 dimensions"...)`
- **Context**: Occurs in tests that may be deliberately testing validation failures
- **Analysis**: Likely expected in test scenarios, but validation could be mocked for targeted testing
- **Status**: ℹ️ VERIFY TEST INTENT

## Testing Infrastructure Issues

### TICKET-011: Package version compatibility warnings
- **Warning**: `package 'hdf5r' was built under R version 4.3.3`
- **Warning**: `package 'withr' was built under R version 4.3.3`
- **Root Cause**: Packages built with newer R version than currently running
- **Impact**: Typically informational, but may cause subtle test behavior differences
- **Status**: ℹ️ INFORMATIONAL

### TICKET-012: Test suite maximum failures
- **Error**: `Maximum number of failures exceeded; quitting at end of file`
- **Analysis**: Too many test failures causing early termination
- **Temporary Solution**: Can increase with `testthat::set_max_fails(Inf)` for debugging
- **Status**: 🔄 DEPENDENT ON OTHER FIXES

## Next Steps

### High Priority
1. Investigate why `close_h5_safely` is causing "attempt to apply non-function" errors (TICKET-007)
2. Fix progressr export issue or provide an alternative (TICKET-008)
3. Review tests that may be expected to fail and ensure they're properly designed

### Medium Priority
1. Consider using mock objects more extensively in tests to avoid HDF5 file errors
2. Add regression tests for the fixed syntax issues

### Low Priority
1. Update testing environment to match package versions
2. Add test coverage for error handling scenarios
