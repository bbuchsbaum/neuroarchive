# Testing Guide for neuroarchive

## Quick Start

**Recommended for reliability:**
```bash
Rscript run-tests-serial.R
```

This runs all 96 test files one-at-a-time with HDF5 cleanup between each file.

## Available Test Runners

### 1. Serial Test Runner (Recommended)
**File:** `run-tests-serial.R`

**Usage:**
```bash
Rscript run-tests-serial.R
```

**What it does:**
- Runs each test file in isolation
- Forces HDF5 garbage collection after each file
- Adds 50ms pause to allow OS to release file handles
- Reports detailed pass/fail status for each file

**When to use:**
- Default for all local development
- When debugging test failures
- When HDF5 state corruption is suspected
- For reliable CI/CD pipelines

**Output example:**
```
=== Running 96 test files serially ===

[1/96] Testing test-api.R...
  ✓ PASSED

[2/96] Testing test-basis.R...
  ✓ PASSED

...

=== SUMMARY ===
Total files: 96
Passed: 96
Failed: 0

All tests passed! ✓
```

### 2. Fast Test Runner
**File:** `run-tests.sh`

**Usage:**
```bash
./run-tests.sh
```

**What it does:**
- Runs all tests in a single testthat session
- Faster execution (~30% faster than serial)
- May encounter HDF5 state corruption around test #153

**When to use:**
- Quick smoke tests during active development
- When modifying non-HDF5-related code
- Only if you're confident the changes won't trigger state issues

**Known issue:** Sequential test execution can leave HDF5 in corrupted state, causing later tests to fail with "unable to synchronously open file" errors.

### 3. Test Corruption Diagnostic
**File:** `diagnose-test-corruption.R`

**Usage:**
```bash
Rscript diagnose-test-corruption.R
```

**What it does:**
- Runs tests #140-160 (the problematic range)
- Checks HDF5 health before/after each test
- Reports which specific test corrupts HDF5 state
- Instruments with detailed diagnostic output

**When to use:**
- When investigating new test failures
- When serial runner also starts failing
- To identify specific culprit tests
- For debugging test suite issues

**Output example:**
```
[H5-HEALTH START] Can create/close H5 files: YES

[140/96] Testing test-transform_xyz.R...
[H5-HEALTH BEFORE #140] Can create/close H5 files: YES
  ✓ PASSED
[H5-HEALTH AFTER #140] Can create/close H5 files: YES

...

[153/96] Testing test-materialise.R...
[H5-HEALTH BEFORE #153] Can create/close H5 files: YES
  ✓ PASSED
[H5-HEALTH AFTER #153] Can create/close H5 files: NO  <-- CULPRIT FOUND
```

## Understanding HDF5 State Corruption

### What happens
When running tests sequentially in a single R session, HDF5 file handles can accumulate or library state can become corrupted. Around test #153, writes start failing silently, causing subsequent read operations to fail with:
```
Error: unable to synchronously open file '<file>'
Reason: Unable to open file
```

### Why it happens
- HDF5 handles not properly closed in some tests
- Operating system file descriptor limits reached
- HDF5 library internal state corruption
- Race conditions in parallel file operations

### How we fix it
The serial test runner prevents this by:
1. Running each test file in complete isolation
2. Forcing `hdf5r::h5garbage_collect()` after each file
3. Adding small delays for OS file handle release
4. Using `withr::defer()` in teardown hooks for aggressive cleanup

## Test Infrastructure

### Setup Hooks
**File:** `tests/testthat/setup-hooks.R`

Contains aggressive HDF5 cleanup that runs after each test:
```r
withr::defer({
  try({
    hdf5r::h5garbage_collect()
    Sys.sleep(0.01)
  }, silent = TRUE)
}, envir = teardown_env())
```

### Helper Fixtures
**File:** `tests/testthat/helper-fixtures.R`

Provides common test data generation functions for creating synthetic neuroimaging data.

## Best Practices

1. **Always use serial runner** for final verification before commits
2. **Use fast runner** for quick iteration during development
3. **Run diagnostic tool** if new test failures appear
4. **Check HDF5 cleanup** if adding new tests that create files
5. **Use `withr::defer()`** for cleanup in new test files
6. **Avoid direct file operations** - prefer in-memory testing where possible

## CI/CD Integration

Recommended CI configuration:
```yaml
test:
  script:
    - Rscript run-tests-serial.R
  timeout: 30 minutes
```

Do NOT use `./run-tests.sh` in CI pipelines as it may produce false negatives due to state corruption.

## Troubleshooting

### All tests fail with "unable to open file"
**Symptom:** Many tests fail with HDF5 file access errors

**Solution:** Use `run-tests-serial.R` instead of `run-tests.sh`

### Specific test fails in suite but passes alone
**Symptom:** Test passes when run individually but fails in full suite

**Diagnosis:** Run `diagnose-test-corruption.R` to find which earlier test corrupts state

**Fix:** Add cleanup to the culprit test or run it in isolation

### Serial runner still fails
**Symptom:** Even serial runner encounters HDF5 errors

**Possible causes:**
1. OS file descriptor limits too low
2. Permissions issue in temp directory
3. Corrupted HDF5 library installation
4. Stale file handles from previous runs

**Solutions:**
- `ulimit -n 4096` to increase file descriptor limit
- Clear temp directory: `rm -rf /tmp/*.h5`
- Restart R session
- Reinstall hdf5r package

## Quick Reference

| Command | Speed | Reliability | Use Case |
|---------|-------|-------------|----------|
| `Rscript run-tests-serial.R` | Slower | High | Default, CI/CD |
| `./run-tests.sh` | Faster | Medium | Quick checks |
| `Rscript diagnose-test-corruption.R` | Slowest | N/A | Debugging |
| `R -e "devtools::test()"` | Fast | Medium | Interactive |

**When in doubt, use `run-tests-serial.R`** ✓
