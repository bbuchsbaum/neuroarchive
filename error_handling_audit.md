# Error Handling Audit for neuroarchive Package

## Summary of Current Error Handling Patterns

### 1. stop() Usage (49 files)
- **Total occurrences**: ~95 instances
- **Common patterns**:
  - Direct `stop()` calls with sprintf formatting
  - `stop()` with `call. = FALSE`
  - Used in various contexts: validation, file operations, type checking

### 2. warning() Usage (29 files)
- **Total occurrences**: ~85 instances
- **Common patterns**:
  - Direct `warning()` calls with `call. = FALSE`
  - Conditional warnings based on validation flags
  - Progress/status warnings

### 3. stopifnot() Usage (16 files)
- **Total occurrences**: ~280 instances
- **Common patterns**:
  - Parameter validation at function entry
  - Type checking
  - Dimension validation
  - Often multiple checks in sequence

### 4. abort_lna() Usage (29 files)
- **Total occurrences**: ~380 instances (already well-adopted!)
- **Common subclasses**:
  - `lna_error_validation`
  - `lna_error_io`
  - `lna_error_descriptor`
  - `lna_error_contract`
  - `lna_error_missing_path`
  - `lna_error_transform_step`
  - `lna_error_not_implemented`
  - `lna_error_runtime`
  - `lna_error_closed_reader`
  - `lna_error_float16_unsupported`
  - `lna_error_dependency`
  - `lna_error_internal`
  - `lna_error_missing_data`
  - `lna_error_invalid_input`
  - `lna_error_sequence`
  - `lna_error_no_method`

### 5. warn_lna() Usage (3 files)
- **Total occurrences**: ~9 instances
- **Locations**:
  - core_write.R
  - hrbf_helpers.R
  - transform_spat_hrbf.R

## Key Findings

1. **Good adoption of abort_lna()**: The package already uses `abort_lna()` extensively (380+ instances), showing good progress toward standardized error handling.

2. **Limited warn_lna() adoption**: Only 9 instances of `warn_lna()` compared to 85 `warning()` calls suggests this is an area for improvement.

3. **stopifnot() prevalence**: With 280+ instances, `stopifnot()` is heavily used for parameter validation. These could be converted to more informative `abort_lna()` calls.

4. **Inconsistent error messages**: Many `stop()` calls use different formatting patterns and don't provide structured error information.

## Recommendations for Standardization

### Priority 1: Convert stop() calls
- Replace all `stop()` calls with appropriate `abort_lna()` calls
- Add proper error subclasses for better error handling
- Include location information for debugging

### Priority 2: Convert warning() calls
- Replace all `warning()` calls with `warn_lna()`
- Add warning subclasses where appropriate
- Maintain consistent warning formatting

### Priority 3: Convert stopifnot() calls
- Replace `stopifnot()` with more informative `abort_lna()` calls
- Group related validations
- Provide clear error messages explaining what validation failed

### Priority 4: Enhance existing abort_lna() calls
- Ensure all calls include appropriate subclasses
- Add location information where missing
- Improve error messages for clarity

## Files Requiring Most Attention

### High Priority (many non-standard calls):
1. **utils_hdf5.R**: 20+ stop() calls, 15+ stopifnot() calls
2. **plan.R**: 10+ stop() calls, 15+ stopifnot() calls
3. **pipeline.R**: Many abort_lna() already, but some warning() calls to convert
4. **transform_quant.R**: Mix of stop() and abort_lna()
5. **handle.R**: Many stopifnot() calls

### Medium Priority:
1. **experimental_api.R**: Mix of stopifnot() and abort_lna()
2. **validate.R**: Some warning() calls to convert
3. **materialise.R**: Mix of stop() and warning() calls
4. **core_write.R**: stopifnot() and some warning() calls

## Error Subclass Usage Analysis

Most commonly used subclasses:
1. `lna_error_validation` (most common)
2. `lna_error_descriptor`
3. `lna_error_io`
4. `lna_error_contract`
5. `lna_error_missing_path`

Consider adding new subclasses for:
- `lna_error_dimension_mismatch`
- `lna_error_type_mismatch`
- `lna_error_parameter_invalid`
- `lna_error_file_exists`
- `lna_error_permission`

## Implementation Strategy

1. **Phase 1**: Convert all stop() calls in critical paths (IO operations, transforms)
2. **Phase 2**: Convert warning() calls to warn_lna()
3. **Phase 3**: Replace stopifnot() with descriptive abort_lna() calls
4. **Phase 4**: Review and enhance existing abort_lna() calls for consistency