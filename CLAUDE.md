# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

neuroarchive is an R package implementing the Latent NeuroArchive (LNA) v2.0 format for compressing and storing neuroimaging data using HDF5. The package provides a flexible transform pipeline for data compression including quantization, basis decomposition (PCA), and specialized transforms like HRBF (Hemodynamic Response Basis Functions).

## Development Commands

### Building and Installation
```bash
# Install dependencies and build
R -e "devtools::install_deps()"
R -e "devtools::build()"
R -e "devtools::install()"

# Document (generate man pages from roxygen comments)
R -e "devtools::document()"

# Check package (comprehensive checks)
R -e "devtools::check()"
```

### Testing
```bash
# Run all tests using the provided script
./run-tests.sh

# Run specific test file
R -e "testthat::test_local('tests/testthat', filter='transform_basis')"

# Run all tests with detailed output
Rscript run_all_tests.R

# Run tests in R session
R -e "devtools::test()"
```

### Development Workflow
```bash
# Load package for interactive development
R -e "devtools::load_all()"

# Build and reload
R -e "devtools::build(); devtools::reload()"
```

## Architecture

### Core Components

1. **API Layer** (`R/api.R`)
   - `write_lna()`: Main entry point for writing LNA files
   - `read_lna()`: Main entry point for reading LNA files
   - `validate_lna()`: Validates LNA file structure and checksums

2. **Core Pipeline** 
   - `core_write.R`: Orchestrates the forward transform pipeline
   - `core_read.R`: Orchestrates the inverse transform pipeline
   - `materialise.R`: Handles actual HDF5 file writing with retry logic

3. **Data Handling**
   - `DataHandle` (R6 class in `handle.R`): Manages state during transforms
   - `Plan` (R6 class in `plan.R`): Tracks datasets and transform descriptors
   - `lna_reader` (`reader.R`): Lazy reading interface with subsetting

4. **Transform System**
   - S3 dispatch via `forward_step()` and `invert_step()` generics
   - Each transform implements these methods (e.g., `forward_step.quant()`)
   - Transform files follow pattern: `transform_*.R`

### Transform Pipeline Flow

1. **Writing**: 
   - Input data → Validate → Create Plan/Handle
   - Apply transforms sequentially via `forward_step()`
   - Each transform updates the stash and plan
   - `materialise_plan()` writes all data to HDF5

2. **Reading**:
   - Open HDF5 → Load transform descriptors
   - Apply inverse transforms via `invert_step()`
   - Support for ROI masking and temporal subsetting
   - Option for lazy reading with `lna_reader`

### Key Design Patterns

- **S3 Dispatch**: Transforms are implemented as S3 methods on type
- **Stash Pattern**: `DataHandle$stash` passes data between transforms
- **Plan Pattern**: Accumulates all datasets/descriptors before writing
- **Parameter Resolution**: Schema defaults → Package options → User params

### Adding New Transforms

1. Create JSON schema in `inst/schemas/<type>.schema.json`
2. Implement `forward_step.<type>()` and `invert_step.<type>()`
3. Optionally add `lna_default.<type>()` for defaults
4. Follow existing patterns in `transform_*.R` files

## Important Notes

- The package uses C++14 with OpenMP for performance-critical operations
- HDF5 operations have retry logic for chunk size and compression failures
- Parameter merging follows specific precedence (see LNA spec v1.4)
- Parallel writes should use temporary files + atomic rename
- Schema validation uses cached compiled validators (fork-safety considerations)
- Transform implementations should handle subsetting (ROI/time)

## Testing Strategy

- Unit tests for each transform's forward/inverse operations
- Integration tests for full pipelines
- Tests use in-memory HDF5 (driver='core') where possible
- ROI and temporal subsetting tested for all transforms
- Multi-run handling tested explicitly

## Dependencies

Core dependencies: hdf5r, jsonlite, Matrix, Rcpp, R6, neuroim2
Optional for specific transforms: sparsepca, irlba, wavelets, multitaper