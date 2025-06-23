pkgname <- "neuroarchive"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('neuroarchive')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("LNAFacade")
### * LNAFacade

flush(stderr()); flush(stdout())

### Name: LNAFacade
### Title: LNAFacade class
### Aliases: LNAFacade

### ** Examples

fac <- LNAFacade$new()
tmp <- tempfile(fileext = ".h5")
fac$write(array(1, dim = c(1, 1, 1)), tmp, transforms = character())
fac$read(tmp)



cleanEx()
nameEx("basis_space_centers")
### * basis_space_centers

flush(stderr()); flush(stdout())

### Name: basis_space_centers
### Title: Generate Poisson-disk sampled spatial basis centers
### Aliases: basis_space_centers

### ** Examples





cleanEx()
nameEx("basis_space_embed")
### * basis_space_embed

flush(stderr()); flush(stdout())

### Name: basis_space_embed
### Title: Embed data using spatial basis
### Aliases: basis_space_embed

### ** Examples





cleanEx()
nameEx("basis_space_reconstruct")
### * basis_space_reconstruct

flush(stderr()); flush(stdout())

### Name: basis_space_reconstruct
### Title: Reconstruct data from spatial basis embedding
### Aliases: basis_space_reconstruct

### ** Examples





cleanEx()
nameEx("basis_space_visualize")
### * basis_space_visualize

flush(stderr()); flush(stdout())

### Name: basis_space_visualize
### Title: Visualize spatial basis functions
### Aliases: basis_space_visualize

### ** Examples





cleanEx()
nameEx("basis_time")
### * basis_time

flush(stderr()); flush(stdout())

### Name: basis_time
### Title: Create temporal basis functions
### Aliases: basis_time

### ** Examples





cleanEx()
nameEx("delta_transform")
### * delta_transform

flush(stderr()); flush(stdout())

### Name: delta_transform
### Title: Apply delta (temporal difference) transform
### Aliases: delta_transform

### ** Examples





cleanEx()
nameEx("fmri_workflow")
### * fmri_workflow

flush(stderr()); flush(stdout())

### Name: fmri_workflow
### Title: Complete fMRI preprocessing and compression workflow
### Aliases: fmri_workflow

### ** Examples





cleanEx()
nameEx("hello")
### * hello

flush(stderr()); flush(stdout())

### Name: hello
### Title: Hello, World!
### Aliases: hello

### ** Examples

hello()



cleanEx()
nameEx("lna_get_transform_report")
### * lna_get_transform_report

flush(stderr()); flush(stdout())

### Name: lna_get_transform_report
### Title: Retrieve a transform report from an LNA file
### Aliases: lna_get_transform_report

### ** Examples

tmp <- tempfile(fileext = ".h5")
arr <- array(runif(6), dim = c(2,3))
write_lna(arr, tmp, transforms = "quant")
rep <- lna_get_transform_report(tmp, 0)



cleanEx()
nameEx("lna_reader")
### * lna_reader

flush(stderr()); flush(stdout())

### Name: lna_reader
### Title: Class for lazy LNA reading
### Aliases: lna_reader

### ** Examples

r <- read_lna("ex.h5", lazy = TRUE)
r$subset(roi_mask = array(TRUE, dim = c(4,4,4)), time_idx = 1)
data <- r$data()
r$close()



cleanEx()
nameEx("quant")
### * quant

flush(stderr()); flush(stdout())

### Name: quant
### Title: Add a quantisation step to a pipeline
### Aliases: quant
### Keywords: manip

### ** Examples


cleanEx()
nameEx("quantize_data")
### * quantize_data

flush(stderr()); flush(stdout())

### Name: quantize_data
### Title: Quantize data matrix with advanced compression options
### Aliases: quantize_data

### ** Examples





cleanEx()
nameEx("rcpp_control")
### * rcpp_control

flush(stderr()); flush(stdout())

### Name: rcpp_control
### Title: Control Rcpp acceleration usage
### Aliases: rcpp_control

### ** Examples





cleanEx()
nameEx("read_lna")
### * read_lna

flush(stderr()); flush(stdout())

### Name: read_lna
### Title: Read data from an LNA file
### Aliases: read_lna

### ** Examples

r <- read_lna("ex.h5", run_id = "run-*", lazy = TRUE,
              roi_mask = array(TRUE, dim = c(4,4,4)), time_idx = 1:2)
r$data()
r$close()



cleanEx()
nameEx("sparse_pca")
### * sparse_pca

flush(stderr()); flush(stdout())

### Name: sparse_pca
### Title: Sparse Principal Component Analysis
### Aliases: sparse_pca

### ** Examples





cleanEx()
nameEx("suggest_dpss_fmri")
### * suggest_dpss_fmri

flush(stderr()); flush(stdout())

### Name: suggest_dpss_fmri
### Title: Suggest DPSS parameters for fMRI applications
### Aliases: suggest_dpss_fmri

### ** Examples


# Resting-state fMRI with TR = 2s, 300 timepoints
params_rest <- suggest_dpss_fmri(TR = 2.0, n_time = 300, study_type = "resting")
basis <- temporal_basis("dpss", n_time = 300, 
                        n_basis = params_rest$n_basis,
                        time_bandwidth_product = params_rest$time_bandwidth_product)

# Task fMRI with faster sampling
params_task <- suggest_dpss_fmri(TR = 1.0, n_time = 400, study_type = "task")

# Custom frequency range
params_custom <- suggest_dpss_fmri(TR = 1.5, n_time = 350, 
                                   study_type = "custom", max_freq = 0.1)




cleanEx()
nameEx("temporal_basis")
### * temporal_basis

flush(stderr()); flush(stdout())

### Name: temporal_basis
### Title: Generate temporal basis matrix
### Aliases: temporal_basis

### ** Examples

# General-purpose DCT basis for compression
dct_basis <- temporal_basis("dct", n_time = 200, n_basis = 50)

# Simple polynomial detrending (linear + quadratic)
poly_basis <- temporal_basis("polynomial", n_time = 200, n_basis = 3)

# DPSS for preserving BOLD frequencies (0-0.08 Hz with TR=2s)
dpss_basis <- temporal_basis("dpss", n_time = 300, n_basis = 12, 
                             time_bandwidth_product = 2.0)




cleanEx()
nameEx("temporal_basis.dpss")
### * temporal_basis.dpss

flush(stderr()); flush(stdout())

### Name: temporal_basis.dpss
### Title: DPSS Temporal Basis for fMRI
### Aliases: temporal_basis.dpss

### ** Examples

## Not run: 
##D # Example 1: Conservative denoising for resting-state fMRI
##D # TR = 2s, 300 TRs (10 minutes), preserve 0-0.08 Hz
##D basis_rest <- temporal_basis("dpss", n_time = 300, n_basis = 12, 
##D                              time_bandwidth_product = 2.0)
##D 
##D # Example 2: Task fMRI with faster sampling  
##D # TR = 1s, 400 TRs, preserve 0-0.12 Hz for task frequencies
##D basis_task <- temporal_basis("dpss", n_time = 400, n_basis = 20,
##D                              time_bandwidth_product = 3.0)
##D 
##D # Example 3: High-resolution temporal filtering
##D # TR = 0.8s, 500 TRs, broader bandwidth for event-related designs
##D basis_event <- temporal_basis("dpss", n_time = 500, n_basis = 25,
##D                               time_bandwidth_product = 4.0)
##D 
##D # Example 4: Artifact removal while preserving BOLD
##D # Remove respiratory (~0.3 Hz) and cardiac (~1 Hz) with TR = 1s
##D # Use NW = 2 to concentrate energy below 0.08 Hz
##D basis_clean <- temporal_basis("dpss", n_time = 600, n_basis = 15,
##D                               time_bandwidth_product = 2.0)
## End(Not run)




cleanEx()
nameEx("temporal_basis.modwt")
### * temporal_basis.modwt

flush(stderr()); flush(stdout())

### Name: temporal_basis.modwt
### Title: MODWT Temporal Basis for fMRI
### Aliases: temporal_basis.modwt

### ** Examples

## Not run: 
##D # Optimal MODWT for resting-state fMRI (300 TRs, TR=2s)
##D basis_rest <- temporal_basis("modwt", n_time = 300, wavelet = "sym8")
##D 
##D # Task fMRI with controlled decomposition levels
##D basis_task <- temporal_basis("modwt", n_time = 400, wavelet = "sym8", levels = 5)
##D 
##D # High-resolution event-related fMRI
##D basis_event <- temporal_basis("modwt", n_time = 500, wavelet = "coif6", levels = 6)
##D 
##D # Comparison with standard DWT (requires power-of-2 length)
##D # MODWT works with any length:
##D basis_flexible <- temporal_basis("modwt", n_time = 347, wavelet = "sym8")
## End(Not run)




cleanEx()
nameEx("temporal_project_denoise")
### * temporal_project_denoise

flush(stderr()); flush(stdout())

### Name: temporal_project_denoise
### Title: Project data onto temporal basis with optional denoising
### Aliases: temporal_project_denoise

### ** Examples





cleanEx()
nameEx("validate_lna")
### * validate_lna

flush(stderr()); flush(stdout())

### Name: validate_lna
### Title: Validate an LNA file
### Aliases: validate_lna

### ** Examples

validate_lna("ex.h5")



cleanEx()
nameEx("write_lna")
### * write_lna

flush(stderr()); flush(stdout())

### Name: write_lna
### Title: Write data to an LNA file
### Aliases: write_lna

### ** Examples

x <- array(rnorm(64), dim = c(4,4,4,1))
write_lna(x, "ex.h5", transforms = "quant")
read_lna("ex.h5")

# Example with a DenseNeuroVec from neuroim2
# library(neuroim2)
# vec <- neuroim2::DenseNeuroVec(array(rnorm(64), dim = c(4,4,4,1)),
#                                space = neuroim2::NeuroSpace(c(4,4,4,1)))
# write_lna(vec, "vec_ex.h5")



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
