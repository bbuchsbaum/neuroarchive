# Transform Registry Demo
# This demonstrates the new transform registry system in neuroarchive

library(neuroarchive)

# The registry is automatically populated on package load
# Let's see what transforms are available

# List all registered transforms
transforms <- list_transforms()
cat("Available transforms:\n")
print(transforms)

# Get detailed information about a specific transform
quant_info <- get_transform_info("quant")
cat("\n\nQuantization transform details:\n")
str(quant_info)

# Check if a transform is registered
cat("\n\nIs 'temporal' registered?", is_transform_registered("temporal"), "\n")
cat("Is 'nonexistent' registered?", is_transform_registered("nonexistent"), "\n")

# Get capabilities of a transform
basis_caps <- get_transform_capabilities("basis")
cat("\n\nBasis transform capabilities:\n")
print(basis_caps)

# List transforms by category
cat("\n\nCompression transforms:\n")
compression_transforms <- list_transforms(category = "compression")
print(compression_transforms)

cat("\n\nDimensionality reduction transforms:\n")
dim_transforms <- list_transforms(category = "dimensionality")
print(dim_transforms)

# List transforms by capability
cat("\n\nTransforms with compression capability:\n")
compression_capable <- list_transforms(capability = "compression")
print(compression_capable)

cat("\n\nTransforms with spatial capability:\n")
spatial_transforms <- list_transforms(capability = "spatial")
print(spatial_transforms)

# Get a summary of all transforms
cat("\n\nTransform summary:\n")
summary_df <- transform_summary()
print(summary_df)

# Get detailed information including all metadata
cat("\n\nDetailed transform information:\n")
detailed <- list_transforms(with_details = TRUE)
cat("Number of transforms with details:", length(detailed), "\n")

# Example: Check specific transform features
if (is_transform_registered("quant")) {
  info <- get_transform_info("quant")
  cat("\n\nQuantization transform features:\n")
  cat("- Has forward step:", info$has_forward, "\n")
  cat("- Has invert step:", info$has_invert, "\n")
  cat("- Has schema:", info$has_schema, "\n")
  cat("- Min dimensions:", info$min_dims, "\n")
  cat("- Is lossy:", info$capabilities$lossy, "\n")
  cat("- Category:", info$metadata$category, "\n")
}

# The registry can be used to discover what transforms are available
# and their capabilities before using them in a pipeline
cat("\n\nThis registry system enables:\n")
cat("- Runtime discovery of available transforms\n")
cat("- Capability-based transform selection\n")
cat("- Validation of transform requirements\n")
cat("- Better error messages when transforms are missing\n")