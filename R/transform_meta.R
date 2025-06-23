#' Minimum input dimensionality for a transform
#'
#' Provides the minimum number of dimensions required by a transform's
#' forward step. Packages can define methods for their own transforms.
#' The default requirement is 3 dimensions if a specific transform is not listed.
#'
#' @param type Character transform type.
#' @export
transform_min_dims <- function(type) {
  # Ensure type is a character string
  if (!is.character(type)) {
    type <- as.character(type)
  }
  
  # Take only the first element if multiple are provided
  if (length(type) > 1) {
    type <- type[1]
  }
  
  # Handle NULL or empty cases
  if (length(type) == 0 || is.null(type)) {
    return(3L)
  }
  
  # message(paste0("[transform_min_dims] called for type: ", type))
  switch(type,
         delta = 1L,
         quant = 1L,
         basis = 2L,
         embed = 2L,
         # myorg.sparsepca would also be 2L if it had its own entry
         3L # Default value if type is not matched
  )
}

# Old S3 methods are now fully removed/commented correctly to avoid NAMESPACE issues.
# No @export tags should remain for these.

# # transform_min_dims.default <- function(type) {
# #   3L
# # }
# 
# # transform_min_dims.quant <- function(type) {
# #  1L
# # }
# # 
# # transform_min_dims.basis <- function(type) {
# #  2L
# # }
# # 
# # transform_min_dims.embed <- function(type) {
# #  2L
# # }
# # 
# # transform_min_dims.delta <- function(type) {
# #  1L
# # }
