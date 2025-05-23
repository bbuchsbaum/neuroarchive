library(neuroarchive)

# define aggregator functions to match tests
aggregator_forward <- function(type, desc, handle) {
  input_key <- desc$inputs[[1]]
  lst <- handle$stash[[input_key]]
  mats <- lapply(lst, function(x) if (is.matrix(x)) x else as.matrix(x))
  aggregated <- do.call(rbind, mats)
  desc$version <- "1.0"
  handle$plan$add_descriptor(handle$plan$get_next_filename(type), desc)
  output_key <- desc$outputs[[1]]
  new_values <- setNames(list(aggregated), output_key)
  handle <- handle$update_stash(keys = input_key, new_values = new_values)
}

aggregator_invert <- function(type, desc, handle) {
  if (!handle$has_key("aggregated_matrix")) return(handle)
  X <- handle$get_inputs("aggregated_matrix")[[1]]
  handle <- handle$update_stash("aggregated_matrix", list(input = X))
}

assign("forward_step.myorg.aggregate_runs", aggregator_forward, envir = .GlobalEnv)
assign("invert_step.myorg.aggregate_runs", aggregator_invert, envir = .GlobalEnv)

set.seed(1)
run1_data <- matrix(rnorm(50), nrow = 10, ncol = 5)
dim(run1_data) <- c(dim(run1_data), 1)
run2_data <- matrix(rnorm(50), nrow = 10, ncol = 5)
dim(run2_data) <- c(dim(run2_data), 1)

tmp <- tempfile(fileext = ".h5")
write_lna(list(`run-01` = run1_data, `run-02` = run2_data), file = tmp,
          transforms = c("myorg.aggregate_runs", "myorg.sparsepca"),
          transform_params = list(myorg.sparsepca = list(n_components = 3)))

cat("Written to:", tmp, "\n")

h5 <- neuroarchive:::open_h5(tmp, "r")
tf_group <- h5[["transforms"]]
obj_names <- tf_group$ls()
print(obj_names)
desc_raw <- tf_group[["01_myorg.sparsepca.json"]]$read()
cat(rawToChar(desc_raw), "\n")

neuroarchive:::close_h5_safely(h5) 