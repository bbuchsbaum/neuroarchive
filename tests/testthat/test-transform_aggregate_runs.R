test_that("forward_step.aggregate_runs aggregates multiple runs correctly", {
  
  # Create test data - multiple 4D arrays
  run1 <- array(1, dim = c(2, 2, 2, 3))  # 2x2x2 spatial, 3 time points
  run2 <- array(2, dim = c(2, 2, 2, 3))
  run3 <- array(3, dim = c(2, 2, 2, 3))
  
  # Create mock handle with list of runs
  plan <- Plan$new()
  handle <- DataHandle$new(
    initial_stash = list(initial_input_list = list(
      run1 = run1,
      run2 = run2, 
      run3 = run3
    )),
    plan = plan
  )
  
  # Test sum aggregation (default)
  desc <- list(params = list())
  result <- forward_step.aggregate_runs("aggregate_runs", desc, handle)
  
  # Check the aggregated matrix is in stash
  expect_true("aggregated_matrix" %in% names(result$stash))
  agg_matrix <- result$stash$aggregated_matrix
  
  # Should be matrix with 3 rows (time) and 8 columns (2*2*2 voxels)
  expect_equal(dim(agg_matrix), c(3, 8))
  
  # Sum of 1+2+3 = 6 for each element
  expect_equal(unique(as.vector(agg_matrix)), 6)
  
  # Check the descriptor was added to the plan
  expect_equal(length(plan$descriptors), 1)
  
  # Get the stored descriptor from the plan
  # First, let's check what the structure actually is
  stored_desc <- plan$descriptors[[1]]
  
  # Check descriptor values in the plan
  expect_equal(stored_desc$outputs, "aggregated_matrix")
  expect_equal(stored_desc$params$runs_included, c("run1", "run2", "run3"))
  expect_equal(stored_desc$params$agg_op, "sum")
  expect_equal(stored_desc$params$orig_dims, "3x8")
})

test_that("forward_step.aggregate_runs handles mean aggregation", {
  
  run1 <- array(2, dim = c(2, 2, 2, 3))
  run2 <- array(4, dim = c(2, 2, 2, 3))
  
  plan <- Plan$new()
  handle <- DataHandle$new(
    initial_stash = list(initial_input_list = list(run1 = run1, run2 = run2)),
    plan = plan
  )
  
  desc <- list(params = list(agg_op = "mean"))
  result <- forward_step.aggregate_runs("aggregate_runs", desc, handle)
  
  agg_matrix <- result$stash$aggregated_matrix
  
  # Mean of 2 and 4 is 3
  expect_equal(unique(as.vector(agg_matrix)), 3)
})

test_that("forward_step.aggregate_runs handles concatenate aggregation", {
  
  run1 <- array(1, dim = c(2, 2, 2, 2))  # 2 time points
  run2 <- array(2, dim = c(2, 2, 2, 3))  # 3 time points
  
  plan <- Plan$new()
  handle <- DataHandle$new(
    initial_stash = list(initial_input_list = list(run1 = run1, run2 = run2)),
    plan = plan
  )
  
  desc <- list(params = list(agg_op = "concatenate"))
  result <- forward_step.aggregate_runs("aggregate_runs", desc, handle)
  
  agg_matrix <- result$stash$aggregated_matrix
  
  # Should have 5 rows (2+3 time points) and 8 columns
  expect_equal(dim(agg_matrix), c(5, 8))
  
  # First 2 rows should be 1s, next 3 should be 2s
  expect_equal(unique(as.vector(agg_matrix[1:2,])), 1)
  expect_equal(unique(as.vector(agg_matrix[3:5,])), 2)
})

test_that("forward_step.aggregate_runs handles single run passthrough", {
  
  single_run <- array(5, dim = c(2, 2, 2, 3))
  
  plan <- Plan$new()
  handle <- DataHandle$new(
    initial_stash = list(input = single_run),
    plan = plan
  )
  
  desc <- list(params = list())
  result <- forward_step.aggregate_runs("aggregate_runs", desc, handle)
  
  agg_matrix <- result$stash$aggregated_matrix
  
  # Should just reshape to matrix
  expect_equal(dim(agg_matrix), c(3, 8))
  expect_equal(unique(as.vector(agg_matrix)), 5)
  
  # Check the descriptor was added to the plan
  expect_equal(length(plan$descriptors), 1)
  stored_desc <- plan$descriptors[[1]]
  expect_equal(stored_desc$params$runs_included, "single_run")
  expect_equal(stored_desc$params$orig_dims, "3x8")
})

test_that("forward_step.aggregate_runs handles matrix inputs", {
  
  # Already matrix inputs
  mat1 <- matrix(1, nrow = 5, ncol = 10)
  mat2 <- matrix(2, nrow = 5, ncol = 10)
  
  plan <- Plan$new()
  handle <- DataHandle$new(
    initial_stash = list(initial_input_list = list(a = mat1, b = mat2)),
    plan = plan
  )
  
  desc <- list(params = list())
  result <- forward_step.aggregate_runs("aggregate_runs", desc, handle)
  
  agg_matrix <- result$stash$aggregated_matrix
  
  # Should sum the matrices
  expect_equal(dim(agg_matrix), c(5, 10))
  expect_equal(unique(as.vector(agg_matrix)), 3)
})

test_that("forward_step.aggregate_runs errors on invalid aggregation operation", {
  
  plan <- Plan$new()
  handle <- DataHandle$new(
    initial_stash = list(initial_input_list = list(a = matrix(1, 2, 2))),
    plan = plan
  )
  
  desc <- list(params = list(agg_op = "invalid"))
  
  expect_error(
    forward_step.aggregate_runs("aggregate_runs", desc, handle),
    "Unsupported aggregation operation"
  )
})

test_that("forward_step.aggregate_runs errors on missing input", {
  
  plan <- Plan$new()
  handle <- DataHandle$new(
    initial_stash = list(),  # Empty stash
    plan = plan
  )
  
  desc <- list(params = list())
  
  expect_error(
    forward_step.aggregate_runs("aggregate_runs", desc, handle),
    "No input data found"
  )
})