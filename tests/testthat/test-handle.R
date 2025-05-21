library(testthat)
library(tibble) # Needed for Plan which might be used

# Load Plan class definition if not running via devtools::test()
# source("../R/plan.R")
# Load DataHandle class definition
# source("../R/handle.R")

# Mock Plan object for testing initialization
mock_plan <- Plan$new() # Assumes Plan class is available
# Mock H5File object
mock_h5 <- list( # Simple list mock
  filename = "test.lna.h5",
  close = function() { TRUE }
)
# Add H5File class for inherits check
class(mock_h5) <- c("H5File", class(mock_h5))

test_that("DataHandle initialization works correctly", {

  # Default initialization
  h_default <- DataHandle$new()
  expect_true(is.list(h_default$stash))
  expect_equal(length(h_default$stash), 0)
  expect_true(is.list(h_default$meta))
  expect_equal(length(h_default$meta), 0)
  expect_null(h_default$plan)
  expect_null(h_default$h5)
  expect_true(is.list(h_default$subset))
  expect_equal(length(h_default$subset), 0)

  # Initialization with values
  init_stash <- list(a = 1, b = "hello")
  init_meta <- list(dim = c(10, 5))
  init_subset <- list(time = 1:5)

  h_init <- DataHandle$new(
    initial_stash = init_stash,
    initial_meta = init_meta,
    plan = mock_plan,
    h5 = mock_h5,
    subset = init_subset
  )
  expect_identical(h_init$stash, init_stash)
  expect_identical(h_init$meta, init_meta)
  expect_identical(h_init$plan, mock_plan)
  expect_identical(h_init$h5, mock_h5)
  expect_identical(h_init$subset, init_subset)

  # Input validation checks
  expect_error(DataHandle$new(initial_stash = "not_a_list"))
  expect_error(DataHandle$new(initial_meta = 123))
  expect_error(DataHandle$new(subset = FALSE))
  expect_error(DataHandle$new(plan = list()), "must be a Plan R6 object or NULL")
  expect_error(DataHandle$new(h5 = list()), "must be an H5File object") # Check error message based on implementation

})

test_that("DataHandle exists works correctly", {
  h <- DataHandle$new(initial_stash = list(a = 1, b = NULL))

  expect_true(h$exists("a"))
  expect_true(h$exists("b")) # Key exists even if value is NULL
  expect_false(h$exists("c"))
  expect_false(h$exists("stash")) # Should not find fields

  # Check error on invalid key type
  expect_error(h$exists(123))
  expect_error(h$exists(c("a", "b")))
  expect_error(h$exists(list()))
})

test_that("DataHandle get_inputs works correctly", {
  h <- DataHandle$new(initial_stash = list(a = 1, b = "hello", c = TRUE))

  # Retrieve single key
  expect_equal(h$get_inputs("a"), list(a = 1))

  # Retrieve multiple keys
  expect_equal(h$get_inputs(c("b", "a")), list(b = "hello", a = 1))

  # Retrieve all keys
  expect_equal(h$get_inputs(c("c", "a", "b")), list(c = TRUE, a = 1, b = "hello"))

  # Error on missing key
  expect_error(
    h$get_inputs("d"),
    class = "lna_error_contract",
    regexp = "Required key\\(s\\) not found in stash: d"
  )

  # Error on partially missing keys
  # Using expect_error to capture the condition and check its fields
  err <- expect_error(
    h$get_inputs(c("a", "d", "e")),
    class = "lna_error_contract"
    # Can add regexp check here too if desired
    # regexp = "Required key\\(s\\) not found"
  )
  # Check the custom data attached to the condition
  expect_true(!is.null(err$missing_keys))
  expect_equal(sort(err$missing_keys), c("d", "e"))

  # Error on invalid key type
  expect_error(h$get_inputs(123))
  expect_error(h$get_inputs(list("a")))
  expect_error(h$get_inputs(character(0))) # Empty vector
})

test_that("DataHandle with method provides immutability", {
  # Initial object
  h1_stash <- list(a = 1)
  h1_meta <- list(orig = TRUE)
  h1 <- DataHandle$new(initial_stash = h1_stash, initial_meta = h1_meta)

  # Create h2 by updating meta
  new_meta <- list(orig = FALSE, added = 1)
  h2 <- h1$with(meta = new_meta)

  # 1. Check h2 is a new object
  expect_false(identical(h1, h2)) # Different objects

  # 2. Check h2 field was updated
  expect_identical(h2$meta, new_meta)

  # 3. Check other fields in h2 are unchanged copies
  expect_identical(h2$stash, h1$stash)
  expect_identical(h2$plan, h1$plan)
  expect_identical(h2$h5, h1$h5)
  expect_identical(h2$subset, h1$subset)
  # Ensure deep copy for lists: modify h2$stash, check h1$stash
  h2$stash$a <- 99
  expect_equal(h1$stash$a, 1) # h1$stash should NOT have changed

  # 4. Check h1 fields are unchanged
  expect_identical(h1$stash, h1_stash) # Should still be the original list
  expect_identical(h1$meta, h1_meta) # Should still be the original list

  # Test updating multiple fields
  h3 <- h1$with(stash = list(b = 2), subset = list(roi = TRUE))
  expect_identical(h3$stash, list(b = 2))
  expect_identical(h3$subset, list(roi = TRUE))
  expect_identical(h3$meta, h1_meta) # Meta should be from h1
  expect_identical(h1$stash, h1_stash) # h1 still unchanged

  # Test warning on unknown field
  expect_warning(
    h4 <- h1$with(unknown_field = 123, meta = list(x=1)),
    "Field 'unknown_field' not found in DataHandle"
  )
  expect_identical(h4$meta, list(x=1)) # meta should be updated
  # Check if unknown_field was added (it shouldn't be by current implementation)
  expect_null(h4$unknown_field)

})

test_that("DataHandle update_stash provides immutability", {
  # Initial object
  h1_stash <- list(a = 1, b = 2, c = 3)
  h1_meta <- list(orig = TRUE)
  h1 <- DataHandle$new(initial_stash = h1_stash, initial_meta = h1_meta)

  # --- Test case 1: Remove 'b', add 'd', update 'a' --- 
  h2 <- h1$update_stash(keys = "b", new_values = list(d = 4, a = 99))

  # 1a. Check h2 is new object
  expect_false(identical(h1, h2))

  # 1b. Check h2 stash is correct
  expected_h2_stash <- list(a = 99, c = 3, d = 4) # Order might vary, use setequal/sort
  expect_equal(sort(names(h2$stash)), sort(names(expected_h2_stash)))
  expect_equal(h2$stash[order(names(h2$stash))], expected_h2_stash[order(names(expected_h2_stash))])

  # 1c. Check h2 other fields are identical copies
  expect_identical(h2$meta, h1$meta)
  # ... check plan, h5, subset if they were initialized ...

  # 1d. Check h1 stash is unchanged!
  expect_identical(h1$stash, h1_stash)
  expect_identical(h1$meta, h1_meta)

  # --- Test case 2: Only remove keys --- 
  h3 <- h1$update_stash(keys = c("a", "c"), new_values = list())
  expect_equal(h3$stash, list(b = 2))
  expect_identical(h1$stash, h1_stash) # h1 unchanged

  # --- Test case 3: Only add/update keys --- 
  h4 <- h1$update_stash(keys = character(0), new_values = list(c = 5, e = 6))
  expected_h4_stash <- list(a = 1, b = 2, c = 5, e = 6)
  expect_equal(sort(names(h4$stash)), sort(names(expected_h4_stash)))
  expect_equal(h4$stash[order(names(h4$stash))], expected_h4_stash[order(names(expected_h4_stash))])
  expect_identical(h1$stash, h1_stash) # h1 unchanged

  # --- Test case 4: Remove non-existent keys --- 
  h5 <- h1$update_stash(keys = c("b", "x"), new_values = list(y = 1))
  expected_h5_stash <- list(a = 1, c = 3, y = 1)
  expect_equal(sort(names(h5$stash)), sort(names(expected_h5_stash)))
  expect_equal(h5$stash[order(names(h5$stash))], expected_h5_stash[order(names(expected_h5_stash))])
  expect_identical(h1$stash, h1_stash) # h1 unchanged

  # --- Test case 5: Empty keys and new_values --- 
  h6 <- h1$update_stash(keys = character(0), new_values = list())
  expect_identical(h6$stash, h1$stash) # Stash should be identical
  expect_false(identical(h1, h6)) # But object should be new (due to clone)
  expect_identical(h1$stash, h1_stash) # h1 unchanged

}) 