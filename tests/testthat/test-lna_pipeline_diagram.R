library(testthat)

# Basic DOT output

test_that("diagram returns dot string", {
pipe <- as_pipeline(array(1:4, dim = c(2,2)))
pipe$add_step(list(type = "quant", params = list(bits = 8)))
dot <- pipe$diagram("dot")
expect_true(is.character(dot))
expect_true(grepl("digraph", dot))
expect_true(grepl("quant", dot))
})

