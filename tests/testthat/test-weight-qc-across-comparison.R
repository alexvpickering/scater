# tests for cell-based pre-processing functions.
# library(scater); library(testthat); source("setup.R"); source("test-sum-across-cells.R")

library(Matrix)
library(DelayedArray)

##########################################################

test_that("internal .process_ids_with_groups method removes ids with less than min.cells in one group", {
  
  x <- rnorm(90, 7)
  ids <- c(rep(LETTERS[1:3], each = 20), rep(LETTERS[4:6], each = 10))
  groups <- rep(1:2, times = length(ids)/2)
  
  new.ids <- .process_ids_with_groups(x=x, ids=ids, groups=groups, subset = NULL, min.cells=10)
  expect_identical(unique(new.ids), c(LETTERS[1:3], NA))
})

test_that("internal .process_ids_with_groups method fails if less than min.cells in one group for each id", {
  
  x <- rnorm(90, 7)
  ids <- c(rep(LETTERS[1:3], each = 20), rep(LETTERS[4:6], each = 10))
  groups <- rep(1:2, times = length(ids)/2)
  
  expect_error(.process_ids_with_groups(x=x, ids=ids, groups=groups, subset = NULL, min.cells=11), '^need at least')
})


