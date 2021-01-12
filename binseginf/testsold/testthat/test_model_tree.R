context("Test model tree")

## .create_node is correct

test_that(".create_node makes a valid node", {
  res <- .create_node(2, 5)
  expect_true(is_valid(res))
})

test_that(".create_node errors if start is larger than end", {
  res <- expect_error(.create_node(5,2))
})

test_that(".create_node creates the right name", {
  res <- .create_node(2, 5)
  expect_true(res$name == "2-5")
})

##########################

## .get_leaves_names is correct

test_that(".get_leaves_names works", {
  data(acme, package = "data.tree")
  res <- .get_leaves_names(acme)
  expect_true(all(res == c("Go agile", "New Accounting Standards", "New Labs",
    "New Product Line", "New Software", "Outsource", "Switch to R")))
})

###############################

## .find_leadingBreakpoint is correct

test_that(".find_leadingBreakpoint will work if there is only one leaf", {
  tree <- .create_node(1, 10)
  tree$cusum <- 10
  
  res <- .find_leadingBreakpoint(tree)
  expect_true(res == "1-10")
})

#####################################

## .split_node is correct

test_that(".split_node correctly returns left and right", {
  tree <- .create_node(1, 10, 5)
  res <- .split_node(tree)
  
  expect_true(length(res) == 2)
  expect_true(res$left$start == 1)
  expect_true(res$left$end == 5)
  expect_true(res$right$start == 6)
  expect_true(res$right$end == 10)
})

#######################################

## .enumerate_splits is correct

test_that(".enumerate_splits is correct", {
  set.seed(10)
  y <- c(rep(0, 10), rep(5, 5), rep(6, 5)) + 0.01*rnorm(20)
  res <- binSeg_fixedSteps(y, 2)
  
  expect_true(all(.enumerate_splits(res$tree) == c("1-20", "11-20")))
})

test_that("the first split is always the full vector", {
  set.seed(10)
  y <- rnorm(100)
  obj <- binSeg_fixedSteps(y, 10)
  res <- .enumerate_splits(obj$tree)
  
  expect_true(res[1] == "1-100")
})

###################################

## .is_valid.bsFs is correct

test_that(".is_valid matches the tree and numSteps", {
  set.seed(10)
  y <- c(rep(0, 10), rep(5, 5), rep(6, 5)) + 0.01*rnorm(20)
  res <- binSeg_fixedSteps(y, 2)
  
  res$numSteps <- 3
  
  expect_error(is_valid(res))
})

########################################

## jumps.tree is correct

test_that("jumps.tree is correct", {
  set.seed(10)
  y <- c(rep(6, 5), rep(5, 5), rep(0, 10)) + 0.01*rnorm(20)
  obj <- binSeg_fixedSteps(y, 2)
  
  res <- jumps(obj$tree)
  
  expect_true(all(res == c(10, 5)))
})

test_that("jumps.tree can sort", {
  set.seed(10)
  y <- c(rep(6, 5), rep(5, 5), rep(0, 10)) + 0.01*rnorm(20)
  obj <- binSeg_fixedSteps(y, 2)
  
  res <- jumps(obj$tree, sorted = T)
  
  expect_true(all(res == c(5, 10)))
})