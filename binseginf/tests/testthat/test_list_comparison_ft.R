context("Test list comparison for fixed threshold")

## .list_comparison.cbsFt is correct

test_that(".list_comparison.cbsFt works", {
  set.seed(10)
  y <- c(rnorm(5), rnorm(5, mean = 50), rnorm(5), rnorm(5, mean = 10), rnorm(5))
  obj <- circularBinSeg_fixedThres(y, thres = 10)
  res <- .list_comparison.cbsFt(obj)
  
  expect_true(is.list(res))
  expect_true(all(as.numeric(sapply(res, function(x){c(ncol(x$winning), ncol(x$losing))})) == 4))
  expect_true(length(res) == c(1 + 3 + 3))
})

test_that(".list_comparison.cbsFt gives the right number of elements with NA winning", {
  set.seed(10)
  y <- c(rnorm(5), rnorm(5, mean = 50), rnorm(5), rnorm(5, mean = 10), rnorm(5))
  obj <- circularBinSeg_fixedThres(y, thres = 10)
  res <- .list_comparison.cbsFt(obj)
  
  vec <- sapply(res, function(x){ifelse(all(is.na(x$winning)), TRUE, FALSE)})
  expect_true(length(which(vec)) == 5)
})

test_that(".list_comparison.cbsFt works when there are no splits", {
  set.seed(10)
  n <- 21
  y <- rnorm(n)
  obj <- circularBinSeg_fixedThres(y,2)
  res <- .list_comparison.cbsFt(obj)
  
  expect_true(is.list(res))
  expect_true(all(as.numeric(sapply(res, function(x){c(ncol(x$winning), ncol(x$losing))})) == 4))
})

test_that(".list_comparison.cbsFt does not give incorrect losing for corner case", {
  set.seed(10)
  n <- 21
  y <- rnorm(n)
  obj <- circularBinSeg_fixedThres(y,0.5)
  res <- .list_comparison.cbsFt(obj)
  
  expect_true(all(sapply(res, function(x){ncol(x$losing)}) == 4))
  
})