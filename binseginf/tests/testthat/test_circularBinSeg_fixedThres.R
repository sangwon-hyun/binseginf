context("Test circular binary segmentation fixed threshold")

## circularBinSeg_fixedThres is correct

test_that("circularBinSeg_fixedThres works", {
  set.seed(10)
  y <- c(rnorm(5), rnorm(5, mean = 10), rnorm(5), rnorm(5, mean = 10), rnorm(5))
  res <- circularBinSeg_fixedThres(y, thres = 10)
  
  expect_true(class(res) == "cbsFt")
  expect_true(all(sort(jumps(res) == c(5,10,15,20))))
})

#########

## jumps.cbsFt works

test_that("jumps.cbsFt works when there are no jumps", {
  set.seed(10)
  y <- rnorm(50)
  obj <- circularBinSeg_fixedThres(y, 50)
  jump <- jumps(obj)
  
  expect_true(length(jump) == 1)
  expect_true(is.na(jump))
})