context("Test changepoint")

## test .computeJumpIdx

test_that("it splits the jump indices evenly", {
  n <- 100
  jump.loc <- seq(0,1,length.out = 5)[-c(1,5)]
  res <- .computeJumpIdx(n, jump.loc)
  d <- diff(diff(res))
  
  expect_true(all(d <= 1))
  expect_true(all(d >= -1))
})

test_that("the number of returned idx are the same number of locs",{
  set.seed(10)
  res <- .computeJumpIdx(100, sample(1:100,10)/100)
  expect_true(length(res) == 10)
})

test_that("it errors when n is too small", {
  expect_error(.computeJumpIdx(10, c(0.1,0.11)))
})

test_that("it returns numeric(0) if jump.loc is NA", {
  expect_true(length(.computeJumpIdx(100, NA)) == 0)
})

##################

## test .formMeanVec

test_that("it forms the mean vector properly", {
  jump.idx <- c(25, 50, 75)
  jump.height <- 1:4
  res <- .formMeanVec(100, jump.height, jump.idx)
  
  expect_true(length(res) == 100)
  
  d <- diff(table(res))
  expect_true(all(d <= 1))
  expect_true(all(d >= -1))
})

test_that("it works when jump.idx is numeric(0)", {
  res <- .formMeanVec(100, 1, numeric(0))
  expect_true(all(res == 1))
  expect_true(length(res) == 100)
})

######################

## test CpVector

test_that("it forms a proper CpVector class", {
  res <- CpVector(100, 1:4, c(.25,.5,.75))
  expect_true(length(res$data) == 100)
  expect_true(class(res) == "CpVector")
  
  expect_true(all(res$jump.height == 1:4))
  expect_true(all(res$jump.idx == c(25,50,75)))
})

test_that("it errors jump.height isn't one longer than jump.loc",{
  expect_error(CpVector(10, 1:4, .5))
})

test_that("it errors when n, jump.loc or jump.height are not numeric", {
  expect_error(CpVector("a", 1:2, .5))
  expect_error(CpVector(10, c("a","b"), .5))
  expect_error(CpVector(10, 1:2, "a"))
  
  res <- CpVector(10, 1:2, .5)
  expect_true(class(res) == "CpVector")
})

test_that("jump.loc must be between 0 and 1, inclusive and exclusive", {
  expect_error(CpVector(10, 1:2, -1))
  expect_error(CpVector(10, 1:2, 1))
})

test_that("when jump.loc is 0, there's a jump.idx of 1", {
  res <- CpVector(100, 1:2, 0)
  expect_true(res$jump.idx == 1)
})

test_that("it errors when jump.loc is not increasing",{
  expect_error(CpVector(100, 1:4, c(.25,.75,.5)))
  expect_error(CpVector(100, 1:4, c(.25, .6, .6)))
})

test_that("it errors when n is not a positive integer", {
  expect_error(CpVector(0, 1:2, .5))
  expect_error(CpVector(-1, 1:2, .5))
  expect_error(CpVector(50.2, 1:2,.5))
})

test_that("it errors if jump.loc has mis-NA's", {
  expect_error(CpVector(100, 1:3, c(.5, NA)))
  expect_error(CpVector(100, 1:3, c(.5, .7, NA)))
})

test_that("it fails if jump.height is not length 1 when jump.idx is NA",{
  expect_error(CpVector(100, numeric(0), NA))
  expect_error(CpVector(100, NA, NA))
  expect_error(CpVector(100, c(0,1), NA))
})

test_that("it works with no jumps", {
  set.seed(10)
  res <- CpVector(10000, 1, NA)
  expect_true(abs(mean(res$data) - 1) < 0.01)
  expect_true(res$jump.height == 1)
  expect_true(length(res$jump.idx) == 0)
})