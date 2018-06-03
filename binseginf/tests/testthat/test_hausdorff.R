context("Test Hausdorff")

## test hausdorff

test_that("it errors when set1 and set2 are not numerics", {
  expect_error(hausdorff(c("a","b"), 1:5))
  expect_error(hausdorff(1:5,c("a","b")))
})

test_that("hausdorff works properly", {
  set1 <- c(1,5,10)
  set2 <- set1 + 1
  
  expect_true(hausdorff(set1, set2) == 1)
})

test_that("hausdorff works on singletons", {
  set1 <- 1:5
  set2 <- 2
  
  expect_true(hausdorff(set1, set2) == 3)
})

test_that("hausdorff is symmetric", {
  set.seed(10)
  set1 <- sample(1:100, 10)
  set2 <- sample(1:100, 10)
  
  expect_true(hausdorff(set1, set2) == hausdorff(set2, set1))
})

test_that("it returns NA if either set is empty",{
  set1 <- numeric(0)
  set2 <- 1:5
  expect_true(is.na(hausdorff(set1, set2)))
  expect_true(is.na(hausdorff(set2, set1)))
})

test_that("it computes one-sided correctly", {
  set1 <- c(1,5,10)
  set2 <- c(2,3,6)
  expect_true(hausdorff(set1, set2, T) == 4)
  expect_true(hausdorff(set2, set1, T) == 2)
})

test_that("the maximum of both one-sided is the two-sided", {
  set.seed(10)
  set1 <- sample(1:100, 10)
  set2 <- sample(1:100, 10)
  
  target <- max(hausdorff(set1, set2, T), hausdorff(set2, set1, T))
  expect_true(hausdorff(set1, set2) == target)
})

test_that("it works on the identity (equal to 0)", {
  set.seed(10)
  set1 <- sample(1:100, 10)
  expect_true(hausdorff(set1, set1, T) == 0)
  expect_true(hausdorff(set1, set1, F) == 0)
})

test_that("it works on the singleton", {
  expect_true(hausdorff(1,5,T) == 4)
  expect_true(hausdorff(5,1) == 4)
  expect_true(hausdorff(1,5,F) == 4)
  expect_true(hausdorff(1,1) == 0)
})

test_that("the hausdorff distance is the maximum distance of intersection",{
  set.seed(5)
  len <- 10
  set1 <- sample(1:100, len)
  set2 <- sample(1:100, len)
  
  res <- hausdorff(set1, set2)
  vec1 <- numeric(len)
  for(i in 1:len){
    vec1[i] <- any((set1[i]-res):(set1[i]+res) %in% set2)
  }
  
  expect_true(all(vec1 == TRUE))
  
  vec2<- numeric(len)
  for(i in 1:len){
    vec2[i] <- any((set2[i]-res):(set2[i]+res) %in% set1)
  }
  
  expect_true(all(vec2 == TRUE))
})

test_that("the one-sided hausdorff is the max dist from set1 to set2", {
  set.seed(20)
  len <- 10
  set1 <- sample(1:100, len)
  set2 <- sample(1:100, len)
  
  res <- hausdorff(set1, set2, T)
  vec1 <- numeric(len)
  for(i in 1:len){
    vec1[i] <- any((set1[i]-res):(set1[i]+res) %in% set2)
  }
  
  expect_true(all(vec1 == TRUE))
})

test_that("it works when one set is singleton", {
  set1 <- c(1,2,8,11)
  set2 <- 4
  
  expect_true(hausdorff(set1, set2, F) == 7)
  expect_true(hausdorff(set2, set1, F) == 7)
  expect_true(hausdorff(set1, set2, T) == 7)
  expect_true(hausdorff(set2, set1, T) == 2)
})

#########################

## test jumps

test_that("it errors when vec is not numeric", {
  expect_error(jumps(c("a","b")))
})

test_that("it returns numeric(0) for vec is length 0", {
  expect_true(length(jumps(numeric(0))) == 0)
})

test_that("jumps works properly", {
  vec <- rep(1:4, each = 10)
  res <- jumps(vec)
  
  expect_true(all(res == c(10,20,30)))
})

test_that("it works when vec is length 1", {
  vec <- 5
  expect_true(length(jumps(vec)) == 0)
})

test_that("it works when vec is all singleton", {
  vec <- 1:10
  res <- jumps(vec)
  
  expect_true(all(res == 1:9))
})

test_that("combining jumps and .formMeanVec works", {
  jump.height <- c(1:4)
  jump.idx <- c(25,50,75)
  vec <- .formMeanVec(100, jump.height, jump.idx)
  res <- jumps(vec)
  
  expect_true(all(res == jump.idx))
})

test_that("jumps are measured at the end of the piecewise segment",{
  set.seed(10)
  vec <- rep(sample(1:10), times = sample(2:11))
  res <- jumps(vec)
  
  expect_true(length(res) == 9)
  for(i in 1:length(res)){
    expect_true(vec[res[i]] == vec[res[i]-1])
    expect_true(vec[res[i]] != vec[res[i]+1])
  }
})