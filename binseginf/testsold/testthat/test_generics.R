context("Test generics")

## jump_sign is correct

test_that("jump_sign works for bsFs", {
  set.seed(10)
  y <- c(rep(0, 10), rep(-5, 10), rep(-10,10)) + 0.01*rnorm(30)
  obj <- binSeg_fixedSteps(y, 2)
  
  res <- jump_sign(obj)
  expect_true(all(dim(res) == c(2,2)))
  expect_true(all(sort(res[,1]) == c(10,20)))
  expect_true(all(res[,2] == -1))
  expect_true(all(colnames(res) == c("Jump", "Sign")))
})

test_that("jump_sign works for flFs", {
  set.seed(10)
  y <- c(rep(0, 10), rep(-5, 10), rep(0,10)) + 0.01*rnorm(30)
  obj <- fLasso_fixedSteps(y,2)
  
  res <- jump_sign(obj)
  expect_true(all(dim(res) == c(2,2)))
  expect_true(all(sort(res[,1]) == c(10,20)))
  expect_true(all(sort(res[,2]) == c(-1,1)))
  expect_true(all(colnames(res) == c("Jump", "Sign")))
})

test_that("jump_sign works for ss", {
  set.seed(10)
  y <- c(rep(0, 10), rep(-5, 10), rep(-10,10)) + 0.01*rnorm(30)
  obj <- sample_splitting(y, method = binSeg_fixedSteps, numSteps = 2)
  
  res <- jump_sign(obj)
  expect_true(all(dim(res) == c(2,2)))
  expect_true(all(sort(res[,1]) == c(10,20)))
  expect_true(all(res[,2] == -1))
  expect_true(all(colnames(res) == c("Jump", "Sign")))
})