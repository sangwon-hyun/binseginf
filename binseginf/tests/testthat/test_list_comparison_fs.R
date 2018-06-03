context("Test list comparison for fixed steps")

test_that(".list_comparison works", {
  set.seed(10)
  y <- c(rep(0, 20), rep(5, 10), rep(10, 5), rep(11, 5)) + 0.01*rnorm(40)
  obj <- binSeg_fixedSteps(y, 3)
  
  res <- .list_comparison(obj)
  
  expect_true(length(res) == 3)
  expect_true(is.list(res))
  
  mat <- matrix(c(1, 20, 40, 21, 30, 40, 31, 35, 40), nrow = 3, ncol = 3, byrow = T)
  for(i in 1:3){
    expect_true(all(res[[i]]$winning == mat[i,]))
    expect_true(length(which(res[[i]]$losing[,1] == mat[i,1])) == 
        mat[i,3] - mat[i,1] - 1)
  }
})

#############

## .list_comparison.cbsFs is correct

test_that(".list_comparison.cbsFs works", {
  set.seed(10)
  y <- c(rnorm(5), rnorm(5, mean = 10), rnorm(5), rnorm(5, mean = 10), rnorm(5))
  obj <- circularBinSeg_fixedSteps(y, 2)
  res <- .list_comparison.cbsFs(obj)
  
  expect_true(is.list(res))
  expect_true(all(c(sapply(res[[1]], ncol), sapply(res[[2]], ncol)) == 4))
})