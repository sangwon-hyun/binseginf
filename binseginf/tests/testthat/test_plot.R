context("test plot")

## test .splitChangepoints

test_that("it makes a list with x and y components",{
  res <- .splitChangepoints(100, c(1:2), 50)
  
  expect_true(length(res) == 2)
  
  nam <- sapply(res, names)
  expect_true(all(as.vector(nam) == c("x","y","x","y")))
  
  expect_true(all(res[[1]]$x == 1:51))
  expect_true(all(res[[2]]$x == 51:100))
  
  expect_true(all(res[[1]]$y == rep(1,51)))
  expect_true(all(res[[2]]$y == rep(2,50)))
})

test_that("it works when one of the segments is a singleton",{
  res <- .splitChangepoints(10, 1:3, c(5,6))
  
  expect_true(length(res) == 3)
  
  expect_true(all(res[[1]]$x == 1:6))
  expect_true(all(res[[2]]$x == c(6,7)))
  expect_true(all(res[[2]]$y == c(2,2)))
})