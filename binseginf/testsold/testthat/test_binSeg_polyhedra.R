context("Test binSeg_polyhedra")

## polyhedra.bsFs is correct

test_that("polyhedra.bsFs works", {
  set.seed(10)
  y <- c(rep(1, 50), rep(10, 25), rep(5, 25)) + rnorm(100)
  obj <- binSeg_fixedSteps(y, 2)

  res <- polyhedra(obj)

  expect_true(class(res) == "polyhedra")
  expect_true(length(res) == 2)
  expect_true(all(res$u == 0))

  expect_true(all(res$gamma %*% y >= res$u))
})

test_that("it is invalid if a few inequalities are flipped",{
  set.seed(10)
  y <- c(rep(1, 50), rep(10, 25), rep(5, 25)) + rnorm(100)
  obj <- binSeg_fixedSteps(y, 2)

  res <- polyhedra(obj)
  gamma <- res$gamma
  idx <- sample(1:nrow(gamma), 5)
  gamma[idx,] <- -gamma[idx,]

  expect_true(any(gamma %*% y < res$u))
})

test_that("having the same model if and only if the inequalities are satisfied", {
  set.seed(5)
  y <- c(rep(0,50), rep(-2,20), rep(-1,30)) + rnorm(100)
  obj <- binSeg_fixedSteps(y,2)

  model.jumps <- jumps(obj)
  model.sign <- sign(jump_cusum(obj))
  poly <- polyhedra(obj)

  expect_true(all(poly$gamma %*% y >= poly$u))

  trials <- 100
  for(i in 1:trials){
    set.seed(i*10)
    y.tmp <- c(rep(0,50), rep(-2,20), rep(-1,30)) + rnorm(100)
    obj.tmp <- binSeg_fixedSteps(y.tmp,2)

    model.jumps.tmp <- jumps(obj.tmp)
    model.sign.tmp <- sign(jump_cusum(obj.tmp))

    bool1 <- (all(model.jumps.tmp == model.jumps) & all(model.sign == model.sign.tmp))
    bool2 <- all(poly$gamma %*% y.tmp >= poly$u)

    expect_true(bool1 == bool2)
  }
})

test_that("same model iff the inequalities are satisfied for no signal model", {
  set.seed(5)
  y <- rnorm(100)
  obj <- binSeg_fixedSteps(y,2)

  model.jumps <- jumps(obj)
  model.sign <- sign(jump_cusum(obj))
  poly <- polyhedra(obj)

  expect_true(all(poly$gamma %*% y >= poly$u))

  trials <- 100
  for(i in 1:trials){
    set.seed(i*10)
    y.tmp <- c(rep(0,50), rep(-2,20), rep(-1,30)) + rnorm(100)
    obj.tmp <- binSeg_fixedSteps(y.tmp,2)

    model.jumps.tmp <- jumps(obj.tmp)
    model.sign.tmp <- sign(jump_cusum(obj.tmp))

    bool1 <- (all(model.jumps.tmp == model.jumps) & all(model.sign == model.sign.tmp))
    bool2 <- all(poly$gamma %*% y.tmp >= poly$u)

    expect_true(bool1 == bool2)
  }
})

test_that("binSeg does not crash in one-jump one-length case but estimate 2 jumps", {
  set.seed(10)
  y <- c(-5, rep(0,9)) + 0.05*rnorm(10)

  obj <- binSeg_fixedSteps(y, 2)
  poly <- polyhedra(obj)

  expect_true(all(poly$gamma %*% y >= poly$u))
})

test_that("binSeg does not crash when there is multiple splits and first one is singleton", {
  set.seed(2)
  vec <- c(0, 0.05)
  dat <- CpVector(20, vec[c(1,2,1)], c(1/3, 2/3))
  y <- dat$data

  obj2 <- binSeg_fixedSteps(y, 2)
  poly <- polyhedra(obj2)

  expect_true(all(poly$gamma %*% y >= poly$u))
})

test_that("binSeg does not crash when there is multiple splits and last one is singleton", {
  set.seed(21)
  vec <- c(0, 0.05)
  dat <- CpVector(100, vec[c(1,2,1)], c(1/3, 2/3))
  y <- dat$data

  obj2 <- binSeg_fixedSteps(y, 2)
  poly <- polyhedra(obj2)

  expect_true(all(poly$gamma %*% y >= poly$u))
})

test_that("binSeg works when there are 2 changepoints", {
  vec <- c(0, 0.05)
  n <- 100
  method <- binSeg_fixedSteps

  set.seed(104)
  dat <- CpVector(n, vec[c(1,2,1)], c(1/3, 2/3))
  y <- dat$data
  obj <- method(y, 2)
  poly <- polyhedra(obj)

  expect_true(class(poly) == "polyhedra")
})

###############################

## .vector_matrix_signedDiff is correct

test_that(".vector_matrix_signedDiff works", {
  res <- .vector_matrix_signedDiff(c(1,2,3), matrix(1:6, ncol = 3, byrow = T),
    1, c(1,-1))

  expect_true(all(res == matrix(c(0,0,0,5,7,9), ncol = 3, nrow = 2, byrow = T)))
})

test_that(".vector_matrix_signedDiff gets signs right", {
  res <- .vector_matrix_signedDiff(rep(0,3), matrix(1:15, ncol = 3),
    0, c(1,1,1,-1,-1))

  expect_true(all(apply(res, 1, function(x){length(unique(sign(x)))}) == 1))
})

test_that(".vector_matrix_signedDiff returns a matrix of the right dim", {
  res <- .vector_matrix_signedDiff(1:10, matrix(1:60, ncol = 10),
    1, rep(c(-1,1), each = 3))

  expect_true(all(dim(res) == c(6, 10)))
})

###################################

## .gammaRows_from_comparisons is correct

test_that(".gammaRows_from_comparisons works", {
  set.seed(10)
  y <- sample(c(1:10))
  vec <- matrix(c(1,5,10), ncol = 3)
  mat <- cbind(1, c(1:9)[-5], 10)

  res <- .gammaRows_from_comparisons(vec, mat, 1, 10)

  expect_true(all(dim(res) == c(16, 10)))
})

test_that(".gammaRow_from_comparisons is fulfilled by y", {
  set.seed(10)
  y <- c(rep(0, 5), rep(10,4), -9) + 0.01*rnorm(10)
  obj <- binSeg_fixedSteps(y, 1)

  expect_true(obj$tree$breakpoint == 9)

  vec <- matrix(c(1,9,10), ncol = 3)
  mat <- cbind(1, 1:8, 10)

  res <- .gammaRows_from_comparisons(vec, mat, sign(obj$tree$cusum), 10)

  expect_true(all(res %*% y >= 0))
})