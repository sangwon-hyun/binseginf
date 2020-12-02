context("Test test statistics")

## next_jump.bsFs is correct

test_that("next_jump.bsFs works", {
  set.seed(10)
  y <- c(rep(0, 10), rep(1, 10)) + 0.01*rnorm(20)
  fit <- binseginf::bsfs(y, 1)

  res <- next_jump(fit, y)

  expect_true(class(res) == class(fit))
})

test_that("next_jump.bsFs gives a model with one more jump", {
  set.seed(10)
  y <- c(rep(0, 10), rep(1, 10)) + 0.01*rnorm(20)
  fit <- binseginf::bsfs(y, 1)

  res <- next_jump(fit, y)

  jump_before <- binseginf::jumps(fit)
  jump_after <- binseginf::jumps(res)

  expect_true(sum(!jump_after %in% jump_before) == 1)
})

######################

## .jump_contrast is correct

test_that(".jump_contrast works", {
  set.seed(10)
  y <- c(rep(0, 10), rep(1, 10)) + 0.01*rnorm(20)
  fit1 <- binseginf::bsfs(y, 1)
  fit2 <- next_jump(fit1, y)

  res <- .jump_contrast(y, binseginf::jumps(fit1), binseginf::jumps(fit2))

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
})

test_that(".jump_contrast computes the correct value", {
  set.seed(10)
  y <- c(rep(0, 5), rep(1, 5), rep(10, 10))
  fit1 <- binseginf::bsfs(y, 1)
  fit2 <- next_jump(fit1, y)

  res <- .jump_contrast(y, binseginf::jumps(fit1), binseginf::jumps(fit2))

  expect_true(res == 1)
})

#################

## next_jump_statistic is correct

test_that("next_jump_statistic works", {
  set.seed(10)
  y <- c(rep(0, 10), rep(1, 10)) + 0.01*rnorm(20)
  fit <- binseginf::bsfs(y, 1)

  res <- next_jump_statistic(y, fit)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
})

########################

## segment_difference is correct

test_that("segment_difference works", {
  set.seed(10)
  y <- c(rep(0, 10), rep(1, 10)) + rnorm(20)
  fit <- binseginf::bsfs(y, 1)

  res <- segment_difference(y, fit, 1)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
})

test_that("segment_difference calculates correctly", {
  set.seed(5)
  y <- c(rep(0, 10), rep(1, 5), rep(0, 5)) + rnorm(20)
  fit <- binseginf::bsfs(y, 2)
  jumps <- binseginf::jumps(fit)

  res1 <- segment_difference(y, fit, 1)
  val1 <- mean(y[(jumps[1]+1):jumps[2]]) - mean(y[1:jumps[1]])

  res2 <- segment_difference(y, fit, 2)
  val2 <- mean(y[(jumps[2]+1):20]) - mean(y[(jumps[1]+1):jumps[2]])

  expect_true(sum(abs(c(res1, res2) - c(val1, val2))) < 1e-6)
})
