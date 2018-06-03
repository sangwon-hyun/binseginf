context("Test contrast vector")

## contrast_vector is correct

test_that("contrast_vector works", {
  set.seed(10)
  y <- c(rep(0, 10), rep(5, 5), rep(6, 5)) + 0.01*rnorm(20)
  obj <- binSeg_fixedSteps(y, 2)
  
  res <- contrast_vector(obj, 1)
  expect_true(length(res) == 20)
  expect_true(sum(res) == 0)
  expect_true(all(res[1:10] == -1/10))
  expect_true(all(res[11:15] == 1/5))
})

test_that("contrast_vector is signed positive", {
  set.seed(10)
  y <- c(rep(0, 10), rep(-5, 5), rep(6, 5)) + 0.01*rnorm(20)
  obj <- binSeg_fixedSteps(y, 2)
  
  res <- contrast_vector(obj, 2, sorted = T)
  expect_true(res %*% y > 0)
})

test_that("contrast_vector is signed negative", {
  set.seed(10)
  y <- c(rep(0, 10), rep(-5, 5), rep(6, 5)) + 0.01*rnorm(20)
  obj <- binSeg_fixedSteps(y, 2)
  
  res <- contrast_vector(obj, 1, sorted = T)
  expect_true(res %*% y > 0)
})

test_that("contrast_vector agrees with jump direction when there's 1 jump", {
  len <- 250

  for(i in 1:len){
    set.seed(i*10)
    y <- rnorm(20)
    obj <- binSeg_fixedSteps(y, 1)
    sign1 <- sign(summary(obj)[,"Cusum"])

    poly <- polyhedra(obj)
    contrast <- contrast_vector(obj, 1)
    
    expect_true(contrast %*% y > 0)
  }
})

################################

## .contrast_vector_segment is correct

test_that(".contrast_vector_segment works", {
  set.seed(10)
  y <- c(rep(0, 10), rep(5, 5), rep(6, 5)) + 0.01*rnorm(20)
  obj <- binSeg_fixedSteps(y, 2)
  
  res <- .contrast_vector_segment(obj, 10, 20)
  
  expect_true(all(res[1:10] == -1/10))
  expect_true(all(res[11:15] == 1/5))
})

#############################################

## .contrast_vector_spike is correct

test_that(".contrast_vector_spike works", {
  set.seed(10)
  y <- c(rep(0, 10), rep(5, 5), rep(6, 5)) + 0.01*rnorm(20)
  obj <- binSeg_fixedSteps(y, 2)
  
  res <- .contrast_vector_spike(obj, 10, 20)
  
  expect_true(all(res[c(10, 11)] == c(-1, 1)))
  expect_true(all(res[1:9] == 0))
  expect_true(all(res[12:20] == 0))
})

###############################################

## .get_length.bsFs is correct

test_that(".get_length.bsFs is correct", {
  set.seed(10)
  y <- c(rep(0, 10), rep(5, 5), rep(6, 5)) + 0.01*rnorm(20)
  obj <- binSeg_fixedSteps(y, 2)
  
  res <- .get_length(obj)
  expect_true(res == 20)
})

##########################################

## .get_length.flFs is correct

test_that(".get_length.bsFs is correct", {
  set.seed(10)
  y <- c(rep(0, 10), rep(5, 5), rep(6, 5)) + 0.01*rnorm(20)
  obj <- fLasso_fixedSteps(y, 2)
  
  res <- .get_length(obj)
  expect_true(res == 20)
})
