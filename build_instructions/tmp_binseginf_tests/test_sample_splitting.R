context("Test sample splitting")

## sample_splitting is correct

test_that("sample_splitting works with binSeg_fixedSteps", {
  set.seed(10)
  y <- c(rep(0, 50), rep(5, 50)) + 0.1*rnorm(100)
  
  res <- sample_splitting(y, binSeg_fixedSteps, numSteps = 1)
  
  expect_true(class(res) == "ss")
  expect_true(res$method == "binSeg_fixedSteps")
  expect_true(res$jumps == 50)
})

test_that("sample_splitting works with fLasso_fixedSteps", {
  set.seed(10)
  y <- c(rep(0, 50), rep(5, 50)) + 0.1*rnorm(100)
  
  res <- sample_splitting(y, fLasso_fixedSteps, numSteps = 1)
  
  expect_true(class(res) == "ss")
  expect_true(res$method == "fLasso_fixedSteps")
  expect_true(res$jumps == 50)
})

###################################

## jumps.ss is correct

test_that("jumps.ss works", {
  set.seed(10)
  y <- c(rep(0, 50), rep(5, 50)) + 0.1*rnorm(100)
  
  obj <- sample_splitting(y, binSeg_fixedSteps, numSteps = 1)
  res <- jumps(obj)
  
  expect_true(res == 50)
})

##################################

## contrast_vector_ss is correct

test_that("contrast_vector_ss works", {
  set.seed(10)
  y <- c(rep(0, 50), rep(5, 50)) + 0.1*rnorm(100)
  
  obj <- sample_splitting(y, binSeg_fixedSteps, numSteps = 1)
  
  res <- contrast_vector_ss(obj, 1)
  expect_true(all(which(res < 0) == seq(3, 49, by = 2)))
  expect_true(all(which(res > 0) == seq(53, 99, by = 2)))
})

test_that("contrast_vector_ss properly handles when n is odd", {
  set.seed(10)
  y <- c(rep(0, 50), rep(5, 49)) + 0.1*rnorm(99)
  
  obj <- sample_splitting(y, binSeg_fixedSteps, numSteps = 1)
  
  res <- contrast_vector_ss(obj, 1)
  expect_true(all(which(res < 0) == seq(3, 49, by = 2)))
  expect_true(all(which(res > 0) == seq(53, 97, by = 2)))
})

test_that("contrast_vector_ss properly when jump is at 4", {
  set.seed(10)
  y <- c(rep(0, 4), rep(100, 96)) + 0.1*rnorm(100)
  
  obj <- sample_splitting(y, binSeg_fixedSteps, numSteps = 1)
  
  res <- contrast_vector_ss(obj, 1)
  expect_true(all(which(res < 0) == 3))
  expect_true(all(which(res > 0) == seq(7, 99, by = 2)))
})

test_that("contrast_vector_ss properly when jump is at 3", {
  set.seed(10)
  y <- c(rep(0, 3), rep(100, 97)) + 0.1*rnorm(100)
  
  obj <- sample_splitting(y, binSeg_fixedSteps, numSteps = 1)
  
  res <- contrast_vector_ss(obj, 1)
  expect_true(all(which(res < 0) == 1))
  expect_true(all(which(res > 0) == seq(5, 99, by = 2)))
})

test_that("contrast_vector_ss properly when jump is at 2", {
  set.seed(10)
  y <- c(rep(0, 2), rep(100, 98)) + 0.1*rnorm(100)
  
  obj <- sample_splitting(y, binSeg_fixedSteps, numSteps = 1)
  
  res <- contrast_vector_ss(obj, 1)
  expect_true(all(which(res < 0) == 1))
  expect_true(all(which(res > 0) == seq(5, 99, by = 2)))
})

test_that("contrast_vector_ss properly when jump is at 96", {
  set.seed(10)
  y <- c(rep(0, 96), rep(100, 4)) + 0.1*rnorm(100)
  
  obj <- sample_splitting(y, binSeg_fixedSteps, numSteps = 1)
  
  res <- contrast_vector_ss(obj, 1)
  expect_true(all(which(res < 0) == seq(3, 95, by = 2)))
  expect_true(all(which(res > 0) == 99))
})

test_that("contrast_vector_ss properly when jump is at 99", {
  set.seed(10)
  y <- c(rep(0, 99), rep(100, 1)) + 0.1*rnorm(100)
  
  obj <- sample_splitting(y, binSeg_fixedSteps, numSteps = 1)
  
  res <- contrast_vector_ss(obj, 1)
  expect_true(all(which(res < 0) == seq(3, 97, by = 2)))
  expect_true(all(which(res > 0) == 99))
})

test_that("contrast_vector_ss properly when jump is at 98", {
  set.seed(10)
  y <- c(rep(0, 98), rep(100, 2)) + 0.1*rnorm(100)
  
  obj <- sample_splitting(y, binSeg_fixedSteps, numSteps = 1)
  
  res <- contrast_vector_ss(obj, 1)
  expect_true(all(which(res < 0) == seq(3, 97, by = 2)))
  expect_true(all(which(res > 0) == 99))
})

test_that("contrast_vector_ss is signed positive", {
  set.seed(10)
  y <- c(rep(0, 50), rep(5, 50)) + 0.1*rnorm(100)
  
  obj <- sample_splitting(y, binSeg_fixedSteps, numSteps = 1)
  
  res <- contrast_vector_ss(obj, 1)
  expect_true(attr(res, "sign") == 1)
})

test_that("contrast_vector_ss is signed negative", {
  set.seed(10)
  y <- c(rep(0, 50), rep(-5, 50)) + 0.1*rnorm(100)
  
  obj <- sample_splitting(y, binSeg_fixedSteps, numSteps = 1)
  
  res <- contrast_vector_ss(obj, 1)
  expect_true(attr(res, "sign") == -1)
})

#####################################

## pvalue_ss is correct

test_that("pvalue_ss works", {
  set.seed(10)
  y <- c(rep(0, 50), rep(5, 50)) + 0.1*rnorm(100)
  
  obj <- sample_splitting(y, binSeg_fixedSteps, numSteps = 1)
  v <- contrast_vector_ss(obj, 1)
  
  res <- pvalue_ss(y, v)
  expect_true(length(res) == 1)
  expect_true(res < 0.05)
})

test_that("pvalue_ss works for the reverse direction", {
  set.seed(10)
  y <- c(rep(0, 50), rep(-5, 50)) + 0.1*rnorm(100)
  
  obj <- sample_splitting(y, binSeg_fixedSteps, numSteps = 1)
  v <- contrast_vector_ss(obj, 1)
  
  res <- pvalue_ss(y, v)
  expect_true(length(res) == 1)
  expect_true(res < 0.05)
})


test_that("pvalue_ss is roughly uniform under no signal", {
  trials <- 500
  res_vec <- numeric(trials)
  
  for(trial in 1:trials){
    set.seed(trial)
    y <- rnorm(100)
    
    obj <- sample_splitting(y, binSeg_fixedSteps, numSteps = 1)
    v <- contrast_vector_ss(obj, 1)
    
    res_vec[trial] <- pvalue_ss(y, v)
  }
  
  expect_true(sum(abs(sort(res_vec) - seq(0, 1, length.out = trials))) < 15)
})

###################################

## confidence_interval_ss is correct

test_that("confidence_interval_ss works", {
  set.seed(10)
  y <- c(rep(0, 50), rep(5, 50)) + 0.1*rnorm(100)
  
  obj <- sample_splitting(y, binSeg_fixedSteps, numSteps = 1)
  v <- contrast_vector_ss(obj, 1)
  
  res <- confidence_interval_ss(y, v, 0.95)
  
  expect_true(length(res) == 2)
  expect_true(is.numeric(res))
  expect_true(res[1] <= 5 & 5 <= res[2])
})

test_that("confidence_interval_ss gets wider if alpha gets bigger", {
  set.seed(10)
  y <- c(rep(0, 50), rep(5, 50)) + 0.1*rnorm(100)
  
  obj <- sample_splitting(y, binSeg_fixedSteps, numSteps = 1)
  v <- contrast_vector_ss(obj, 1)
  
  res1 <- confidence_interval_ss(y, v, 0.95)
  res2 <- confidence_interval_ss(y, v, 0.5)
  
  expect_true(diff(res1) >= diff(res2))
})

test_that("confidence_interval_ss has the right coverage", {
  trials <- 1000
  res_vec <- numeric(trials)
  
  for(trial in 1:trials){
    set.seed(trial)
    y <- rnorm(100)
    
    obj <- sample_splitting(y, binSeg_fixedSteps, numSteps = 1)
    v <- contrast_vector_ss(obj, 1)
    
    tmp <- confidence_interval_ss(y, v, alpha = 0.95)
    res_vec[trial] <- ifelse(tmp[1] <= 0 & 0 <= tmp[2], 1, 0)
  }
  
  res <- sum(res_vec)
  
  sd_est <- sqrt(trials * 0.95 * 0.05)
  expect_true(0.95*trials - 3*sd_est <= res)
  expect_true(res <= 0.95*trials + 3*sd_est)
})

########################

## .refit_sample_split is correct

test_that(".refit_sample_split works", {
  set.seed(10)
  y <- rnorm(100)
  jump <- c(4,10,58)
  res <- .refit_sample_split(y, jump)
  
  expect_true(length(unique(res)) == 4)
  expect_true(all(res[1:4] == mean(y[c(2,4)])))
  expect_true(all(res[5:10] == mean(y[seq(6,10,by=2)])))
  expect_true(all(res[11:58] == mean(y[seq(12,58,by=2)])))
  expect_true(all(res[59:100] == mean(y[seq(60,100,by=2)])))
})

test_that(".refit_sample_split works on the boundary", {
  set.seed(10)
  y <- rnorm(100)
  jump <- c(2,10,98)
  res <- .refit_sample_split(y, jump)
  
  expect_true(length(unique(res)) == 4)
  expect_true(all(res[1:2] == y[2]))
  expect_true(all(res[3:10] == mean(y[seq(4,10,by=2)])))
  expect_true(all(res[11:98] == mean(y[seq(12,98,by=2)])))
  expect_true(all(res[99:100] == y[100]))
})