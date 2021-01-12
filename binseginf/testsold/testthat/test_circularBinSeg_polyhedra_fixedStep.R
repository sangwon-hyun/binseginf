context("Test circular binSeg fixed Step polyhedra")

## polyhedra.cbsFs is correct

test_that("polyhedra.cbsFs works", {
  set.seed(10)
  y <- c(rnorm(5), rnorm(5, mean = 10), rnorm(5), rnorm(5, mean = 10), rnorm(5))
  obj <- circularBinSeg_fixedSteps(y, 2)
  res <- polyhedra(obj)
  
  expect_true(is.list(res))
  expect_true(length(res) == 2)
  expect_true(all(sort(names(res)) == c("gamma", "u")))
  expect_true(ncol(res$gamma) == length(y))
  expect_true(length(res$u) == nrow(res$gamma))
  expect_true(all(res$u == 0))
})

test_that("polyhedra.cbsFs gives the right number of rows", {
  set.seed(10)
  y <- c(rnorm(5), rnorm(5, mean = 10), rnorm(5), rnorm(5, mean = 10), rnorm(5))
  obj <- circularBinSeg_fixedSteps(y, 2)
  res <- polyhedra(obj)
  
  n <- length(y)
  jump_vec <- jumps(obj, sorted = F)
  start <- jump_vec[1]; end <- jump_vec[2]
  
  ans <- 2*((nrow(.enumerate_breakpoints_cbs(n))-1) + 
    (nrow(.enumerate_breakpoints_cbs(start)) + nrow(.enumerate_breakpoints_cbs(n-end)) +
    nrow(.enumerate_breakpoints_cbs(end - start)) - 1))
  
  expect_true(nrow(res$gamma) == ans)
})

test_that("polyhedra.cbsFs satisfies polyhedra requirement", {
  set.seed(10)
  y <- c(rnorm(5), rnorm(5, mean = 10), rnorm(5), rnorm(5, mean = 10), rnorm(5))
  obj <- circularBinSeg_fixedSteps(y, 2)
  res <- polyhedra(obj)
  
  expect_true(all(res$gamma %*% y >= res$u))
})


test_that("polyhedra.cbsFs having the same model if and only if the inequalities are satisfied", {
  set.seed(5)
  y <- c(rep(0,10), rep(-2,10), rep(-1,5)) + rnorm(25)
  obj <- circularBinSeg_fixedSteps(y,2)
  
  model <- jumps(obj, sorted = F)
  poly <- polyhedra(obj)
  
  expect_true(all(poly$gamma %*% y >= poly$u))
  
  trials <- 100
  for(i in 1:trials){
    set.seed(i*10)
    y_tmp <- c(rep(0,10), rep(-2,10), rep(-1,5)) + rnorm(25)
    obj_tmp <- circularBinSeg_fixedSteps(y_tmp,2)
    
    model_tmp <- jumps(obj_tmp, sorted = F)
    
    bool1 <- all(length(model) == length(model_tmp) && all(model == model_tmp))
    bool2 <- all(poly$gamma %*% y_tmp >= poly$u)
    
    expect_true(bool1 == bool2)
  }
})

test_that("polyhedra.cbsFs having the same model if and only if the inequalities are satisfied, wrong model", {
  set.seed(5)
  y <- rnorm(25)
  obj <- circularBinSeg_fixedSteps(y,2)
  
  model <- obj$model[,1:2]
  poly <- polyhedra(obj)
  
  expect_true(all(poly$gamma %*% y >= poly$u))
  
  trials <- 100
  for(i in 1:trials){
    set.seed(i*10)
    y_tmp <- c(rnorm(5), rnorm(5, mean = -10), rnorm(5), rnorm(5, mean = 10), rnorm(5))
    obj_tmp <- circularBinSeg_fixedSteps(y_tmp,2)
    
    model_tmp <- jumps(obj_tmp, sorted = F)
    
    bool1 <- all(length(model) == length(model_tmp) && all(model == model_tmp))
    bool2 <- all(poly$gamma %*% y_tmp >= poly$u)
    
    expect_true(bool1 == bool2)
  }
})

test_that("polyhedra.cbsFs leads to uniform p values", {
  trials <- 100
  n <- 21
  null_vec <- rep(NA, trials); alt_vec <- rep(NA, trials)
  
  form_contrast <- function(obj, n){
    contrast <- rep(-1, n)
    jump_vec <- jumps(obj, sorted = T)
    if(length(jump_vec) == 1) {jump_vec <- c(0, jump_vec)}
    contrast[(jump_vec[1]+1):jump_vec[2]] <- 1
    contrast[contrast < 0] <- -1/sum(contrast < 0)
    contrast[contrast > 0] <- 1/sum(contrast > 0)
    
    contrast
  }
  
  for(i in 1:trials){
    set.seed(i*10)
    y <- rnorm(n)
    obj <- circularBinSeg_fixedSteps(y,1)
    
    poly <- polyhedra(obj)
    contrast <- form_contrast(obj, n)
  
    null_vec[i] <- pvalue(y, poly, contrast)
  }
  
  for(i in 1:trials){
    set.seed(i*10)
    y <- c(rnorm(n/3), rnorm(n/3, mean = 1), rnorm(n/3))
    obj <- circularBinSeg_fixedSteps(y,1)
    
    poly <- polyhedra(obj)
    contrast <- form_contrast(obj, n)
    
    alt_vec[i] <- pvalue(y, poly, contrast)
  }
  
  quant <- seq(0, 1, length.out = 11)
  expect_true(sum(abs(quantile(null_vec, prob = quant) - quant)) <= 
                sum(abs(quantile(alt_vec, prob = quant) - quant)) )
  
})

######################

## .cusum_cbs_contrast_full is correct

test_that(".cusum_cbs_contrast_full works", {
  res <- .cusum_cbs_contrast_full(1, c(5,10), 15, 15)
  
  expect_true(length(res) == 15)
})

test_that(".cusum_cbs_contrast_full gives the correct answer", {
  set.seed(10)
  y <- rnorm(15)
  con <- .cusum_cbs_contrast_full(1, c(5,10), 15, 15)
  cusum <- .cusum_cbs(c(5,10), cumsum(y))
  
  expect_true(abs(con%*%y - cusum) <= 1e-6)
})

test_that(".cusum_cbs_contrast_full gives the correct answer with shift", {
  set.seed(10)
  y <- rnorm(15)
  con <- .cusum_cbs_contrast_full(5, c(7,9), 10, 15)
  cusum <- .cusum_cbs(c(3,5), cumsum(y[5:10]))
  
  expect_true(abs(con%*%y - cusum) <= 1e-6)
})

######################

## .gammaRows_from_comparisons_cbsfs is correct

test_that(".gammaRows_from_comparisons_cbsfs works", {
  set.seed(10)
  y <- c(rep(0,10), rep(5,10), rep(0,10))
  
  vec <- c(1, c(11,20), 30)
  mat <- matrix(c(1,10,20,30, 1,11,19,30), 2, 4, byrow = T)
  res <- .gammaRows_from_comparisons_cbsfs(vec, mat, 1, 30)
  
  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(4,30)))
  expect_true(all(res %*% y >= 0))
})

test_that(".gammaRows_from_comparisons_cbsfs works when vec is NA", {
  set.seed(10)
  y <- c(rep(0,10), rep(5,10), rep(0,10))
  
  vec <- rep(NA, 4)
  mat <- matrix(c(1,10,20,30, 1,11,19,30), 2, 4, byrow = T)
  res <- .gammaRows_from_comparisons_cbsfs(vec, mat, 1, 30)
  
  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(4, 30)))
  
  expect_true(length(unique(res[1,1:9])) == 1)
  expect_true(length(unique(res[1,10:20])) == 1)
  expect_true(length(unique(res[1,21:30])) == 1)
  expect_true(length(unique(res[2,1:10])) == 1)
  expect_true(length(unique(res[2,11:19])) == 1)
  expect_true(length(unique(res[2,20:30])) == 1)
})

test_that(".gammaRows_from_comparisons_cbsfs works with add is T", {
  y <- c(rep(0,10), rep(5,10), rep(0,10))
  
  vec <- c(1, c(11,20), 30)
  mat <- matrix(c(1,10,20,30, 1,11,19,30), 2, 4, byrow = T)
  res <- .gammaRows_from_comparisons_cbsfs(vec, mat, 1, 30)
  
  res2 <- .gammaRows_from_comparisons_cbsfs(vec, mat, 1, 30, add = T)
  
  expect_true(all(res == res2[1:4,]))
  expect_true(length(unique(res2[5,1:10])) == 1)
  expect_true(length(unique(res2[5,11:20])) == 1)
  expect_true(length(unique(res2[5,21:30])) == 1)
})

test_that(".gammaRows_from_comparisons_cbsfs gives positive win_contrast when add is T", {
  y <- c(rep(0,10), rep(-5,10), rep(0,10))
  obj <- circularBinSeg_fixedThres(y, 5)
  comp_lis <- .list_comparison(obj)
  
  res <- .gammaRows_from_comparisons_cbsfs(comp_lis[[1]]$winning, 
                                           comp_lis[[1]]$losing, -1, 30, T)
  expect_true(nrow(res) == 2*nrow(comp_lis[[1]]$losing) + 1)
  expect_true(res[nrow(res),]%*%y >= 5)
})

test_that(".gammaRows_from_comparisons_cbsfs works when y is length 2", {
  y <- c(0,10)
  obj <- circularBinSeg_fixedThres(y, 5)
  comp_lis <- .list_comparison(obj)
  
  res <- .gammaRows_from_comparisons_cbsfs(comp_lis[[1]]$winning, 
                                           comp_lis[[1]]$losing, .get_signs_cbsFt(obj)[1], 2, T)
  expect_true(res %*% y >= 5)
  
  y <- c(0,-10)
  obj <- circularBinSeg_fixedThres(y, 5)
  comp_lis <- .list_comparison(obj)
  
  res <- .gammaRows_from_comparisons_cbsfs(comp_lis[[1]]$winning, 
                                           comp_lis[[1]]$losing, .get_signs_cbsFt(obj)[1], 2, T)
  expect_true(res %*% y >= 5)
})

test_that(".gammaRows_from_comparisons_cbsfs works when mat is NA", {
  res <- .gammaRows_from_comparisons_cbsfs(matrix(c(1,3,4,5), ncol = 4), NA, 1, 5, T)
  expect_true(is.matrix(res))
  expect_true(ncol(res) == 5)
})
