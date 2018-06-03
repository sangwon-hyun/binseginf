context("Test circular binSeg fixed Threshold polyhedra")

## .get_signs_cbsFt is correct

test_that(".get_signs_cbsFt works", {
  set.seed(10)
  y <- c(rnorm(5), rnorm(5, mean = 10), rnorm(5), rnorm(5, mean = 10), rnorm(5))
  obj <- circularBinSeg_fixedThres(y, thres = 10)
  res <- .get_signs_cbsFt(obj)
  
  expect_true(length(res) == 7)
  expect_true(all(res[3:7] == 0))
  expect_true(all(res[1:2] != 0))
})

############

## polyhedra.cbsFt is correct

test_that("polyhedra.cbsFt works", {
  set.seed(10)
  y <- c(rnorm(5), rnorm(5, mean = 10), rnorm(5), rnorm(5, mean = 10), rnorm(5))
  obj <- circularBinSeg_fixedThres(y, thres = 10)
  res <- polyhedra(obj)
  
  expect_true(class(res) == "polyhedra")
  expect_true(all(names(res) == c("gamma", "u")))
  expect_true(ncol(res$gamma) == length(y))
  expect_true(nrow(res$gamma) == length(res$u))
})

test_that("polyhedra.cbsFt satisfies the polyhedra requirement for 1 jump", {
  set.seed(10)
  y <- c(rnorm(5), rnorm(5, mean = 10), rnorm(5))
  obj <- circularBinSeg_fixedThres(y, thres = 10)
  res <- polyhedra(obj)

  expect_true(all(res$gamma %*% y >= res$u))
})

test_that("polyhedra.cbsFt satisfies the polyhedra requirement", {
  set.seed(10)
  y <- c(rnorm(5), rnorm(5, mean = 10), rnorm(5), rnorm(5, mean = 10), rnorm(5))
  obj <- circularBinSeg_fixedThres(y, thres = 10)
  res <- polyhedra(obj)

  expect_true(all(res$gamma %*% y >= res$u))
})

test_that("polyhedra.cbsFt gives the right number of rows", {
  set.seed(10)
  n <- 10
  y <- c(rnorm(n/2), rnorm(n/2, mean = 10))
  obj <- circularBinSeg_fixedThres(y,10)
  poly <- polyhedra(obj)

  expect_true(nrow(poly$gamma) == (nrow(.enumerate_breakpoints_cbs(n))*2 + 1 - 2) +
                nrow(.enumerate_breakpoints_cbs(n/2))*4)
})


test_that("polyhedra.cbsFt having the same model if and only if the inequalities are satisfied", {
  set.seed(5)
  y <- c(rep(0,10), rep(-2,10), rep(-1,5)) + rnorm(25)
  obj <- circularBinSeg_fixedThres(y,2)
  
  model <- jumps(obj, sorted = F)
  poly <- polyhedra(obj)
  
  expect_true(all(poly$gamma %*% y >= poly$u))
  
  trials <- 100
  for(i in 1:trials){
    set.seed(i*10)
    y_tmp <- c(rep(0,10), rep(-2,10), rep(-1,5)) + rnorm(25)
    obj_tmp <- circularBinSeg_fixedThres(y_tmp,2)
    
    model_tmp <- jumps(obj_tmp, sorted = F)
    
    bool1 <- all(length(model) == length(model_tmp) && all(model == model_tmp))
    bool2 <- all(poly$gamma %*% y_tmp >= poly$u)
    
    expect_true(bool1 == bool2)
  }
})

test_that("polyhedra.cbsFt having the same model if and only if the inequalities are satisfied, wrong model", {
  set.seed(5)
  y <- rnorm(25)
  obj <- circularBinSeg_fixedThres(y,2)
  
  model <- obj$model[,1:2]
  poly <- polyhedra(obj)
  
  expect_true(all(poly$gamma %*% y >= poly$u))
  
  trials <- 100
  for(i in 1:trials){
    set.seed(i*10)
    y_tmp <- c(rnorm(5), rnorm(5, mean = -10), rnorm(5), rnorm(5, mean = 10), rnorm(5))
    obj_tmp <- circularBinSeg_fixedThres(y_tmp,2)
    
    model_tmp <- jumps(obj_tmp, sorted = F)
    
    bool1 <- all(length(model) == length(model_tmp) && all(model == model_tmp))
    bool2 <- all(poly$gamma %*% y_tmp >= poly$u)
    
    expect_true(bool1 == bool2)
  }
})

test_that("polyhedra.cbsFt works when there are too many splits", {
  set.seed(10)
  n <- 21
  y <- rnorm(n)
  obj <- circularBinSeg_fixedThres(y,0.5)
  poly <- polyhedra(obj)

  expect_true(class(poly) == "polyhedra")
  expect_true(all(dim(poly$gamma) == c(length(poly$u), length(y))))
  expect_true(all(poly$gamma %*% y >= poly$u))
})

test_that("polyhedra.cbsFt works when there are no splits", {
  set.seed(10)
  n <- 21
  y <- rnorm(n)
  obj <- circularBinSeg_fixedThres(y,2)
  poly <- polyhedra(obj)

  expect_true(class(poly) == "polyhedra")
  expect_true(all(dim(poly$gamma) == c(nrow(.enumerate_breakpoints_cbs(n))*2, n)))
  expect_true(all(poly$gamma %*% y >= poly$u))
})


test_that("polyhedra.cbsFt leads to uniform p values", {
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
    obj <- circularBinSeg_fixedThres(y,2)

    poly <- polyhedra(obj)
    
    if(all(is.na(jumps(obj)))){
      null_vec[i] <- runif(1)
    } else {
      contrast <- form_contrast(obj, n)
      null_vec[i] <- pvalue(y, poly, contrast)
    }
  }

  for(i in 1:trials){
    set.seed(i*10)
    y <- c(rnorm(n/3), rnorm(n/3, mean = 1), rnorm(n/3))
    obj <- circularBinSeg_fixedThres(y,2)

    poly <- polyhedra(obj)
    if(all(is.na(jumps(obj)))){
      alt_vec[i] <- runif(1)
    } else {
      contrast <- form_contrast(obj, n)
      alt_vec[i] <- pvalue(y, poly, contrast)
    }
  }

  quant <- seq(0, 1, length.out = 11)
  expect_true(sum(abs(quantile(null_vec, prob = quant) - quant)) <=
                sum(abs(quantile(alt_vec, prob = quant) - quant)) )

})