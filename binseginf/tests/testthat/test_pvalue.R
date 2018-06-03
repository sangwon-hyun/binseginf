context("Test pvalues from truncated Gaussians")

## pvalue is correct

test_that("p value is high power for correct changepoint", {
  set.seed(10)
  y <- c(rep(0, 10), rep(50, 10)) + rnorm(20)
  obj <- binSeg_fixedSteps(y, 1)
  
  poly <- polyhedra(obj)
  contrast <- contrast_vector(obj, 1)
  
  res <- pvalue(y, poly, contrast)
  expect_true(res < .1)
})

test_that("p value are roughly uniform", {
  len <- 250
  pvalue_null.vec <- numeric(len)
  
  for(i in 1:len){
    set.seed(i*10)
    y <- rnorm(20)
    obj <- binSeg_fixedSteps(y, 1)
    
    poly <- polyhedra(obj)
    contrast <- contrast_vector(obj, 1)
    
    pvalue_null.vec[i] <- pvalue(y, poly, contrast)
  }
  
  expect_true(sum(abs(sort(pvalue_null.vec) - seq(0, 1, length.out = 250))) <= 7)
})

test_that("p value are roughly uniform compared to alt", {
  len <- 50
  pvalue_null.vec <- numeric(len)
  pvalue_alt.vec <- numeric(len)
  
  for(i in 1:len){
    set.seed(i*10)
    y <- rnorm(20)
    obj <- binSeg_fixedSteps(y, 1)
    
    poly <- polyhedra(obj)
    contrast <- contrast_vector(obj, 1)
    
    pvalue_null.vec[i] <- pvalue(y, poly, contrast)
  }
  
  for(i in 1:len){
    set.seed(i*10)
    y <-  c(rep(0, 10), rep(3, 10)) + rnorm(20)
    obj <- binSeg_fixedSteps(y, 1)
    
    poly <- polyhedra(obj)
    contrast <- contrast_vector(obj, 1)
    
    pvalue_alt.vec[i] <- pvalue(y, poly, contrast)
  }
  
  quant <- c(0, 0.25, 0.5, 0.75, 1)
  expect_true(sum(abs(quantile(pvalue_null.vec, probs = quant, na.rm = T) - quant))
    < sum(abs(quantile(pvalue_alt.vec, probs = quant, na.rm = T) - quant)))
})

test_that("p value one-sided and two-sided are related for binseg", {
  set.seed(10)
  y <- c(rep(0, 50), rep(0.1, 50)) + rnorm(100)
  obj <- binSeg_fixedSteps(y, 1)
  
  poly <- polyhedra(obj)
  contrast <- contrast_vector(obj, 1)
  
  res.pos.onesided <- pvalue(y, poly, contrast)
  res.two.sided <- pvalue(y, poly, contrast, alternative = "two.sided")
  res.neg.onesided <- pvalue(y, poly, -contrast)
  
  expect_true(contrast%*%y > 0)
  expect_true(res.neg.onesided > res.pos.onesided)
  expect_true(abs(res.pos.onesided*2 - res.two.sided) < 1e-4)
})

test_that("p value one-sided and two-sided are related for fused lasso", {
  set.seed(10)
  y <- c(rep(0, 50), rep(0.1, 50)) + rnorm(100)
  obj <- fLasso_fixedSteps(y, 1)
  
  poly <- polyhedra(obj)
  contrast <- contrast_vector(obj, 1)
  
  res.pos.onesided <- pvalue(y, poly, contrast)
  res.two.sided <- pvalue(y, poly, contrast, alternative = "two.sided")
  res.neg.onesided <- pvalue(y, poly, -contrast)
  
  expect_true(contrast%*%y > 0)
  expect_true(res.neg.onesided > res.pos.onesided)
  expect_true(abs(res.pos.onesided*2 - res.two.sided) < 1e-4)
})

test_that("pvalue is not 1 when the signal is extremely large", {
  set.seed(1)
  y <- CpVector(100, c(0,5), 0.5)$data
  
  obj <- binSeg_fixedSteps(y, 1)
  
  poly <- polyhedra(obj, y)
  contrast <- contrast_vector(obj, 1)
  
  res <- pvalue(y, poly, contrast)
  expect_true(res < 1e-4)
})

test_that("p value is correct in reverse for binSeg", {
  set.seed(10)
  y <- c(rep(0, 10), rep(-50, 10)) + rnorm(20)
  obj <- binSeg_fixedSteps(y, 1)
  
  poly <- polyhedra(obj)
  contrast <- contrast_vector(obj, 1)
  
  res <- pvalue(y, poly, contrast)
  expect_true(res < .1)
})

test_that("p value is correct in reverse for fLasso", {
  set.seed(10)
  y <- c(rep(0, 10), rep(-50, 10)) + rnorm(20)
  obj <- fLasso_fixedSteps(y, 1)
  
  poly <- polyhedra(obj)
  contrast <- contrast_vector(obj, 1)
  
  res <- pvalue(y, poly, contrast)
  expect_true(res < .1)
})

test_that("p value is correct for a bump signal", {
  set.seed(10)
  y <- c(rep(0, 10), rep(50, 10), rep(0, 10)) + rnorm(30)
  obj <- binSeg_fixedSteps(y, 2)
  poly <- polyhedra(obj)
  
  contrast1 <- contrast_vector(obj, 1, sorted = T)
  res1 <- pvalue(y, poly, contrast1)
  contrast2 <- contrast_vector(obj, 2, sorted = T)
  res2 <- pvalue(y, poly, contrast2)
  
  expect_true(all(c(res1, res2) < .1))
})

############################

## .truncated_gauss_cdf is correct

test_that(".truncated_gauss_cdf gives 1 when out of bounds", {
  res <- .truncated_gauss_cdf(10, 0, 1, 0, 9)
  expect_true(res == 1)
})

test_that(".truncated_gauss_cdf returns 0 for extreme values", {
  res <- .truncated_gauss_cdf(10, 0, 1, 5, Inf)
  expect_true(res == 0)
})

###################################

## .compute_truncGaus_terms is correct for binSeg

test_that(".compute_truncGaus_terms preserves vlo and vup correctly for binseg", {
  set.seed(5)
  y <- c(rep(0,50), rep(-2,20), rep(-1,30)) + rnorm(100)
  obj <- binSeg_fixedSteps(y,2)

  poly <- polyhedra(obj)
  contrast <- contrast_vector(obj, 1)
  
  res <- .compute_truncGaus_terms(y, poly, contrast, 1)
  expect_true(res$a <= contrast%*%y & 
                res$b >= contrast%*%y)
  
  trials <- 100
  for(i in 1:trials){
    set.seed(i*10)
    y.tmp <- c(rep(0,50), rep(-2,20), rep(-1,30)) + rnorm(100)
    res2 <- .compute_truncGaus_terms(y.tmp, poly, contrast, 1)
    
    bool1 <- (res2$a <= contrast%*%y.tmp & 
                res2$b >= contrast%*%y.tmp)
    bool2 <- all(poly$gamma %*% y.tmp >= poly$u)
    expect_true(bool1 == bool2)
  }
})

test_that(".compute_truncGaus_terms preserves vlo and vup correctly for flasso", {
  set.seed(5)
  y <- c(rep(0,50), rep(-2,20), rep(-1,30)) + rnorm(100)
  obj <- fLasso_fixedSteps(y,2)
  
  poly <- polyhedra(obj)
  contrast <- contrast_vector(obj, 1)
  
  res <- .compute_truncGaus_terms(y, poly, contrast, 1)
  expect_true(res$a <= contrast%*%y & 
                res$b >= contrast%*%y)
  
  trials <- 100
  for(i in 1:trials){
    set.seed(i*10)
    y.tmp <- c(rep(0,50), rep(-2,20), rep(-1,30)) + rnorm(100)
    res2 <- .compute_truncGaus_terms(y.tmp, poly, contrast, 1)
    
    bool1 <- (res2$a <= contrast%*%y.tmp & 
                res2$b >= contrast%*%y.tmp)
    bool2 <- all(poly$gamma %*% y.tmp >= poly$u)
    expect_true(bool1 == bool2)
  }
})