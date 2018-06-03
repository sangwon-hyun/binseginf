context("Test comparison matrix for circular binseg")

## .form_comparison_cbs is correct

test_that(".form_comparison_cbs works", {
  y <- c(rep(0,10), rep(1,10))
  obj <- circularBinSeg_fixedSteps(y, 1)
  
  res <- .form_comparison_cbs(obj$tree, "1-10", c(3,5))
  
  expect_true(is.list(res))
  expect_true(length(res) == 2)
  expect_true(all(sort(names(res)) == c("losing", "winning")))
  expect_true(all(sapply(res, ncol) == 4))
})

test_that(".form_comparison_cbs has the right number of rows", {
  y <- c(rep(0,5), rep(1,5))
  obj <- circularBinSeg_fixedSteps(y, 1)
  
  res <- .form_comparison_cbs(obj$tree, "1-5", c(3,5))
  
  expect_true(all(as.numeric(res$winning) == c(1,3,5,5)))
  expect_true(nrow(res$losing) == c(2*nrow(.enumerate_breakpoints_cbs(5))-1))
})
  
test_that(".form_comparison_cbs gives unique rows", {
  y <- c(rep(0,5), rep(1,5))
  obj <- circularBinSeg_fixedSteps(y, 1)
  
  res <- .form_comparison_cbs(obj$tree, "1-5", c(3,5))
  
  losing <- apply(res$losing, 1, function(x){x[1]*11^3 + x[2]*11^2 + x[3]*11 + x[4]})
  expect_true(length(unique(losing)) == length(losing))
})

test_that(".form_comparison_cbs can set excluding to false", {
  y <- c(rep(0,5), rep(1,5))
  obj <- circularBinSeg_fixedSteps(y, 1)
  
  res <- .form_comparison_cbs(obj$tree, "1-5", c(3,5), excluding = F)
  
  expect_true(nrow(res$losing) == nrow(.enumerate_breakpoints_cbs(5))-1)
})

test_that(".form_comparison_cbs works for vectors of length 2", {
  y <- c(0,1)
  obj <- circularBinSeg_fixedThres(y, 5)
  
  res <- .form_comparison_cbs(obj$tree, "1-2", c(2,2), excluding = F)
  
  expect_true(ncol(res$losing) == 4)
  expect_true(all(is.na(res$losing)))
})

######################

## .fourColumnMatrix_from_nodeMatrix is correct

test_that(".fourColumnMatrix_from_nodeMatrix works", {
  mat <- matrix(c(1,5,3,7), 2, 2)
  res <- .fourColumnMatrix_from_nodeMatrix(mat)
  
  expect_true(is.list(res))
  expect_true(length(res) == 2)
  expect_true(all(sapply(res, ncol)== 4))
})

################

## .fourColumnMatrix_from_nodeVec is correct

test_that(".fourColumnMatrix_from_nodeVec works", {
  vec <- c(3,8)
  res <- .fourColumnMatrix_from_nodeVec(vec)
  
  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(nrow(.enumerate_breakpoints_cbs(6)), 4)))
})

test_that(".fourColumnMatrix_from_nodeVec allows exclude", {
  vec <- c(3,8)
  res <- .fourColumnMatrix_from_nodeVec(vec, exclude = c(4,7))

  expect_true(nrow(res) == nrow(.enumerate_breakpoints_cbs(6)) - 1)
  expect_true(!any(apply(res, 1, function(x){ifelse(x[2] == 4 & x[3] == 7, TRUE, FALSE)})))
  
  res <- .fourColumnMatrix_from_nodeVec(vec)
  bool <- apply(res, 1, function(x){ifelse(x[2] == 4 & x[3] == 7, TRUE, FALSE)})
  expect_true(length(which(bool)) == 1)
})

test_that(".fourColumnMatrix_from_nodeVec returns all NA's if vec[2] == vec[1]+1", {
  res <- .fourColumnMatrix_from_nodeVec(c(2,3),c(3,3))
  expect_true(all(is.na(res)))
  expect_true(all(dim(res) == c(1,4)))
  
  res <- .fourColumnMatrix_from_nodeVec(c(2,3),c(2,2))
  expect_true(all(is.na(res)))
  expect_true(all(dim(res) == c(1,4)))
  
  res <- .fourColumnMatrix_from_nodeVec(c(2,3),NA)
  expect_true(res[1,1] == 2)
  expect_true(res[1,4] == 3)
})