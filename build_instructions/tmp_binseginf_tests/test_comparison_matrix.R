context("Test comparison_matrix")

## .comparison_mat is correct

test_that(".comparison_mat enforces winning and losing requirements", {
  expect_error(.comparison_mat(1:2, matrix(1:15, ncol = 3, nrow = 5)))
  expect_error(.comparison_mat(1:3, matrix(1:15, ncol = 5, nrow = 3)))
})

test_that(".comparison_mat works", {
  res <- .comparison_mat(1:3, matrix(1:15, ncol = 3))
  
  expect_true(length(res) == 2)
  expect_true(all(names(res) == c("winning", "losing")))
  expect_true(all(res$winning == c(1:3)))
  expect_true(all(res$losing == matrix(1:15, ncol = 3)))
})

test_that(".comparison_mat enforces ordering on losing", {
  expect_error(.comparison_mat(1:3, matrix(15:1, ncol = 3, nrow = 5)))
  
  mat <-  matrix(1:15, ncol = 3, nrow = 5)
  mat2 <- mat
  mat2[2,1] <- 8
  expect_error(.comparison_mat(1:3, mat2))
  
  mat3 <- mat
  mat3[3,3] <- 4
  expect_error(.comparison_mat(1:3, mat3))
})

####################################

## .get_startEnd is correct

test_that(".get_startEnd works", {
  res <- .get_startEnd("1-10")
  expect_true(all(res == c(1,10)))
})

#####################################

## .get_leaves_matrix_excluding is correct

test_that(".get_leaves_matrix_excluding works", {
  set.seed(10)
  y <- c(rep(0, 20), rep(5, 10), rep(10, 5), rep(11, 5)) + 0.01*rnorm(40)
  obj <- binSeg_fixedSteps(y, 3)
  res <- .get_leaves_matrix_excluding(obj$tree, "1-20")
  
  expect_true(all(res == matrix(c(21, 30, 31, 35, 36, 40), nrow = 2, ncol =)))
})

test_that(".get_leaves_matrix_excluding errors if node is not a split", {
  set.seed(10)
  y <- c(rep(0, 20), rep(5, 10), rep(10, 5), rep(11, 5)) + 0.01*rnorm(40)
  obj <- binSeg_fixedSteps(y, 3)
  
  expect_error(.get_leaves_matrix_excluding(obj$tree, "1-23"))
})

test_that(".get_leaves_matrix_excluding returns NA if only one leaf", {
  obj <- .create_node(1, 30)
  
  expect_true(is.na(.get_leaves_matrix_excluding(obj, "1-30")))
})

test_that(".get_leaves_matrix_excluding works for singleton leaves", {
  set.seed(10)
  y <- c(-5, rep(0,9)) + 0.05*rnorm(10)
  
  obj <- binSeg_fixedSteps(y, 1)
  
  expect_true(is.na(.get_leaves_matrix_excluding(obj$tree, 
                                                 .get_leaves_names(obj$tree)[2])))
})
#########################################

## .threeColumnMatrix_from_nodeVec is correct

test_that(".threeColumnMatrix_from_nodeVec works", {
  res <- .threeColumnMatrix_from_nodeVec(c(1,5))
  
  expect_true(all(res == cbind(1, 1:4, 5)))
})

test_that(".threeColumnMatrix_from_nodeVec can exclude a value", {
  res <- .threeColumnMatrix_from_nodeVec(c(1,5), 4)
  
  expect_true(all(dim(res) == c(3,3)))
  expect_true(all(res == cbind(1, 1:3, 5)))
})

test_that(".threeColumnMatrix_from_nodeVec has a midpoint even if excluded", {
  res <- .threeColumnMatrix_from_nodeVec(c(1,2), 1)
  
  expect_true(all(dim(res) == c(1,3)))
  expect_true(all(res == c(1,1,2)))
})

############################################

## .threeColumnMatrix_from_nodeMatrix is correct

test_that(".threeColumnMatrix_from_nodeMatrix works", {
  mat <- matrix(1:10, ncol = 5, nrow = 2, byrow = T)
  res <- .threeColumnMatrix_from_nodeMatrix(mat)
  
  expect_true(length(res) == 5)
  expect_true(is.list(res))
  expect_true(all(sapply(res, ncol) == 3))
  expect_true(all(sapply(res, nrow) == 5))
  
  for(i in 1:5){
    expect_true(all(res[[i]] == cbind(i, i:(i+4), i+5)))
  }
})

#############################################

## .form_comparison is correct

test_that(".form_comparison works", {
  set.seed(10)
  y <- c(rep(0, 20), rep(5, 10), rep(6, 10)) + 0.01*rnorm(40)
  obj <- binSeg_fixedSteps(y, 1)
  
  obj2 <- binSeg_fixedSteps(y, 2)
  node <- .enumerate_splits(obj2$tree)[2]
  breakpoint <- data.tree::FindNode(obj2$tree, node)$breakpoint
  
  expect_true(node == "21-40")
  expect_true(breakpoint == 30)
  
  res <- .form_comparison(obj$tree, "21-40", 30)
  
  expect_true(all(names(res) == c("winning", "losing")))
  expect_true(all(res$winning == c(21, 30, 40)))
  expect_true(nrow(res$losing) == (40 - 21 - 1) + (20 - 1))
  
  expect_true(length(intersect(which(res$losing[,1] == 21), which(res$losing[,2] == 30))) == 0)
})

test_that(".form_comparison still works when there is no leaf", {
  tree <-  .create_node(1, 10)
  res <- .form_comparison(tree, "1-10", 5)
  
  expect_true(nrow(res$losing) == 8)
  mat <- cbind(1, c(1:9)[-5], 10)
  expect_true(all(res$losing == mat))
})

###################################

## .get_comparisonSigns is correct

test_that(".get_comparisonSigns works", {
  y <- c(rep(0, 5), rep(1, 5), rep(-1, 5))
  mat <- matrix(c(1, 5, 10, 6, 10, 15), ncol = 3, nrow = 2, byrow = T)
  
  res <- .get_comparisonSigns(y, mat)
  expect_true(all(res == c(1, -1)))
})