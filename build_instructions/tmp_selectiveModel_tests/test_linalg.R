context("Test linear algebra")

## .projection is correct

test_that(".projection works", {
  set.seed(10)
  vec1 <- rnorm(10); vec2 <- rnorm(10)
  res <- .projection(vec1, vec2)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 10)
})

test_that(".projection preserves orthogonal vectors", {
  set.seed(1)
  vec1 <- rnorm(10); vec2 <- rnorm(10)
  res <- .projection(vec1, vec2)
  res2 <- .projection(res, vec2)

  expect_true(sum(abs(res - res2)) < 1e-6)
})

test_that(".projection removes parallel vectors", {
  set.seed(1)
  vec1 <- rnorm(10); vec2 <- 2*vec1
  res <- .projection(vec1, vec2)

  expect_true(sum(abs(res)) < 1e-6)
})

test_that(".projection reduces the norm", {
  set.seed(1)
  vec1 <- rnorm(10); vec2 <- rnorm(10)
  res <- .projection(vec1, vec2)

  expect_true(.l2norm(vec1) >= .l2norm(res))
})

################################

## .projection_matrix is correct

test_that(".projection_matrix works", {
  set.seed(10)
  mat <- .segments(10, c(3, 7))
  vec <- rnorm(10)
  res <- .projection_matrix(vec, mat)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 10)
})

test_that(".projection_matrix returns a vector that is orthogonal to rows of mat", {
  set.seed(10)
  mat <- .segments(10, c(3, 7))
  vec <- rnorm(10)
  res <- .projection_matrix(vec, mat)

  for(i in 1:nrow(mat)){
    expect_true(abs(res %*% mat[i,]) < 1e-6)
  }
})

test_that(".projection_matrix reduces the norm", {
  set.seed(10)
  mat <- matrix(rnorm(30), ncol = 10, nrow = 3)
  vec <- rnorm(10)
  res <- .projection_matrix(vec, mat)

  expect_true(.l2norm(vec) >= .l2norm(res))
})

test_that(".projection_matrix preserves orthogonal vectors", {
  set.seed(10)
  mat <- matrix(rnorm(30), ncol = 10, nrow = 3)
  vec <- rnorm(10)
  res <- .projection_matrix(vec, mat)
  res2 <- .projection_matrix(res, mat)

  expect_true(sum(abs(res - res2)) < 1e-6)
})

test_that(".projection_matrix removes parallel vectors", {
  set.seed(10)
  mat <- matrix(rnorm(30), ncol = 10, nrow = 3)
  vec <- 2*mat[1,]
  res <- .projection_matrix(vec, mat)

  expect_true(sum(abs(res)) < 1e-6)
})

#########################

## .sample_matrix_space is correct

test_that(".sample_matrix_space works", {
  set.seed(10)
  mat <- .segments(10, c(3, 7))
  res <- .sample_matrix_space(mat, 2)

  expect_true(is.numeric(res))
  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(10,2)))
})

test_that(".sample_matrix_space gives vectors that are orthogonal to mat", {
  set.seed(10)
  mat <- .segments(10, c(3, 7))
  res <- .sample_matrix_space(mat, 3)

  for(i in 1:3){
    for(j in 1:3){
      expect_true(abs(res[,i]%*%mat[j,]) < 1e-6)
    }
  }
})

test_that(".sample_matrix_space gives vectors that are orthogonal when num_vec is not 1", {
  set.seed(10)
  mat <- .segments(10, c(3, 7))
  vec_mat <- .sample_matrix_space(mat, 3)
  res <- t(vec_mat)%*%vec_mat

  expect_true(sum(abs(res - diag(3))) < 1e-6)
})

test_that(".sample_matrix_space gives vectors that are orthogonal", {
  set.seed(10)
  mat <- .segments(10, c(3, 7))
  vec_mat <- .sample_matrix_space(mat, NA)
  res <- t(vec_mat)%*%vec_mat

  expect_true(sum(abs(res - diag(7))) < 1e-6)
})

test_that(".sample_matrix_space samples from rowspace properly", {
  set.seed(10)
  mat <- .segments(10, c(4,8))
  res <- .sample_matrix_space(mat, null = F)

  #check to see if nothing remains when projecting against res
  for(i in 1:3){
    tmp <- .projection_matrix(mat[i,], t(res))
    expect_true(sum(abs(tmp)) < 1e-6)
  }
})

test_that(".sample_matrix_space similary gives proper outputs for rowspace", {
  set.seed(10)
  mat <- .segments(10, c(1,3,9))
  res <- .sample_matrix_space(mat, null = F)

  expect_true(all(dim(res) == c(10,4)))

  tmp <- t(res) %*% res
  expect_true(sum(abs(tmp - diag(4))) < 1e-6)
})

test_that(".sample_matrix_space is properly uniform", {
  set.seed(10)
  trials <- 100
  mat <- .segments(3, 1, ignore_jump = 1)
  v_mat <- sapply(1:trials, function(x){
    .sample_matrix_space(mat, num_vec = 1, null = T)
  })

  basis1 <- c(1,1,1); basis1 <- basis1/.l2norm(basis1)
  basis2 <- c(0,-1,1); basis2 <- basis2/.l2norm(basis2)
  basis3 <- c(-2,1,1); basis3 <- basis3/.l2norm(basis3)
  basis <- rbind(basis1, basis2, basis3)

  v_mat <- apply(v_mat, 2, function(x){
    solve(t(basis), x)[2:3]
  })

  angle_vec <- apply(v_mat, 2, function(x){
    if(x[1] >= 0 & x[2] >= 0) {
      atan(x[2]/x[1])
    } else if (x[1] <= 0 & x[2] >= 0){
      atan(x[2]/x[1])+pi
    } else if (x[1] <= 0 & x[2] <= 0){
      atan(x[2]/x[1])+pi
    } else {
      atan(x[2]/x[1])+2*pi
    }
  })

  angle_vec2 <- rnorm(trials, mean = pi, sd = pi/2)

  expect_true(sum(abs(sort(angle_vec) - seq(0, 2*pi, length.out = trials))) <
                sum(abs(sort(angle_vec) - sort(angle_vec2))))
})

test_that(".sample_matrix_space gives the proper segments_full", {
  n <- 10
  set.seed(10)
  y <- rnorm(n)
  fit <- binseginf::bsfs(y, 1)
  polyhedra <- binseginf::polyhedra(fit)
  gaussian <- .gaussian(rep(0, n), diag(n))

  jump <- binseginf::jumps(fit)
  segments <- .segments(n, jump)
  nullspace_mat <- .sample_matrix_space(segments)
  mean_val <- as.numeric(segments%*%y)
  segments_full <- rbind(t(nullspace_mat), segments)

  res <- segments_full %*% y

  expect_true(abs(res[9] - mean(y[1:jump])) < 1e-6)
  expect_true(abs(res[10] - mean(y[(jump+1):n])) < 1e-6)
})
