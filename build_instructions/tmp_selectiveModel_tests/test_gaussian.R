context("Test gaussian")

## .conditional_gaussian is correct

test_that(".conditional_gaussian works", {
  set.seed(10)
  cov_mat <- diag(10)
  cov_mat[1:5,1:5] <- 0.9
  cov_mat[6:10,6:10] <- 0.9
  diag(cov_mat) <- 1
  gaussian <- .gaussian(rep(0, 10), cov_mat)

  res <- .conditional_gaussian(gaussian, c(5,10))

  expect_true(class(res) == "gaussian")
  expect_true(length(res$mean) == 8)
  expect_true(all(dim(res$covariance) == 8))
})

# source: https://onlinecourses.science.psu.edu/stat505/node/43
test_that(".conditional_gaussian performs the correct calculation", {
  gaussian <- .gaussian(mean = c(175, 71), covariance = matrix(c(550, 40, 40, 8), nrow = 2))
  res <- .conditional_gaussian(gaussian, 70)

  expect_true(abs(res$mean - 170) < 1e-6)
  expect_true(abs(res$covariance - 350) < 1e-6)
})
