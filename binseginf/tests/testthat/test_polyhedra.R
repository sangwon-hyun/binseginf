context("Test polyhedra")


test_that("polyhedra forms a correct class", {
  res <- polyhedra(matrix(1:25,5,5), rep(1,5))

  expect_true(class(res) == "polyhedra")
})



test_that("is_valid.polyhedra works", {
  res <- polyhedra(matrix(1:25,5,5), rep(1,5))

  expect_true(is_valid(res))

  res$u <- rep(1,4)

  expect_error(is_valid(res))
})


test_that("combine.polyhedra works on lists containing empty polyhedra", {

    ## Check that it handles empty polyhedra well.
    g1 = matrix(runif(20),nrow=2)
    u1 = runif(2)
    g2 = rbind(rep(NA,10))[-1,]
    u2 = c()
    b = polyhedra.matrix(g1,u1)
    a = polyhedra.matrix(g2,u2)

    correctgamma = rbind(g1,g1,g2,g1)
    correctu = c(u1,u1,u2,u1)

    expect_equal(correctgamma, combine(b,b,a,b)$gamma)
})



