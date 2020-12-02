context("Test plane C++")
source("test_R_source.R")

#########################

test_that("Plane can be generated properly", {
  trials <- 100
  bool <- sapply(1:trials, function(x){
    set.seed(10*x)
    a <- rnorm(10); b <- rnorm(1)

    res1 <- .plane(a = a, b = b)
    res2 <- new(Plane, a, b)

    sum(abs(res1$a - res2$a)) <= 1e-6 & abs(res1$b - res2$b) <= 1e-6
  })

  expect_true(all(bool))
})

########

test_that("c_intersect_plane_basis works properly", {
  trials <- 100
  bool <- sapply(1:trials, function(x){
    set.seed(10*x)
    a <- rnorm(10); b <- rnorm(1)
    y <- rnorm(10); v <- rnorm(10); w <- rnorm(10)

    res1 <- .plane(a = a, b = b)
    res1 <- .intersect_plane_basis(res1, y, v, w)

    res2 <- new(Plane, a, b)
    res2$c_intersect_basis(y, v, w)

    sum(abs(res1$a - res2$a)) <= 1e-6 & abs(res1$b - res2$b) <= 1e-6
  })

  expect_true(all(bool))
})

############

test_that("c_point_on_plane works properly", {
  trials <- 100
  bool <- sapply(1:trials, function(x){
    set.seed(10*x)
    a <- rnorm(10); b <- rnorm(1)

    plane <- new(Plane, a, b)
    res <- plane$c_point_on_plane()

    sum(abs(res %*% plane$a - plane$b)) <= 1e-6
  })

  expect_true(all(bool))
})

###########

test_that("c_distance_point_to_plane works properly", {
  trials <- 100
  bool <- sapply(1:trials, function(x){
    set.seed(10*x)
    a <- rnorm(10); b <- rnorm(1)
    point <- rnorm(10)

    plane1 <- .plane(a = a, b = b)
    res1 <- .distance_point_to_plane(point, plane1)

    plane2 <- new(Plane, a, b)
    res2 <- plane2$c_distance_point_to_plane(point)

    abs(res1 - res2) <= 1e-6
  })

  expect_true(all(bool))
})

#################

test_that("c_intersect_circle works properly", {
  trials <- 100

  bool <- sapply(1:trials, function(x){
    set.seed(10*x)
    a <- rnorm(2); b <- 1
    plane <- .plane(a, b)
    center <- rnorm(2); radius <- 5
    circle <- .circle(center, radius)
    dis <- .distance_point_to_plane(center, plane)

    if(dis >= radius){
      return(TRUE)
    } else {
      res1 <- .intersect_circle_line(plane, circle, full = T)

      res2 <- c_intersect_circle_tester(a, b, center, radius)

      sum(abs(res1 - res2)) < 1e-6
    }
  })

  expect_true(all(bool))
})

test_that("c_intersect_circle works", {
  res <- c_intersect_circle_tester(c(0,1), 0, c(0,0), 1)

  expect_true(is.matrix(res))
  expect_true(is.numeric(res))
  expect_true(all(dim(res) == c(2,2)))
})

test_that("c_intersect_circle gives a proper point", {
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    a <- rnorm(2)
    points <- c_intersect_circle_tester(a, 0, c(0,0), 1)

    #check to see points are in circle
    res1 <- apply(points, 2, function(x){sum(abs(sqrt(x[1]^2+x[2]^2) - 1)) < 1e-6})

    #check to see points are on plane
    res2 <- apply(points, 2, function(x){abs(a%*%x) < 1e-6})

    all(c(res1, res2))
  })

  expect_true(all(bool_vec))
})

test_that("c_intersect_circle gives a proper point with b", {
  trials <- 100
  radius <- 100
  b <- 3
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    a <- rnorm(2)
    points <- c_intersect_circle_tester(a, b, c(0,0), radius)

    #check to see points are in circle
    res1 <- apply(points, 2, function(x){sum(abs(sqrt(x[1]^2+x[2]^2) - radius)) < 1e-6})

    #check to see points are on plane
    res2 <- apply(points, 2, function(x){abs(a%*%x - b) < 1e-6})

    all(c(res1, res2))
  })

  expect_true(all(bool_vec))
})

test_that("c_intersect_circle gives a proper point with small values of a[1]", {
  trials <- 100
  radius <- 100
  b <- 2
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    a <- c(0, rnorm(1))
    points <- c_intersect_circle_tester(a, b, c(0,0), radius*b)

    if(!is.matrix(points)) points <- as.matrix(points, nrow = 1)

    #check to see points are in circle
    res1 <- apply(points, 2, function(x){sum(abs(sqrt(x[1]^2+x[2]^2) - b*radius)) < 1e-6})

    #check to see points are on plane
    res2 <- apply(points, 2, function(x){abs(a%*%x - b) < 1e-6})

    all(c(res1, res2))
  })

  expect_true(all(bool_vec))
})

test_that(".intersect_circle_line gives two points with small values of a[1] always gives 2 points", {
  trials <- 100
  radius <- 100
  b <- 2
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    a_vec <- c(0, rnorm(1))
    points <- c_intersect_circle_tester(a_vec, b, c(0,0), radius*b)

    ncol(points) == 2
  })

  expect_true(all(bool_vec))
})

test_that(".intersect_circle_line gives a proper point with small values of a[1] and offcenter circle", {
  trials <- 100
  radius <- 100
  b <- 2
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    a <- c(0, rnorm(1))
    center <- c(1,2)
    points <- c_intersect_circle_tester(a, b, center, radius*b)

    if(!is.matrix(points)) points <- as.matrix(points, nrow = 1)

    #check to see points are in circle
    res1 <- apply(points, 2, function(x){sum(abs(sqrt((x[1]-center[1])^2+(x[2]-center[2])^2) - radius*b)) < 1e-6})

    #check to see points are on plane
    res2 <- apply(points, 2, function(x){abs(a%*%x - b) < 1e-6})

    all(c(res1, res2))
  })

  expect_true(all(bool_vec))
})

test_that(".intersect_circle_line gives a proper point with small values of a[2] and offcenter circle", {
  trials <- 100
  radius <- 100
  b <- 2
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    a <- c(rnorm(1), 0)
    center <- c(1,2)
    points <- c_intersect_circle_tester(a, b, center, radius*b)

    if(!is.matrix(points)) points <- as.matrix(points, nrow = 1)

    #check to see points are in circle
    res1 <- apply(points, 2, function(x){sum(abs(sqrt((x[1]-center[1])^2+(x[2]-center[2])^2) - radius*b)) < 1e-6})

    #check to see points are on plane
    res2 <- apply(points, 2, function(x){abs(a%*%x - b) < 1e-6})

    all(c(res1, res2))
  })

  expect_true(all(bool_vec))
})

test_that(".intersect_circle_line can give NA", {
  res <- c_intersect_circle_tester(c(1, -1), 10, c(0,0), 1)

  expect_true(all(dim(res) == c(2,2)))
  expect_true(all(is.na(res)))
})

test_that(".intersect_circle_line can give one point", {
  res <- c_intersect_circle_tester(c(1, 0), 1, c(0,0), 1)

  expect_true(is.matrix(res))
  expect_true(is.numeric(res))
  expect_true(all(dim(res) == c(2,2)))
  expect_true(all(is.na(res[,2])))
  expect_true(all(res[,1] == c(1,0)))
})

#############

## c_intersect_basis is correct

test_that("c_intersect_basis is correct", {
  trials <- 100
  bool <- sapply(1:trials, function(x){
    set.seed(10*x)
    a <- rnorm(10); b <- rnorm(1);
    y <- rnorm(10); v <- rnorm(10); w <- rnorm(10)

    plane <- .plane(a, b)
    res1 <- .intersect_plane_basis(plane, y, v, w)
    res2 <- c_intersect_basis_tester(a, b, y, v, w)

    sum(abs(res1$a - res2$a)) < 1e-6 & abs(res1$b - res2$b) < 1e-6
  })

  expect_true(all(bool))
})

test_that("c_intersect_basis can give NA", {
  res <- c_intersect_basis_tester(c(1,0,0), 1, 1:3, c(0,1,0), c(0,0,1))

  expect_true(all(is.na(res$a)))
  expect_true(is.na(res$b))
})


