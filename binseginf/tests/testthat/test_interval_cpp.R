context("Test interval C++")
source("test_R_source.R")

####################################

## c_euclidean_to_radian is correct

test_that("c_euclidean_to_radian gives the correct answer", {
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    center <- rnorm(2)
    circle <- .circle(center = center, radius = sqrt(sum(center^2)))
    rad <- runif(1, -2*pi, 2*pi)
    point <- 2*(sin(rad)*center[1] + cos(rad)*center[2])*c(sin(rad), cos(rad))

    res1 <- .euclidean_to_radian(circle, point)
    res2 <- c_euclidean_to_radian_tester(center, sqrt(sum(center^2)),
                                          point)

    ifelse(abs(res1 - res2) < 1e-6, TRUE, FALSE)
  })

  expect_true(all(bool_vec))
})

test_that("c_euclidean_to_radian works", {
  point <- c(1,1+sqrt(2))
  res <- c_euclidean_to_radian_tester(c(1,1), sqrt(2), point)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
  expect_true(res <= pi/2)
  expect_true(res >= -pi/2)
})

test_that("c_euclidean_to_radian works with the origin", {
  point <- c(0,0)
  res <- c_euclidean_to_radian_tester(c(1,1), sqrt(2), point)

  expect_true(abs(2*(sin(res)*1+ cos(res)*1)) < 1e-6)
})

test_that("c_euclidean_to_radian returns the correct theta", {
  trials <- 100
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    center <- rnorm(2)
    rad <- runif(1, -2*pi, 2*pi)
    point <- 2*(sin(rad)*center[1] + cos(rad)*center[2])*c(sin(rad), cos(rad))

    res <- c_euclidean_to_radian_tester(center, sqrt(sum(center^2)), point)
    point2 <- 2*(sin(res)*center[1] + cos(res)*center[2])*c(sin(res), cos(res))

    ifelse(.l2norm(point - point2) < 1e-6, TRUE, FALSE)
  })

  expect_true(all(bool_vec))
})

test_that("c_euclidean_to_radian works with c_intersect_circle_line", {
  set.seed(10)
  a <- rnorm(2); b <- 0
  center <- rnorm(2)
  radius <- sqrt(sum(center^2))
  points <- c_intersect_circle_tester(a, b, center, radius)

  theta <- apply(points, 2, function(x){
    c_euclidean_to_radian_tester(center, radius, x)
  })

  bool_vec <- apply(points, 2, function(x){ #check to be on plane.
    abs(a %*%x - b) < 1e-16
  })
  expect_true(all(bool_vec))

  bool_vec <- apply(points, 2, function(x){ #check to be on plane
    abs(.l2norm(x - center) - radius) < 1e-6
  })
  expect_true(all(bool_vec))

  expect_true(all(theta <= pi/2))
  expect_true(all(theta >= -pi/2))
  expect_true(is.numeric(theta))
  expect_true(length(theta) == 2)
})

####################

## c_partition_interval is correct

test_that("c_partition_interval gives the correct answer", {
  trials <- 1000

  bool <- sapply(1:trials, function(x){
    set.seed(10*x)
    rad <- runif(3, min = 0, max = 2*pi)
    center <- rep(1/sqrt(2), 2)
    point_mat <- sapply(rad, function(x){
      c(cos(x), sin(x)) + center
    })
    circ <- .circle(center, 1)
    rad_vec <- apply(point_mat, 2, .euclidean_to_radian, circle = circ)
    interval <- .basic_interval(rad_vec[1:2], rad_vec[3])
    res1 <- .partition_interval(interval)

    res2 <- c_partition_interval(interval)

    all(res1 == res2)
  })

  expect_true(all(bool))
})


test_that(".partition_interval works", {
  res <- c_partition_interval(c(pi/4, pi/3))

  expect_true(is.matrix(res))
  expect_true(is.numeric(res))
  expect_true(ncol(res) == 2)
})

test_that(".partition_interval works when endpoints are different signs, no wrap-around", {
  res <- c_partition_interval(c(-pi/3, pi/4))

  expect_true(all(dim(res) == c(2,2)))
  expect_true(all(res[1,] == c(-pi/3, 0)))
  expect_true(all(res[2,] == c(0, pi/4)))
})

test_that(".partition_interval works when endpoints are different signs, yes wrap-around", {
  res <- c_partition_interval(c(pi/3, 3*pi/4))

  expect_true(all(dim(res) == c(2,2)))
  expect_true(all(res[1,] == c(-pi/2, -pi/4)))
  expect_true(all(res[2,] == c(pi/3, pi/2)))

  res <- c_partition_interval(c(-3*pi/4, -pi/3))

  expect_true(all(dim(res) == c(2,2)))
  expect_true(all(res[1,] == c(-pi/2, -pi/3)))
  expect_true(all(res[2,] == c(pi/4, pi/2)))
})

test_that(".partition_interval works when endpoints are same signs, yes wrap-around", {
  res <- c_partition_interval(c(-7*pi/6, -pi/3))

  expect_true(all(dim(res) == c(3,2)))
  expect_true(sum(abs(res[1,] - c(-pi/2, -pi/3))) < 1e-6)
  expect_true(sum(abs(res[2,] - c(-pi/6, 0))) < 1e-6)
  expect_true(sum(abs(res[3,] - c(0, pi/2))) < 1e-6)

  res <- c_partition_interval(c(pi/3, 7*pi/6))
  expect_true(all(dim(res) == c(3,2)))
  expect_true(sum(abs(res[1,] - c(-pi/2, 0))) < 1e-6)
  expect_true(sum(abs(res[2,] - c(0, pi/6))) < 1e-6)
  expect_true(sum(abs(res[3,] - c(pi/3, pi/2))) < 1e-6)
})

####################

## c_initial_theta is correct

# test_that("c_initial_theta works", {
#   set.seed(10)
#   y <- rnorm(10)
#   v <- rnorm(10); w <- rnorm(10)
#   v <- v/.l2norm(v)
#   w <- .projection(w, v); w <- w/.l2norm(w)
# 
#   res <- c_initial_theta(y, v, w)
# 
#   expect_true(is.numeric(res))
#   expect_true(!is.matrix(res))
#   expect_true(length(res) == 1)
#   expect_true(res <= pi/2)
#   expect_true(res >= -pi/2)
# })

# test_that(".initial_theta gives a radius of 0", {
#   set.seed(10)
#   y <- rnorm(10)
#   v <- rnorm(10); w <- rnorm(10)
#   v <- v/.l2norm(v)
#   w <- .projection(w, v); w <- w/.l2norm(w)
#   theta <- c_initial_theta(y, v, w)
# 
#   radius <- .radius(theta, y, v, w)
# 
#   expect_true(abs(radius) < 1e-6)
# })

# test_that(".initial_theta gives the proper theta", {
#   set.seed(10)
#   y <- rnorm(10)
#   v <- rnorm(10); w <- rnorm(10)
#   v <- v/.l2norm(v)
#   w <- .projection(w, v); w <- w/.l2norm(w)
#   theta <- c_initial_theta(y, v, w)
# 
#   res <- .radians_to_data(theta, y, v, w)
# 
#   expect_true(sum(abs(y - res)) < 1e-6)
# })

#################

## c_interval is correct

test_that("c_interval gives the correct answer", {
  trials <- 1000

  bool <- sapply(1:trials, function(x){
    set.seed(11*x)
    rad <- runif(3, min = 0, max = 2*pi)
    center <- rep(1/sqrt(2), 2)
    point_mat <- sapply(rad, function(x){
      c(cos(x), sin(x)) + center
    })
    circ <- .circle(center, 1)
    rad_vec <- apply(point_mat, 2, .euclidean_to_radian, circle = circ)
    res1 <- .interval(rad_vec[1:2], rad_vec[3])

    res2 <- c_interval(rad_vec[1:2], rad_vec[3])

    all(res1 == res2)
  })

  expect_true(all(bool))
})

###################

## c_form_interval is correct

# test_that("c_form_interval works on one particular case", {
#   set.seed(10)
#   y <- rnorm(10)
#   obj <- binseginf::bsfs(y, 2)
#   poly <- binseginf::polyhedra(obj)
# 
#   v <- rnorm(10); w <- rnorm(10)
#   v <- v/.l2norm(v)
#   w <- .projection(w, v); w <- w/.l2norm(w)
# 
#   res <- c_form_interval(poly$gamma[1,], poly$u[1], y, v, w)
# 
#   expect_true(is.matrix(res))
#   expect_true(ncol(res) == 2)
# })

# test_that("c_form_interval is correct", {
#   set.seed(10)
#   y <- rnorm(10)
#   obj <- binseginf::bsfs(y, 2)
#   poly <- binseginf::polyhedra(obj)
# 
#   v <- rnorm(10); w <- rnorm(10)
#   v <- v/.l2norm(v)
#   w <- .projection(w, v); w <- w/.l2norm(w)
# 
#   res_list1 <- sapply(1:length(poly$u), function(x){
#     plane <- .plane(poly$gamma[x,], poly$u[x])
#     plane <- .intersect_plane_basis(plane, y, v, w)
#     if(any(is.na(plane))) return(matrix(c(-pi/2, pi/2), ncol = 2))
#     center <- c(-y%*%v, -y%*%w)
#     radius <- sqrt(sum(center^2))
#     circle <- .circle(center, radius)
#     dis <- .distance_point_to_plane(center, plane)
# 
#     if(dis >= radius){
#       matrix(c(-pi/2, pi/2), ncol = 2)
#     } else {
#       mat <- .intersect_circle_line(plane, circle)
#       stopifnot(nrow(mat) == 2)
#       vec <- apply(mat, 2, .euclidean_to_radian, circle = circle)
#       init_theta <- .initial_theta(y, v, w)
#       .interval(vec, init_theta)
#     }
#   })
# 
#   res_list2 <- sapply(1:length(poly$u), function(x){
#     c_form_interval(poly$gamma[x,], poly$u[x], y, v, w)
#   })
# 
#   expect_true(length(res_list1) == length(res_list2))
# 
#   bool <- sapply(1:length(res_list1), function(x){
#     sum(abs(res_list1[[x]] - res_list2[[x]])) < 1e-6
#   })
# 
#   expect_true(all(bool))
# })


