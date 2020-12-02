context("Test declutter")

## declutter is correct

test_that("declutter works", {
  jump_vec <- c(40,41,42,50,52,60,63,80,120,121,122)
  sign_vec <- rep(1,length(jump_vec))
  sign_vec[2] <- -1

  res <- declutter(jump_vec, sign_vec)

  expect_true(is.list(res))
  expect_true(length(res) == 3)
  expect_true(all(names(res) == c("jump_vec", "sign_mat", "target_bool")))
})

test_that("declutter works with all singletons", {
  set.seed(10)
  jump_vec <- seq(1,100, by = 10)
  sign_vec <- sample(c(-1,1), length(jump_vec), replace = T)

  res <- declutter(jump_vec, sign_vec)

  expect_true(is.list(res))
  expect_true(length(res) == 3)
  expect_true(all(names(res) == c("jump_vec", "sign_mat", "target_bool")))
})

test_that("declutter works with one big block", {
  set.seed(10)
  jump_vec <- 1:10
  sign_vec <- sample(c(-1,1), length(jump_vec), replace = T)

  res <- declutter(jump_vec, sign_vec)

  expect_true(is.list(res))
  expect_true(length(res) == 3)
  expect_true(all(names(res) == c("jump_vec", "sign_mat", "target_bool")))
})

test_that("declutter works for one jump", {
  jump_vec <- 5
  sign_vec <- 1

  res <- declutter(jump_vec, sign_vec)

  expect_true(is.list(res))
  expect_true(length(res) == 3)
  expect_true(all(names(res) == c("jump_vec", "sign_mat", "target_bool")))
})

test_that("declutter gets smaller if how_close increases", {
  trials <- 100

  bool_vec <- sapply(1:trials, function(x){
    set.seed(10*x)
    jump_vec <- sort(sample(1:100, 10))
    sign_vec <- sample(c(-1,1), length(jump_vec), replace = T)

    res1 <- declutter(jump_vec, sign_vec, how_close = 1)
    res2 <- declutter(jump_vec, sign_vec, how_close = 10)

    length(res1$jump_vec) >= length(res2$jump_vec)
  })

  expect_true(all(bool_vec))
})

test_that("declutter spaces the jumps properly", {
  trials <- 100

  bool_vec <- sapply(1:trials, function(x){
    set.seed(10*x)
    jump_vec <- sort(sample(1:100, 15))
    sign_vec <- sample(c(-1,1), length(jump_vec), replace = T)

    res <- declutter(jump_vec, sign_vec, how_close = 5)

    all(diff(res$jump_vec) > 5)
  })

  expect_true(all(bool_vec))
})

test_that("declutter forms the jump centers correctly", {
  trials <- 100
  dis <- 5

  bool_vec <- sapply(1:trials, function(x){
    set.seed(10*x)
    jump_vec <- sort(sample(1:100, 15))
    sign_vec <- sample(c(-1,1), length(jump_vec), replace = T)

    res <- declutter(jump_vec, sign_vec, how_close = dis)

    jump_center <- res$jump_vec
    jump_list <- vector("list", length(jump_center))

    jump_tmp <- jump_vec
    len <- length(jump_center)
    # initialize
    for(y in length(jump_tmp):1) {
      idx <- which(abs(jump_tmp[y] - jump_center) <= dis)
      if(length(idx) >= 1){
        jump_list[[idx]] <- sort(c(jump_list[[idx]], jump_tmp[y]))
        jump_tmp <- jump_tmp[-y]
      }
    }
    # fill in
    while(length(jump_tmp) > 0){
      for(y in length(jump_tmp):1){
        idx <- which(sapply(1:len, function(z){
          any(abs(jump_tmp[y] - jump_list[[z]]) <= dis)
        }))
        if(length(idx) >= 1){
          jump_list[[idx]] <- sort(c(jump_list[[idx]], jump_tmp[y]))
          jump_tmp <- jump_tmp[-y]
        }
      }
    }

    res2 <- floor(sapply(jump_list, median))

    all(res$jump_vec == res2)
  })

  expect_true(all(bool_vec))
})

test_that("declutter groups the signs correctly", {
  trials <- 100
  dis <- 5

  bool_vec <- sapply(1:trials, function(x){
    set.seed(10*x)
    jump_vec <- sort(sample(1:100, 15))
    sign_vec <- sample(c(-1,1), length(jump_vec), replace = T)

    res <- declutter(jump_vec, sign_vec, how_close = dis)

    jump_center <- res$jump_vec
    jump_list <- vector("list", length(jump_center))

    jump_tmp <- jump_vec
    len <- length(jump_center)
    # initialize
    for(y in length(jump_tmp):1) {
      idx <- which(abs(jump_tmp[y] - jump_center) <= dis)
      if(length(idx) >= 1){
        jump_list[[idx]] <- sort(c(jump_list[[idx]], jump_tmp[y]))
        jump_tmp <- jump_tmp[-y]
      }
    }
    # fill in
    while(length(jump_tmp) > 0){
      for(y in length(jump_tmp):1){
        idx <- which(sapply(1:len, function(z){
          any(abs(jump_tmp[y] - jump_list[[z]]) <= dis)
        }))
        if(length(idx) >= 1){
          jump_list[[idx]] <- sort(c(jump_list[[idx]], jump_tmp[y]))
          jump_tmp <- jump_tmp[-y]
        }
      }
    }

    res2 <- sapply(jump_list, function(x){
      idx <- which(jump_vec %in% x)
      tmp <- sign_vec[idx]
      c(length(which(tmp == -1)), length(which(tmp == 1)))
    })

    all(res$sign_mat == res2)
  })

  expect_true(all(bool_vec))
})

test_that("declutter determines the desired_locations correctly", {
  trials <- 100
  dis <- 5

  bool_vec <- sapply(1:trials, function(x){
    set.seed(10*x)
    jump_vec <- sort(sample(1:100, 15))
    sign_vec <- sample(c(-1,1), length(jump_vec), replace = T)

    res <- declutter(jump_vec, sign_vec, how_close = dis)

    jump_center <- res$jump_vec
    jump_list <- vector("list", length(jump_center))

    jump_tmp <- jump_vec
    len <- length(jump_center)
    # initialize
    for(y in length(jump_tmp):1) {
      idx <- which(abs(jump_tmp[y] - jump_center) <= dis)
      if(length(idx) >= 1){
        jump_list[[idx]] <- sort(c(jump_list[[idx]], jump_tmp[y]))
        jump_tmp <- jump_tmp[-y]
      }
    }
    # fill in
    while(length(jump_tmp) > 0){
      for(y in length(jump_tmp):1){
        idx <- which(sapply(1:len, function(z){
          any(abs(jump_tmp[y] - jump_list[[z]]) <= dis)
        }))
        if(length(idx) >= 1){
          jump_list[[idx]] <- sort(c(jump_list[[idx]], jump_tmp[y]))
          jump_tmp <- jump_tmp[-y]
        }
      }
    }

    res2 <- floor(sapply(jump_list, median))

    all(res$jump_vec == res2)
  })

  expect_true(all(bool_vec))
})

######################

## .consecutive_boolean is correct

test_that(".consecutive_boolean works on an example case", {
  res <- .consecutive_boolean(c(T, F, F))

  expect_true(nrow(res) == 3)
  expect_true(ncol(res) == 2)
  expect_true(all(res[,1] <= res[,2]))
})

#######################

## contrast_from_cluster is correct

test_that("contrast_from_cluster works", {
  jump_vec <- c(40,41,42,50,52,60,63,80,120,121,122)
  sign_vec <- rep(1,length(jump_vec))
  sign_vec[2] <- -1

  cluster_list <- declutter(jump_vec, sign_vec)
  res <- contrast_from_cluster(cluster_list, 200, 1)

  expect_true(is.numeric(res))
  expect_true(any(res > 0))
  expect_true(any(res < 0))
  expect_true(length(res) == 200)
})

test_that("contrast_from_cluster is actually right mean minus left mean", {
  set.seed(20)
  n <- 100
  jump_vec <- sort(sample(1:n, 10))
  sign_vec <- sample(c(-1,1), length(jump_vec), replace = T)

  cluster_list <- declutter(jump_vec, sign_vec)
  contrast <- contrast_from_cluster(cluster_list, n, 1)

  y <- rnorm(n)
  res <- as.numeric(contrast %*% y)
  res2 <- mean(y[(cluster_list$jump_vec[1]+1):cluster_list$jump_vec[2]]) -
    mean(y[1:cluster_list$jump_vec[1]])

  expect_true(abs(res - res2) <= 1e-6)
})

test_that("contrast_from_cluster works when 1 is the first jump", {
  set.seed(10)
  n <- 100
  jump_vec <- c(1, 10,11,12)
  sign_vec <- sample(c(-1,1), length(jump_vec), replace = T)

  cluster_list <- declutter(jump_vec, sign_vec)
  location <- 1
  contrast <- contrast_from_cluster(cluster_list, n, location)

  median_jumps <- c(0, cluster_list$jump_vec, n-1)

  y <- rnorm(n)
  res <- as.numeric(contrast %*% y)
  res2 <- mean(y[(median_jumps[location+1]+1):median_jumps[location+2]]) -
    mean(y[(median_jumps[location]+1):median_jumps[location+1]])

  expect_true(abs(res - res2) <= 1e-6)
})

test_that("contrast_from_cluster works when n-1 is the last jump", {
  set.seed(10)
  n <- 100
  jump_vec <- c(10, 11, 12, n-1)
  sign_vec <- sample(c(-1,1), length(jump_vec), replace = T)

  cluster_list <- declutter(jump_vec, sign_vec)
  location <- length(cluster_list$jump_vec)
  contrast <- contrast_from_cluster(cluster_list, n, location)

  median_jumps <- c(0, cluster_list$jump_vec, n)

  y <- rnorm(n)
  res <- as.numeric(contrast %*% y)
  res2 <- mean(y[(median_jumps[location+1]+1):median_jumps[location+2]]) -
    mean(y[(median_jumps[location]+1):median_jumps[location+1]])

  expect_true(abs(res - res2) <= 1e-6)
})

test_that("contrast_from_cluster is actually right mean minus left mean for many cases", {
  trials <- 100
  n <- 100

  bool_vec <- sapply(1:trials, function(x){
    set.seed(10*x)
    jump_vec <- sort(sample(1:(n-1), 10))
    sign_vec <- sample(c(-1,1), length(jump_vec), replace = T)

    cluster_list <- declutter(jump_vec, sign_vec)
    location <- sample(1:length(cluster_list$jump_vec), 1)
    contrast <- contrast_from_cluster(cluster_list, n, location)

    median_jumps <- c(0, cluster_list$jump_vec, n)

    y <- rnorm(n)
    res <- as.numeric(contrast %*% y)
    res2 <- mean(y[(median_jumps[location+1]+1):median_jumps[location+2]]) -
      mean(y[(median_jumps[location]+1):median_jumps[location+1]])

    abs(res - res2) <= 1e-6
  })

  expect_true(all(bool_vec))
})

