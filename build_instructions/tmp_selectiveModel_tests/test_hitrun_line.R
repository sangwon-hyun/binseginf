context("Test hit and run, Line")

## .remove_nullspace_gaussian is correct

test_that(".remove_nullspace_gaussian works", {
  set.seed(10)
  n <- 10
  gaussian <- .gaussian(rep(0, n), diag(n))
  y <- rnorm(n)
  segments <- .segments(n, 5)

  nullspace_mat <- .sample_matrix_space(segments)
  mean_val <- as.numeric(segments%*%y)
  segments_full <- rbind(t(nullspace_mat), segments)

  res <- .remove_nullspace_gaussian(gaussian, segments_full, mean_val)

  expect_true(class(res) == "gaussian")
  expect_true(length(res$mean) + nrow(segments) == length(gaussian$mean))
  expect_true(all(dim(res$covariance) + nrow(segments) == dim(gaussian$covariance)))
})

#######

## .remove_nullspace_polyhedra is correct

test_that(".remove_nullspace_polyhedra works", {
  set.seed(10)
  y <- c(rep(0, 10), rep(1, 10)) +rnorm(20)
  fit <- binseginf::bsfs(y, 1)
  polyhedra <- binseginf::polyhedra(fit)

  segments <- .segments(20, binseginf::jumps(fit))

  nullspace_mat <- .sample_matrix_space(segments)
  mean_val <- as.numeric(segments%*%y)
  segments_full <- rbind(t(nullspace_mat), segments)

  res <- .remove_nullspace_polyhedra(polyhedra, segments_full, mean_val)

  expect_true(class(res) == "polyhedra")
  expect_true(all(length(res$u) == length(polyhedra$u)))
  expect_true(all(dim(res$gamma) + c(0, nrow(segments)) == dim(polyhedra$gamma)))
})


########

## .remove_nullspace is correct

test_that(".remove_nullspace works", {
  n <- 10
  set.seed(10)
  y <- c(rep(0, n/2), rep(1, n/2)) + rnorm(n)
  fit <- binseginf::bsfs(y, 1)
  polyhedra <- binseginf::polyhedra(fit)
  gaussian <- .gaussian(rep(0, n), diag(n))

  segments <- .segments(n, binseginf::jumps(fit))
  nullspace_mat <- .sample_matrix_space(segments)
  mean_val <- as.numeric(segments%*%y)
  segments_full <- rbind(t(nullspace_mat), segments)

  res <- .remove_nullspace(gaussian, polyhedra, segments_full, mean_val)

  expect_true(is.list(res))
  expect_true("forward_translation" %in% names(res))
  expect_true("backward_translation" %in% names(res))
  expect_true(class(res$gaussian) == "gaussian")
  expect_true(class(res$polyhedra) == "polyhedra")
})

test_that(".remove_nullspace preserves the set of allowed vectors", {
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

  res <- .remove_nullspace(gaussian, polyhedra, segments_full, mean_val)

  trials <- 1000
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    vec <- rnorm(n)
    vec[1:jump] <- vec[1:jump] + (mean_val[1] - mean(vec[1:jump]))
    vec[(jump+1):n] <- vec[(jump+1):n] + (mean_val[2] - mean(vec[(jump+1):n]))

    bool1 <- all(polyhedra$gamma %*% vec >= polyhedra$u)
    bool2 <- all(res$polyhedra$gamma %*% res$forward_translation(vec) >= res$polyhedra$u)

    bool1 == bool2
  })

  expect_true(all(bool_vec))
})

test_that(".remove_nullspace has correct forward and backward translations", {
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

  res <- .remove_nullspace(gaussian, polyhedra, segments_full, mean_val)

  trials <- 1000
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    vec <- rnorm(n)
    vec[1:jump] <- vec[1:jump] + (mean_val[1] - mean(vec[1:jump]))
    vec[(jump+1):n] <- vec[(jump+1):n] + (mean_val[2] - mean(vec[(jump+1):n]))

    vec2 <- res$backward_translation(res$forward_translation(vec))

    sum(abs(vec - vec2)) < 1e-6
  })

  expect_true(all(bool_vec))
})

########

## .whiten_polyhedra is correct

test_that(".whiten_polyhedra works", {
  n <- 10
  set.seed(10)
  y <- c(rep(0, n/2), rep(1, n/2)) + rnorm(n)
  fit <- binseginf::bsfs(y, 1)
  polyhedra <- binseginf::polyhedra(fit)

  sqrt_cov <- diag(n)
  mean_vec <- rep(0, n)

  res <- .whiten_polyhedra(polyhedra, sqrt_cov, mean_vec)

  expect_true(class(res) == "polyhedra")
  expect_true(all(length(res$u) == length(polyhedra$u)))
  expect_true(all(dim(res$gamma) == dim(polyhedra$gamma)))
})

########

## .whiten is correct

test_that(".whiten works", {
  n <- 10
  set.seed(10)
  y <- c(rep(0, n/2), rep(1, n/2)) + rnorm(n)
  fit <- binseginf::bsfs(y, 1)
  polyhedra <- binseginf::polyhedra(fit)
  gaussian <- .gaussian(rep(0, n), diag(n))

  res <- .whiten(gaussian, polyhedra)

  expect_true(is.list(res))
  expect_true("forward_translation" %in% names(res))
  expect_true("backward_translation" %in% names(res))
  expect_true(class(res$gaussian) == "gaussian")
  expect_true(class(res$polyhedra) == "polyhedra")
})

test_that(".whiten preserves the polyhedra", {
  n <- 10
  set.seed(10)
  y <- rnorm(n)
  fit <- binseginf::bsfs(y, 1)
  polyhedra <- binseginf::polyhedra(fit)
  cov_mat <- matrix(1, n, n)
  diag(cov_mat) <- 2
  gaussian <- .gaussian(rep(0, n), cov_mat)

  res <- .whiten(gaussian, polyhedra)

  trials <- 1000

  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    vec <- rnorm(n)

    bool1 <- all(polyhedra$gamma %*% vec >= polyhedra$u)

    vec2 <- res$forward_translation(vec)
    bool2 <- all(res$polyhedra$gamma %*% vec2 >= res$polyhedra$u)

    bool1 == bool2
  })

  expect_true(all(bool_vec))
})

test_that(".remove_nullspace has correct forward and backward translations", {
  n <- 10
  set.seed(10)
  y <- rnorm(n)
  fit <- binseginf::bsfs(y, 1)
  polyhedra <- binseginf::polyhedra(fit)
  cov_mat <- matrix(1, n, n)
  diag(cov_mat) <- 2
  gaussian <- .gaussian(rep(0, n), cov_mat)

  res <- .whiten(gaussian, polyhedra)

  trials <- 1000
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    vec <- rnorm(n)
    vec2 <- res$backward_translation(res$forward_translation(vec))

    sum(abs(vec - vec2)) < 1e-6
  })

  expect_true(all(bool_vec))
})

############

## .generate_directions is correct

test_that(".generate_directions works", {
  res <- .generate_directions(10)

  expect_true(all(dim(res) == c(20,10)))
})

############

## .sampler_hit_run_line is correct

test_that(".sampler_hit_run_line works", {
  set.seed(10)
  n <- 10
  y <- c(rep(0, n/2), rep(1, n/2)) + rnorm(n)
  fit <- binseginf::bsfs(y, 1)
  polyhedra <- binseginf::polyhedra(fit)
  gaussian <- .gaussian(rep(0, n), diag(n))

  segments <- .segments(n, binseginf::jumps(fit), 1)

  res <- .sampler_hit_run_line(y, gaussian, segments, polyhedra, num_samp = 100)

  expect_true(all(dim(res) == c(10, 100)))
})

test_that(".sampler_hit_run_line gives the correct mean for no jumps", {
  set.seed(20)
  n <- 10
  y <- c(rep(0, n/2), rep(1, n/2)) + rnorm(n)
  fit <- binseginf::bsfs(y, 1)
  polyhedra <- binseginf::polyhedra(fit)
  gaussian <- .gaussian(rep(0, n), diag(n))

  segments <- .segments(n, binseginf::jumps(fit), 1)

  res <- .sampler_hit_run_line(y, gaussian, segments, polyhedra, num_samp = 100)

  bool_vec <- apply(res, 2, function(x){
    abs(mean(x) - mean(y)) < 1e-6
  })

  expect_true(all(bool_vec))
})

test_that(".sampler_hit_run_line gives points all in polyhedra", {
  set.seed(30)
  n <- 10
  y <- c(rep(0, n/2), rep(1, n/2)) + rnorm(n)
  fit <- binseginf::bsfs(y, 1)
  polyhedra <- binseginf::polyhedra(fit)
  gaussian <- .gaussian(rep(0, n), diag(n))

  segments <- .segments(n, binseginf::jumps(fit), 1)

  res <- .sampler_hit_run_line(y, gaussian, segments, polyhedra, num_samp = 100)

  bool_vec <- apply(res, 2, function(x){
    all(polyhedra$gamma %*% x >= polyhedra$u)
  })

  expect_true(all(bool_vec))
})
