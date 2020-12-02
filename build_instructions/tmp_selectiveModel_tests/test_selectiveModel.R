context("Test selective model")

## .segments is correct

test_that(".segments works", {
  res <- .segments(10, c(3, 7))

  expect_true(is.numeric(res))
  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(3, 10)))
})

test_that(".segments forms the proper matrix", {
  res <- .segments(10, c(3, 7))

  expect_true(all(rowSums(res) == 1))
  for(i in 1:3){
    expect_true(length(unique(res[i,])) == 2)
  }
  expect_true(all(res[1,1:3] != 0))
  expect_true(all(res[1,4:10] == 0))
  expect_true(all(res[2,4:7] != 0))
  expect_true(all(res[2,c(1:3,8:10)] == 0))
  expect_true(all(res[3,8:10] != 0))
  expect_true(all(res[3,1:7] == 0))
})

test_that(".segments can ignore properly", {
  res <- .segments(10, c(3, 7), ignore_jump = 1)

  expect_true(all(dim(res) == c(2, 10)))
  expect_true(all(res[1,] == c(rep(1/7,7), rep(0,3))))
})

########################

## .segment_means is correct

test_that(".segment_means works", {
  set.seed(10)
  y <- rnorm(20)
  segments <- .segments(length(y), c(5, 10, 15))
  res <- .segment_means(y, segments)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 4)
})

test_that(".segment_means returns the proper means", {
  set.seed(10)
  y <- rnorm(20)
  segments <- .segments(length(y), c(5, 10, 15))
  res <- .segment_means(y, segments)

  target <- c(mean(y[1:5]), mean(y[6:10]), mean(y[11:15]), mean(y[16:20]))

  expect_true(sum(abs(res - target)) < 1e-6)
})

############################

## selected_model_inference is correct

test_that("selected_model_inference works", {
  set.seed(10)
  y <- c(rep(0, 10), rep(1, 10)) + 0.01*rnorm(20)
  fit_method <- function(y){binseginf::bsfs(y, 1)}
  res <- selected_model_inference(y, fit_method, num_samp = 10, verbose = F,
                                  ignore_jump = 1, param = list(burn_in = 10))

  expect_true(length(res) == 3)
  expect_true(all(names(res) == c("pval", "test_stat", "null_stat")))
  expect_true(is.numeric(res$pval))
  expect_true(!is.matrix(res$pval))
  expect_true(length(res$pval) == 1)
  expect_true(res$pval >= 0)
  expect_true(res$pval <= 1)
})

test_that("selected_model_inference works for known sigma", {
  set.seed(10)
  y <- c(rep(0, 10), rep(1, 10)) + 0.01*rnorm(20)
  fit_method <- function(y){binseginf::bsfs(y, 1)}
  res <- selected_model_inference(y, fit_method, sigma = 1, num_samp = 10, verbose = F,
                                  param = list(burn_in = 10))

  expect_true(length(res) == 3)
  expect_true(all(names(res) == c("pval", "test_stat", "null_stat")))
  expect_true(is.numeric(res$pval))
  expect_true(!is.matrix(res$pval))
  expect_true(length(res$pval) == 1)
  expect_true(res$pval >= 0)
  expect_true(res$pval <= 1)
})

test_that("selected_model_inference works reasonably relatively for hit run", {
  fit_method <- function(x){binseginf::bsfs(x, numSteps = 1)}
  test_func <- selectiveModel::segment_difference
  num_samp <- 500
  cores <- NA

  set.seed(10)
  y1 <- c(rep(0, 3), rep(5, 3)) + rnorm(6)
  res1 <- selected_model_inference(y1, fit_method = fit_method, test_func = test_func,
                                   num_samp = num_samp, ignore_jump = 1,
                                   cores = cores, verbose = F, param = list(burn_in = 2000, lapse = 1))


  set.seed(10)
  y2 <- rep(0, 6) + rnorm(6)
  res2 <- selected_model_inference(y2, fit_method = fit_method, test_func = test_func,
                                   num_samp = num_samp, ignore_jump = 1,
                                   cores = cores, verbose = F, param = list(burn_in = 2000, lapse = 1))

  expect_true(res1$pval < res2$pval)
})

test_that("selected_model_inference works for one problem case", {
  set.seed(77)
  y <- rnorm(20)
  fit_method <- function(y){binseginf::bsfs(y, 1)}

  res <- selected_model_inference(y, fit_method, verbose = F, cores = NA,
                                  num_samp = 50,
                                  param = list(burn_in = 10, time_limit = 600))

  expect_true(length(res) == 3)
})

test_that("selected_model_inference works for one problem case", {
  set.seed(238)
  y <- rnorm(20)
  fit_method <- function(y){binseginf::bsfs(y, 1)}

  res <- selected_model_inference(y, fit_method, verbose = F, cores = NA,
                                  num_samp = 50,
                                  param = list(burn_in = 10, time_limit = 600))

  expect_true(length(res) == 3)
})

test_that("selected_model_inference works for one-sided p-values", {
  set.seed(10)
  y <- c(rnorm(10), rnorm(10)+10)
  fit_method <- function(y){binseginf::bsfs(y, 1)}

  res <- selected_model_inference(y, fit_method, verbose = F, cores = NA,
                                  num_samp = 50, direction = 1,
                                  param = list(burn_in = 10, time_limit = 600))

  expect_true(length(res) == 3)
})

test_that("selected_model_inference gives one-sided p-values that are
          smaller than or equal to two-sided p-values", {
  trials <- 100

  res <- sapply(1:trials, function(x){
    set.seed(10*x)
    y <- rnorm(10)
    y <- y - mean(y)
    fit_method <- function(y){binseginf::bsfs(y, 1)}
    fit <- fit_method(y)

    set.seed(10*x)
    res1 <- selected_model_inference(y, fit_method, verbose = F, cores = NA,
                                    num_samp = 50, direction = sign(fit$cp.sign),
                                    sigma = 1, ignore_jump = 1,
                                    param = list(burn_in = 200, time_limit = 600))

    set.seed(10*x)
    res2 <- selected_model_inference(y, fit_method, verbose = F, cores = NA,
                                    num_samp = 50, direction = NA,
                                    sigma = 1, ignore_jump = 1,
                                    param = list(burn_in = 200, time_limit = 600))

    res1$pval <= res2$pval
  })

  expect_true(all(res))
})

test_that("selected_model_inference still perserves mean when decluttering",{
  set.seed(10)
  middle_mutation <- function(lev, n){
    mn <- rep(0,n)
    mn[seq(from=n/2+1, to=n/2+round(.2*n))] <- lev
    mn
  }
  true_jumps <- c(10, 14)

  test_func_closure <- function(contrast){
    function(y, fit = NA, jump = NA){
      as.numeric(contrast %*% y)
    }
  }
  declutter_func <- function(x){selectiveModel::declutter(x, sign_vec = rep(1, length(x)),
                                                          how_close = 2)$jump_vec}
  num_samp <- 100
  burn_in <- 10
  numSteps <- 4

  n <- 20
  dat <- middle_mutation(lev = 2, n = n) + stats::rnorm(n)
  fit_method <- function(x){binseginf::bsfs(x, numSteps = numSteps)}

  fit <- fit_method(dat)
  sign_mat <- binseginf::jump_sign(fit)
  cluster_list <- selectiveModel::declutter(jump_vec = sign_mat[,1], sign_vec = sign_mat[,2],
                                            how_close = 2,
                                            desired_jumps = true_jumps)

  i <- 2
  set.seed(10)
  contrast <- contrast_from_cluster(cluster_list, n, i)
  test_func <- test_func_closure(contrast)
  if(cluster_list$sign_mat["sign:-1",i] == 0){
    direction <- 1
  } else if(cluster_list$sign_mat["sign:+1",i] == 0){
    direction <- -1
  } else {
    direction <- NA
  }

  tmp <- selectiveModel::selected_model_inference(dat, fit_method = fit_method,
                                                  test_func = test_func,
                                                  declutter_func = declutter_func,
                                                  num_samp = num_samp,
                                                  direction = direction,
                                                  ignore_jump = i, sigma = 1,
                                                  verbose = F, param = list(burn_in = burn_in,
                                                                            lapse = 1),
                                                  return_samples = T)

  null_samples <- tmp$null_samples
  intended_jumps <- c(0,cluster_list$jump_vec[-i],n)
  null_means <- apply(null_samples, 2, function(x){
    sapply(1:(length(intended_jumps)-1), function(k){
      mean(x[(intended_jumps[k]+1):intended_jumps[k+1]])
    })
  })

  target_means <- sapply(1:(length(intended_jumps)-1), function(k){
    mean(dat[(intended_jumps[k]+1):intended_jumps[k+1]])
  })

  expect_true(all(sapply(1:length(target_means), function(x){
    sum(abs(null_means[i,] - target_means[i])) < 1e-6
  })))
})


test_that("selected_model_inference gives one-sided p-values that are
          smaller than or equal to two-sided p-values", {
            trials <- 100

            res <- sapply(1:trials, function(x){
              set.seed(10*x)
              y <- rnorm(10)
              y <- y - mean(y)
              fit_method <- function(y){binseginf::bsfs(y, 1)}
              fit <- fit_method(y)

              set.seed(10*x)
              res1 <- selected_model_inference(y, fit_method, verbose = F, cores = NA,
                                               num_samp = 50, direction = sign(fit$cp.sign),
                                               sigma = 1, ignore_jump = 1,
                                               param = list(burn_in = 200, time_limit = 600))

              set.seed(10*x)
              res2 <- selected_model_inference(y, fit_method, verbose = F, cores = NA,
                                               num_samp = 50, direction = NA,
                                               sigma = 1, ignore_jump = 1,
                                               param = list(burn_in = 200, time_limit = 600))

              res1$pval <= res2$pval
            })

            expect_true(all(res))
            })

test_that("selected_model_inference works for null means not equal to 0, known sigma",{
  set.seed(10)
  y <- c(rep(0, 10), rep(1, 10)) + rnorm(20)
  fit_method <- function(y){binseginf::bsfs(y, 3)}
  fit <- fit_method(y)
  sign_mat <- binseginf::jump_sign(fit)

  cluster_list <- selectiveModel::declutter(jump_vec = sign_mat[,1], sign_vec = sign_mat[,2],
                                            how_close = 2)

  contrast <- contrast_from_cluster(cluster_list, 20, 2)

  val <- contrast %*% c(rep(0, 10), rep(1, 10))
  null_mean <- rep(0, 20)
  null_mean[which(contrast > 0)] <- val

  stopifnot(abs(as.numeric(contrast %*% null_mean) - val) <= 1e-6)

  res <- selected_model_inference(y, fit_method, sigma = 1, num_samp = 10,
                                  verbose = F, null_mean = null_mean,
                                  param = list(burn_in = 10))

  expect_true(length(res) == 3)
  expect_true(all(names(res) == c("pval", "test_stat", "null_stat")))
  expect_true(is.numeric(res$pval))
  expect_true(!is.matrix(res$pval))
  expect_true(length(res$pval) == 1)
  expect_true(res$pval >= 0)
  expect_true(res$pval <= 1)
})

test_that("selected_model_inference works for null means not equal to 0 gives lower pvalue",{
  set.seed(10)
  y <- c(rep(0, 10), rep(1, 10)) + rnorm(20)
  fit_method <- function(y){binseginf::bsfs(y, 3)}
  fit <- fit_method(y)
  sign_mat <- binseginf::jump_sign(fit)

  cluster_list <- selectiveModel::declutter(jump_vec = sign_mat[,1], sign_vec = sign_mat[,2],
                                            how_close = 2)

  contrast <- contrast_from_cluster(cluster_list, 20, 2)

  val <- contrast %*% c(rep(0, 10), rep(1, 10))
  null_mean <- rep(0, 20)
  null_mean[which(contrast > 0)] <- val

  stopifnot(abs(as.numeric(contrast %*% null_mean) - val) <= 1e-6)

  declutter_func <- function(x){selectiveModel::declutter(x, sign_vec = rep(1, length(x)),
                                                          how_close = 2)$jump_vec}
  test_func_closure <- function(contrast){
    function(y, fit = NA, jump = NA){
      as.numeric(contrast %*% y)
    }
  }

  test_func <- test_func_closure(contrast)

  res <- selectiveModel::selected_model_inference(y, fit_method = fit_method,
                                                  test_func = test_func,
                                                  declutter_func = declutter_func,
                                                  null_mean = null_mean,
                                                  num_samp = 1000,
                                                  ignore_jump = 2, sigma = 1,
                                                  verbose = F, param = list(burn_in = 10,
                                                                            lapse = 1))

  res2 <- selectiveModel::selected_model_inference(y, fit_method = fit_method,
                                                  test_func = test_func,
                                                  declutter_func = declutter_func,
                                                  null_mean = rep(0, 20),
                                                  num_samp = 1000,
                                                  ignore_jump = 2, sigma = 1,
                                                  verbose = F, param = list(burn_in = 10,
                                                                            lapse = 1))

  expect_true(res$pval >= res2$pval)
})





