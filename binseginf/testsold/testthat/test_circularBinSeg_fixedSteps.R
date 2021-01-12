context("Test circular binary segmentation fixed steps")

## .cusum_cbs is correct

test_that(".cusum_cbs works", {
  vec <- cumsum(1:10)
  x <- c(3,5)
  res <- .cusum_cbs(x, vec)

  answer <- sqrt(1/(1/3 + 1/7)) * (mean(3:5) - mean(c(1:2,6:10)))
  expect_true(abs(res - answer) <= 1e-6)
})

test_that(".cusum_cbs works when the left side is length 1", {
  vec <- cumsum(1:10)
  x <- c(2,5)
  res <- .cusum_cbs(x, vec)

  answer <- sqrt(1/(1/4 + 1/6)) * (mean(2:5) - mean(c(1,6:10)))
  expect_true(abs(res - answer) <= 1e-6)
})

test_that(".cusum_cbs works when the right side is length 1", {
  vec <- cumsum(1:10)
  x <- c(6,9)
  res <- .cusum_cbs(x, vec)

  answer <- sqrt(1/(1/4 + 1/6)) * (mean(6:9) - mean(c(1:5,10)))
  expect_true(abs(res - answer) <= 1e-6)
})

test_that(".cusum_cbs works when the right side is length 1", {
  vec <- cumsum(1:10)
  x <- c(6,9)
  res <- .cusum_cbs(x, vec)

  answer <- sqrt(1/(1/4 + 1/6)) * (mean(6:9) - mean(c(1:5,10)))
  expect_true(abs(res - answer) <= 1e-6)
})

test_that(".cusum_cbs works when there is no left side", {
  vec <- cumsum(1:10)
  x <- c(1,5)
  res <- .cusum_cbs(x, vec)

  answer <- sqrt(1/(1/5 + 1/5)) * (mean(1:5) - mean(c(6:10)))
  expect_true(abs(res - answer) <= 1e-6)
})

test_that(".cusum_cbs works when there is no right side", {
  vec <- cumsum(1:10)
  x <- c(7,10)
  res <- .cusum_cbs(x, vec)

  answer <- sqrt(1/(1/4 + 1/6)) * (mean(7:10) - mean(c(1:6)))
  expect_true(abs(res - answer) <= 1e-6)
})

test_that(".cusum_cbs fails when the bounds are the entire interval", {
  vec <- cumsum(1:10)
  x <- c(1,10)
  expect_error(.cusum_cbs(x, vec))
})

test_that(".cusum_cbs computes correctly when there are negative values", {
  set.seed(10)
  y <- c(rnorm(10), rnorm(10, mean = 10), rnorm(10))
  vec <- cumsum(y)

  x <- c(11, 11)
  res <- .cusum_cbs(x, vec)
  answer <- sqrt(1/(1/1 + 1/29))*mean(y[11] - y[c(1:10, 12:30)])
  expect_true(abs(res - answer) <= 1e-6)

  x <- c(11, 20)
  res <- .cusum_cbs(x, vec)
  answer <- sqrt(1/(1/10 + 1/20))*mean(y[11:20] - y[c(1:10, 21:30)])
  expect_true(abs(res - answer) <= 1e-6)
})

######################

## .enumerate_breakpoints_cbs is correct

test_that(".enumerate_breakpoints_cbs works", {
  n <- 5
  res <- .enumerate_breakpoints_cbs(n)
  expect_true(is.matrix(res))
  expect_true(ncol(res) == 2)
  expect_true(all(res[,1] <= res[,2]))

  map <- res[,1]*6+res[,2]
  expect_true(length(map) == length(unique(map))) #all unique rows
})

test_that(".enumerate_breakpoints_cbs outputs the right rows", {
  n <- 4
  res <- .enumerate_breakpoints_cbs(n)
  answer <- matrix(c(1,1, 1,2, 1,3, 2,2, 2,3, 2,4, 3,3, 3,4, 4,4), ncol = 2, byrow = T)

  expect_true(nrow(res) == nrow(answer))

  map_res <- res[,1]*5+res[,2]
  map_answer <- answer[,1]*5+answer[,2]
  expect_true(all(sort(map_res) == sort(map_answer)))
})

test_that(".enumerate_breakpoints_cbs outputs the right rows", {
  n <- 10
  res <- .enumerate_breakpoints_cbs(n)
  answer <- t(cbind(combn(1:(n-1), 2), rbind(2:(n-1), rep(n, n-2)), rbind(1:n, 1:n)))

  expect_true(nrow(res) == nrow(answer))

  map_res <- res[,1]*11+res[,2]
  map_answer <- answer[,1]*11+answer[,2]
  expect_true(all(sort(map_res) == sort(map_answer)))
})

test_that(".enumerate_breakpoints_cbs works on singleton", {
  res <- .enumerate_breakpoints_cbs(1)

  expect_true(is.matrix(res))
  expect_true(all(dim(res) == c(1,2)))
  expect_true(all(as.numeric(res) == c(1,1)))
})

####################

## .find_breakpoint_cbs is correct

test_that(".find_breakpoint_cbs works", {
  y <- 1:10
  res <- .find_breakpoint_cbs(y)

  expect_true(is.list(res))
  expect_true(length(res) == 2)
  expect_true(length(res$breakpoint) == 2)
  expect_true(res$breakpoint[1] >= 1)
  expect_true(res$breakpoint[2] <= 10)
  expect_true(res$breakpoint[1] <= res$breakpoint[2])
  expect_true(is.numeric(res$cusum))
})

test_that(".find_breakpoint_cbs finds the right breakpoint", {
  set.seed(10)
  y <- c(rnorm(10), rnorm(10, mean = 10), rnorm(10))
  res <- .find_breakpoint_cbs(y)

  expect_true(all(res$breakpoint == c(11,20)))
})

test_that(".find_breakpoint_cbs returns NA when y is length 1", {
  y <- 1
  res <- .find_breakpoint_cbs(y)
  expect_true(is.list(res))
  expect_true(is.na(res$breakpoint))
  expect_true(res$cusum == 0)
})

####################

## circularBinSeg_fixedSteps is correct

test_that("circularBinSeg_fixedSteps works", {
  set.seed(10)
  y <- c(rnorm(10), rnorm(10, mean = 10), rnorm(10))
  res <- circularBinSeg_fixedSteps(y, 1)

  expect_true(class(res) == "cbsFs")
})

test_that("circularBinSeg_fixedSteps works with two jumps", {
  set.seed(10)
  y <- c(rnorm(5), rnorm(5, mean = 10), rnorm(5), rnorm(5, mean = 10), rnorm(5))
  res <- circularBinSeg_fixedSteps(y, 2)

  expect_true(class(res) == "cbsFs")
})

test_that("circularBinSeg_fixedSteps works when jump is at the edge", {
  set.seed(10)
  y <- c(rnorm(10), rnorm(10, mean = 10))
  res <- circularBinSeg_fixedSteps(y, 1)

  expect_true(class(res) == "cbsFs")
})

test_that("circularBinSeg_fixedSteps will work in a strange circumstance where start=end", {
  set.seed(70)
  y <- c(rep(0,10), rep(-2,10), rep(-1,5)) + rnorm(25)
  obj <- circularBinSeg_fixedSteps(y,2)

  expect_true(class(obj) == "cbsFs")
})

test_that("circularBinSeg_fixedSteps works with repetitions", {
  set.seed(10)
  y <- c(rep(0,10), rep(5,10), rep(0,10))
  res <- circularBinSeg_fixedSteps(y, 1)

  expect_true(all(jumps(res, sorted = T) == c(10,20)))
})

####################

## jumps.cbsFs is correct

test_that("jumps.cbsFs works", {
  set.seed(10)
  y <- c(rnorm(5), rnorm(5, mean = 10), rnorm(5), rnorm(5, mean = 10), rnorm(5))
  obj <- circularBinSeg_fixedSteps(y, 2)
  res <- jumps(obj)

  expect_true(all(res == c(5,10,15,20)))
})

test_that("jumps.cbsFs works when one of the jumps is at the edge", {
  set.seed(10)
  y <- c(rnorm(10), rnorm(10, mean = 10))
  obj <- circularBinSeg_fixedSteps(y, 1)
  res <- jumps(obj)

  expect_true(res == 10)
})

test_that("jumps.cbsFs will not report the same jump if sorted is TRUE", {
  set.seed(10)
  y <- c(rnorm(10), rnorm(5, mean = 25), rnorm(5, mean = 20), rnorm(10))
  obj <- circularBinSeg_fixedSteps(y, 2)
  res <- jumps(obj, sorted = T)

  expect_true(all(res == c(10,15,20)))

  res <- jumps(obj, sorted = F)

  expect_true(all(res == c(10, 20, 10, 15)))
})

test_that("jumps.cbsFs works with sorted is FALSE and there are not two shoulders", {
  set.seed(10)
  y <- c(rep(0, 10), rep(5, 10))
  obj <- circularBinSeg_fixedSteps(y, 1)
  res <- jumps(obj, sorted = F)

  expect_true(all(is.na(res[1]), res[2] == 10) | all(res == c(10, 20)))
})


####################
## Test inference for CBS

## Synopsis: trying out if cbs inference works well
test_that("cbsFs saturated tests after two steps gives uniform null p-values", {
    nsim=3000
    onejump <- function(lev,n){c(rep(0,n/2),rep(lev,n/2))}
    results = mclapply(1:nsim, function(isim){
        printprogress(isim,nsim)
        meanfun = onejump

        set.seed(isim)
        ## Generate some data
        n=60
        lev=0
        meanfun=onejump
        mn = meanfun(lev,n)
        sigma=1
        y = mn + rnorm(n, 0, sigma)

        ## Fit CBS
        numSteps=2
        g = circularBinSeg_fixedSteps(y,  numSteps=numSteps)
        poly = polyhedra(g)
        vlist <- make_all_segment_contrasts_from_cp(cp=g$cp, cp.sign = g$cp.sign, n=n)

        ## Get the p-values
        pvs = sapply(vlist, function(v){
               return(poly.pval2(y=y, poly=poly, v=v, sigma=sigma)$pv)
        })
        names(pvs) = (g$cp*g$cp.sign)
        return(data.frame(pvs=pvs,cp=g$cp*g$cp.sign))
    },mc.cores=3)
    mytable = do.call(rbind,results)
    ## qqunif(mytable[,"pvs"])
    expect_equal(ks.test(mytable[,"pvs"], punif)$p.value>0.05, TRUE)
}
