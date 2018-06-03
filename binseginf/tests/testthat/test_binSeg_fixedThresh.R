context("Test binary segmentation with fixed threshold")

## Data settings
n = 20
threshold = 2
lev = 0
sigma = 1

test_that("Algorithm output is the same as wbs:sbs() each time", {

    for(seed in 1:1000){

        ## Generate some data
        mn <- rep(c(0,lev), each=n/2)
        y0 <- mn + rnorm(n, 0, sigma)

        ## Generate some random threshold
        thresh = runif(1,0,sigma)

        ## From wbs() from wbs package
        w <- wbs::sbs(y0,th=thresh)
        w.cpt <- wbs::changepoints(w,th=thresh)
        c1 = unlist(w.cpt$cpt.th)

        ## From our function
        b = binSeg_fixedThresh(y0, thresh, verbose=FALSE, return.env=FALSE)
        c2 = b$cp

        if(all(is.na(c1))){
            expect_equal(c2, integer())
        } else {
            expect_equal(sort(c1),sort(c2))
        }
    }
})



context("Test binary segmentation with fixed threshold")

## Data settings
numIntervals = 10
n = 20
threshold = 2
lev = 0
sigma = 1

test_that("bsFt polyhedra is correct.", {
  
  ## Generate data, run algorithm
  ## set.seed(0) gives single changepoint at -5
  mn <- rep(c(0,lev), each=n/2)
  y0 <- mn + rnorm(n, 0, sigma)
  thresh = 1
  obj = binSeg_fixedThresh(y0, thresh, verbose=FALSE, return.env=FALSE)
  
  ## Check that
  mypoly = polyhedra(obj, verbose=TRUE)
  nsim = 1000
  for(isim in 1:nsim){
    ynew = y0 + rnorm(n,0,0.1)
    objnew = binSeg_fixedThresh(ynew, thresh, verbose=FALSE,
                                return.env=FALSE)
    if(all(mypoly$gamma %*% cbind(ynew) >= mypoly$u)){
      expect_equal(obj$cp*obj$cp.sign, objnew$cp*obj$cp.sign)
    } else {
      expect_true(all.equal(obj$cp*obj$cp.sign, objnew$cp*obj$cp.sign)==FALSE | length(objnew$cp)==0)
    }
  }
  ## obj$cp*obj$cp.sign
  ## objnew$cp*obj$cp.sign
})





