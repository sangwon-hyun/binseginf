context("Test wildBinSeg.R and some helper functions.")

test_that("Intervals are collected correctly", {

    ## Test settings
    n=100
    numIntervals = 100

    ## See if it always gives the right number of elements.
    expect_equal(length(generate_intervals(n, numIntervals)$se),
                 numIntervals)

    ## See if it produces NULL elements.
    intervals = generate_intervals(n, numIntervals)
    expect_false(any(unlist(lapply(intervals, is.null))))
})


test_that("get_which_qualify() is correctly functioning", {

    n=60
    intervals = generate_intervals(n,10, seed=0)
    m = which(.get_which_qualify(15,45,intervals))
    expect_equal(sort(m), c(2,3,4,10))
})


test_that("make_se() is correctly functioning",{

    ## Test setting.
    s=1
    e=60
    y = rnorm(rep(0,4,each=30),0,1)
    intervals = generate_intervals(length(y),10,seed=0)
    thresh=0.5

    m = which(.get_which_qualify(s,e,intervals))
    semat = make_semat(m=m,s=s,e=e,intervals=intervals,y=y,thresh=thresh)

    ## Separately obtain all cusums.
    all.s = lapply(intervals$se[m], function(se) se[1])
    all.e = lapply(intervals$se[m], function(se) se[2])
    all.max.cusums.manual = unlist(Map(function(s,b,e,y)cusum(s=s,e=e,b=b,y=y),
                                       all.s, semat[,"b"], all.e,
                                       rep(list(y),length(m))))

    expect_equal(semat[,"maxcusum"], all.max.cusums.manual)
    ## See if the maximizing is correctly done, with respect to internal and
    ## externally created max cusums.
    expect_equal(as.numeric(abs(semat[which(semat[,"maxhere"]==1), "maxcusum"])),
                 max(abs(semat[,"maxcusum"])))
    expect_equal(abs(all.max.cusums.manual[which(semat[,"maxhere"]==1)]),
                 max(abs(semat[,"maxcusum"])))
})


test_that("wildBinSeg_fixedThresh() doesn't produce environment |env| whose tree element env$tree has null elements",{


    ## Test setting
    s=1
    e=60
    y = rnorm(rep(0,4,each=30),0,1)
    augment=FALSE
    env = wildBinSeg_fixedThresh(y,1,10,return.env=TRUE, augment)

    ## See if the tree has any empty (NULL) elements
    expect_true(all(lapply(env$tree, length)>0))

})


test_that("wildBinSeg_fixedThresh() gives the correct changepoint in a strong-signal case.", {

    ## Test setting
    thresh=10
    seed=0
    set.seed(seed)
    y = rep(c(0,10),each=30) + rnorm(60,0,1)
    numInterval=100

    ## Run WBS
    output = wildBinSeg_fixedThresh(y, thresh, numInterval,seed=seed,augment=FALSE)

    ## See if single upward changepoint is detected
    expect_equal(as.numeric(output$cp), 30)
    expect_equal(as.numeric(output$cp.sign), +1)
})


test_that("wildBinSeg_fixedThresh() handles the silly case of numIntervals=NULL with grace. ", {
    ## Test setting
    thresh=10
    seed=0
    set.seed(seed)
    y = rep(c(0,10),each=30) + rnorm(60,0,1)
    numInterval=0

    ## Run WBS
    output = wildBinSeg_fixedThresh(y, thresh, numInterval,seed=seed,augment=FALSE)

    ## See if single upward changepoint is detected
    expect_equal(as.numeric(output$cp), numeric(0))
    expect_equal(as.numeric(output$cp.sign), numeric(0))
  })

test_that("wildBinSeg_fixedThresh() with augment option and zero intervals has the same output as binSeg_fixedThresh().", {

    ## Test setting
    thresh=10
    seed=0
    set.seed(seed)
    y = rep(c(0,10),each=30) + rnorm(60,0,1)
    numInterval=0

    ## Run WBS
    output = wildBinSeg_fixedThresh(y, thresh, numInterval,seed=seed,augment=TRUE)
    output2 = binSeg_fixedThresh(y, thresh)

    ## See if single upward changepoint is detected
    expect_equal(as.numeric(output$cp), output2$cp)
    expect_equal(as.numeric(output$cp * output$cp.sign),
                 as.numeric(output2$cp * output2$cp.sign))

  })



test_that("make_semat() finds breakpoints that are between the start and end points of the appropriate interval",
{
    ## Test setting
    thresh=10
    n=60
    s=1
    n=60
    seed=0
    set.seed(seed)
    y = rep(c(0,10),each=30) + rnorm(60,0,1)
    numInterval=10
    intervals = generate_intervals(n,numInterval)

    ## Pick some indices and check breakpoints
    set.seed(0)
    M = sample(1:numInterval,10,replace=FALSE)
    semat = make_semat(m=M,s=s,e=e,intervals=intervals,y=y,thresh=thresh)

    ## Check that the breakpoints are all between start/ends
    for(m in M){
        se = intervals$se[[m]]
        expect_true(se[1] <= semat[semat[,"m"]==m,"b"])
        expect_true(semat[semat[,"m"]==m,"b"] <= se[2])
    }
})



test_that("combine.polyhedra() combines simple polyhedra as expected.", {
    polyobj1 = polyobj2 = polyhedra.matrix(matrix(rep(0,12),nrow=2),c(0,0))
    newpolyobj = combine.polyhedra(polyobj1, polyobj2)
    expect_equal(newpolyobj$gamma, matrix(rep(0,24),nrow=4))
    expect_equal(newpolyobj$u, rep(0,4))
})


test_that("Single polyhedron is correct",{

    ## Test settings 1
    n = 6
    lev = 1
    sigma = 1
    thresh = 1
    numInterval = 10
    mn <- rep(c(0,lev), each=n/2)
    set.seed(0)
    y <- mn + rnorm(n, 0, sigma)

    ## Run method, collect things
    obj = wildBinSeg_fixedThresh(y,1,numInterval=numInterval)

    ## Collect a polyhedron
    poly <- polyhedra(obj)

    ## Check if it is correct.
    expect_true(all(poly$gamma %*% y >= poly$u))
})



test_that("Fixed-threshold WBS Polyhedron is exactly correct. (with and without augmentation)",{

    inds1 = c()
    inds2 = c()

    for(lev in c(1,10,100)){
        for(augment in c(FALSE,TRUE)){
            ## Make this a test:
            numIntervals = 1
            n = 10 ## 4
            sigma = 1
            mn <- rep(c(0,lev), each=n/2)
            seed = 48
            set.seed(seed)
            thresh = 0
            y0 <- mn + rnorm(n, 0, sigma)

            ## Run method on original data |y0|, collect things.
            intervals = generate_intervals(n,numIntervals,seed)
            obj = wildBinSeg_fixedThresh(y0, thresh, intervals=intervals, verbose=FALSE,augment=augment)
            poly <- polyhedra.wbsFt(obj)

            ## Generate many new datasets from your polyhedron, see if they /all/ give
            ## you the same fit. No need to do Gaussian generation of ynew.
            ## for(jj in 100000:1001){
            for(jj in 100000:99001){
                ynew = y0 + rnorm(n,0,0.5)
                if(all(poly$gamma%*% (ynew) >= poly$u)){
                    if(!augment)inds1<-c(inds1,jj)
                    if(augment) inds2<-c(inds2,jj)
                    objnew =  wildBinSeg_fixedThresh(y=ynew, thresh=thresh, intervals=intervals,augment=augment)
                    expect_true(all((objnew$cp * objnew$cp.sign) %in% (obj$cp * obj$cp.sign)))
                }
            }
        }
    }
})

test_that("get_vup_vlo() produces numerator and denominator consistent with external p-value functions", {

    ## Sample settings
    lev=3
    thresh=2
    numIntervals=10
    n=10
    sigma=1

    ## Generate data
    mn <- rep(c(0,lev), each=n/2)
    set.seed(2)
    y <- mn + rnorm(n, 0, sigma)

    ## Run method /once/, collect things
    set.seed(3)
    obj = wildBinSeg_fixedThresh(y,thresh, numIntervals)
    if(length(obj$cp)==0)return(NULL)

    ## Collect a polyhedron
    poly <- polyhedra(obj)

    ## Make one segment contrast
    v = make_all_segment_contrasts(obj)[[1]]

    ## Partition the TG statistic
    a = partition_TG(y, poly=poly, v=v, sigma=sigma, nullcontrast=0, bits=50)

    ## Separately make the p-value
    pv1 = poly.pval(y=y,v=v, G=poly$gamma, u=poly$u, sigma=sigma)$pv

    ## Check equality
    expect_equal(a$numer/a$denom, pv1)
    ## expect_equal(a$numer/a$denom, pv2)

})



## TODO: add back after reduction is written into wbsFt functionality.
## test_that("Reduced WBS p-values still returns the same p-values.",{

##     ## Simulation settings
##     set.seed(0)
##     mn <- function(lev,n){c(rep(0,n/2),rep(lev,n/2))}
##     lev = 0
##     sigma = 1
##     n = 20
##     y <- mn(lev,n) + rnorm(n,0,sigma)
##     nsim.is = 1
##     numSteps = 2
##     numIntervals = 30
##     n.levs = 1

##     ## Fit WBS-FS and form contrast
##     method <- wildBinSeg_fixedSteps
##     set.seed(1)
##     intervals <- generate_intervals(n=length(y),numIntervals=numIntervals)
##     obj <- method(y, numSteps=numSteps, intervals=intervals)
##     v <- make_all_segment_contrasts(obj)[[1]]

##     ## New test
##     poly1 <- polyhedra(obj, reduce=TRUE,v=v,sigma=sigma)
##     poly2 <- polyhedra(obj, reduce=FALSE)

##     expect_equal(poly.pval2(y=y,vup=poly1$vup, vlo=poly1$vlo, sigma=sigma, v=v),
##                  poly.pval(y=y, G=poly2$gamma, u=poly2$u, v=v, sigma=sigma, bits=100))
## })


