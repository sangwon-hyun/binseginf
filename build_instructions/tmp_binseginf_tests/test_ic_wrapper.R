context("Test IC wrapper functions")

test_that("IC minimization combined with fixed (nonrandomized) binseg model selection event has uniform p-values.",
{
    ## Settings
    ## lev = 0
    n = 20
    consec = 2
    sigma = 1
    numIntervals = n
    onejump = function(lev,n){c(rep(0,n/2),rep(lev,n/2))}
    meanfun = onejump
    nsim = 5000
    dosim <- function(lev){
        results = lapply(1:nsim, function(isim){

            ## Generate Data
            mn = meanfun(lev,n)
            y = mn + rnorm(n, 0, sigma)

            ## Conduct nonrandomized polyhedron
            h = binSeg_fixedSteps(y, numSteps = 8)
            ic_obj = get_ic(h$cp, h$y, consec=consec, sigma=sigma, type="bic")
            if(ic_obj$flag!="normal") return(NULL)
            stoptime = ic_obj$stoptime
            ic.poly = ic_obj$poly
            poly = polyhedra(h, numSteps=stoptime+consec)

            ## Get combined polyhedron
            cp = h$cp[1:stoptime]
            cp.sign = h$cp.sign[1:stoptime]
            poly$gamma = rbind(poly$gamma, ic.poly$gamma)
            poly$u = c(poly$u, ic.poly$u)

            vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n)
            retain = which(sapply(vlist, function(v){all.equal(sum(v*mn),0)==TRUE}))
            if(length(retain)==0) return(NULL)
            vlist = vlist[retain]

            pvs = lapply(vlist, function(v){
                poly.pval2(y=y,v=v,poly=poly, sigma=sigma)$pv
            })
            return(pvs)
        })
        return(results)
    }
    lev0.results = unlist(dosim(0))
    lev1.results = unlist(dosim(1))

    ## Expect uniform
    expect_equal(ks.test(lev1.results,"punif")$p.value > 0.05, TRUE)
    expect_equal(ks.test(lev0.results,"punif")$p.value > 0.05, TRUE)
})


test_that("IC polyhedron gives exact up and downs.", {

    ## Draw original data.
    n = 10
    mn = rep(0,n)
    sigma = 1
    sigma.add = 0.2
    nsim = 500
    consec=2
    max.numSteps = n-2
    set.seed(11)
    y.orig = mn + rnorm(n,0,sigma)

    ## Nonrandomized inference
    h = binSeg_fixedSteps(y.orig, numSteps = max.numSteps)
    ic_obj = get_ic(h$cp, h$y, consec=consec, sigma=sigma, type="bic")
    if(ic_obj$flag!="normal") return(ic_obj$flag)
    stoptime = ic_obj$stoptime
    orig.seqdirs = c(.getorder(ic_obj$ic))
    ic.poly = ic_obj$poly

    ## Basic check: is my original y even covered?
    poly = polyhedra(h, numSteps=stoptime+consec)
    poly$gamma = rbind(poly$gamma, ic.poly$gamma)
    poly$u = c(poly$u, ic.poly$u)
    expect_true(all(poly$gamma%*%y.orig>poly$u))

    nsim = 200
    results = mclapply(1:nsim, function(isim){

        ## Generate new data
        printprogress(isim,nsim)
        ## set.seed(isim)
        ## y.new = mn + rnorm(n,0,sigma)
        y.new = y.orig + rnorm(n,0,0.1)
        h.new = binSeg_fixedSteps(y.new, numSteps = stoptime+consec)

        ## Get the IC sequence
        ic_obj.new = get_ic(h.new$cp, h.new$y, consec=consec, sigma=sigma, type="bic")
        new.seqdirs = c(.getorder(ic_obj.new$ic))

        ## See if the signs matching make the data be in the polyhedron
        if(all(h.new$cp*h.new$cp.sign  == (h$cp*h$cp.sign)[1:(stoptime+consec)] )){
            if((all.equal(new.seqdirs, orig.seqdirs[1:(stoptime+consec+1)])==TRUE)) {
                expect_true(all(poly$gamma %*%y.new > poly$u))
            }}
    },mc.cores=4)


    nsim=200
    results = mclapply(1:nsim, function(isim){

        ## Generate new data
        printprogress(isim,nsim)
        ## set.seed(isim)
        ## y.new = mn + rnorm(n,0,sigma)
        y.new = y.orig + rnorm(n,0,0.1)
        h.new = binSeg_fixedSteps(y.new, numSteps = stoptime+consec)

        ## Get the IC sequence
        ic_obj.new = get_ic(h.new$cp, h.new$y, consec=consec, sigma=sigma, type="bic")
        new.seqdirs = c(.getorder(ic_obj.new$ic))

        ## See if the signs matching make the data be in the polyhedron
        if(all(poly$gamma %*%y.new > poly$u)){
            expect_equal(ic_obj.new$stoptime, ic_obj$stoptime)
            expect_equal(ic_obj.new$seqdirs, ic_obj$seqdirs[1:(stoptime+consec+1)])
        }
    },mc.cores=4)
}


test_that("IC minimization combined with WBS randomization has uniform p-values.",
{

    ## Source in helper
    source("../main/artificial/artif-helpers.R")

    ## Check uniformity
    n = 10
    mn = rep(0,n)
    sigma = 1
    sigma.add = 0.2
    nsim = 1000
    results = mclapply(1:nsim, function(isim){
        printprogress(isim,nsim)
        y = mn + rnorm(n,0,sigma)
        ## Randomized
        ## pvs = do_rbs_inference(y=y, max.numSteps=8,
        ##                        consec=2, sigma=sigma, postprocess=TRUE,
        ##                        locs=1:length(y), numIS=100, sigma.add = sigma.add, bits=1000,
        ##                        inference.type="rows")

        expect_equal(coords, declutter(coords=coords, how.close=1))

        pvs = do_rwbs_inference(y=y, max.numSteps=10, numIntervals=length(y),
                                consec=2, sigma=sigma, postprocess=TRUE,
                                better.segment=FALSE, locs=1:length(y),
                                numIS=100, inference.type="pre-multiply",
                                improve.nomass.problem=TRUE)
        return(pvs)
    }, mc.cores=4)

    ## Check uniformity
    res = results[sapply(results, function(myresult){length(myresult)>1})]
    all.pvs=unlist(sapply(res,function(a)a$pv))
    qqunif(all.pvs)
    expect_equal(ks.test(all.pvs,"punif")$p.value > 0.05, TRUE)
})


test_that("Helper functions for get_ic() are working correctly", {

    ## Check tDB forming funciton
    n = 15
    cps = c(5,10)
    tDb = make.tDb(cps=cps,n=n)
    expect_equal(which(tDb[,1]!=0), c(1:5))
    expect_equal(which(tDb[,2]!=0), c(6:10))
    expect_equal(which(tDb[,3]!=0), c(11:15))

    ## Check basis forming function
    n = 15
    cps_old = c(5,10)
    cp_new = 7
    d = get_basis(cps_old, cp_new, n)
    expect_true(all(d[6:7] < 0))
    expect_true(all.equal(d[7], d[6]))
    expect_true(all(d[8:10] > 0))
    expect_true(all(d[8]==c(d[9:10])))

    cps_old=10
    cp_new=15
    n=20
    d = get_basis(cps_old, cp_new, n)
    expect_true(all(d[16:20] == 0.2))
    expect_true(all(d[11:15] == -0.2))
    expect_true(all(d[1:10] == 0))
})




test_that("IC minimization combined with additive noise randomization has uniform p-values.",
{
    ## Warning! Even on 8 cores, this takes about 5-10 minutes to run
    n = 30
    lev = 0
    mn = c(rep(0,n/2), rep(lev,n/2))
    sigma = 1
    sigma.add = 0.2
    nsim = 1000
    max.numSteps = 8
    locs = 1:n
    inference.type = "pre-multiply"
    consec=2
    results = mclapply(1:nsim, function(isim){
        printprogress(isim,nsim)
        set.seed(isim)
        y = mn + rnorm(n,0,sigma)
        cumsum.y = cumsum(y)

        ## Fit model and get IC information
        new.noise = rnorm(length(y),0,sigma.add)
        h.fudged = binSeg_fixedSteps(y + new.noise, numSteps=max.numSteps)
        ic_obj = get_ic(h.fudged$cp, h.fudged$y, consec=consec, sigma=sigma+sigma.add, type="bic")
        stoptime = ic_obj$stoptime
        if(ic_obj$flag!="normal"){return(ic_obj$flag)}

        ## Collect stopped model and postprocess
        cp = h.fudged$cp[1:stoptime]
        cp.sign = h.fudged$cp.sign[1:stoptime]
        vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=length(y))

        ## Retain only the changepoints we want results from:
        retain = which((cp %in% locs))
        if(length(retain)==0) return(list(pvs=c(), null.true=c()))
        vlist = vlist[retain]

        ## Do noise-added inference
        pvs = sapply(vlist, function(v){
            pv = randomize_addnoise(y= y, v=v, sigma=sigma, numIS=numIS,
                                    sigma.add=sigma.add, orig.fudged.obj = h.fudged,
                                    numSteps = stoptime + consec,
                                    ic.poly = ic_obj$poly, bits=bits,
                                    inference.type=inference.type)
            return(pv)
        })
        names(pvs) = names(vlist)

        return(pvs)
    }, mc.cores=8)
    res = results[sapply(results, class)!="character"]
    all.pvs = (unlist(res))
    expect_equal(ks.test(all.pvs,"punif")$p.value > 0.05, TRUE)
})

