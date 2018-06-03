<<<<<<< Updated upstream
test_that("Reduction of inference procedure into collecting Gy and Gv correctly works",{

    n=60
    lev = 0
    mn = c(rep(0,n/2), rep(lev,n/2))
    sigma = 1
    y = mn + rnorm(n, 0, sigma)
    numIntervals = 10
    numSteps = 3
    g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps, cumsum.y=cumsum.y, cumsum.v=cumsum.v, inference.type='rows')
    vlist <- make_all_segment_contrasts(g)
    v=vlist[[1]]
    cumsum.y = cumsum(y)
    cumsum.v = cumsum(v)

    ## ## See if the results are equal
    ## all.equal(sort(as.numeric(g$gamma %*% y)), sort(as.numeric( g$Gy)))
    ## all.equal(sort(as.numeric(g$gamma %*% v)), sort(as.numeric( g$Gv)))
    ## poly = polyhedra(obj=g$gamma, u=g$u)

    ## See if there is speedup
    g1 = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps, cumsum.y=cumsum.y, cumsum.v=cumsum.v, inference.type="pre-multiply")
    g1 = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps, cumsum.y=cumsum.y, cumsum.v=cumsum.v, inference.type="rows")
    microbenchmark({g2=wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps, cumsum.y=cumsum.y, cumsum.v=cumsum.v, inference.type="pre-multiply")}, times=100)
    microbenchmark({g1=wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps, cumsum.y=cumsum.y, cumsum.v=cumsum.v, inference.type="rows")}, times=100)



test_that("WBS.FT gives uniform p-values",{
=======
test_that("WBS has uniform null p-values and power when signal is present.",{
>>>>>>> Stashed changes

    mysim = function(n, lev){

        nsim = 1000
        numSteps = 3
        sigma = 1

        pvs = mclapply(1:nsim,function(isim){

            ## Generate some data
            mn = c(rep(0,n/2), rep(lev,n/2))
            y = mn + rnorm(n, 0, sigma)

            ## Fit WBS
            numIntervals = 30
            g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps)
            poly = polyhedra(obj=g$gamma, u=g$u)
            vlist <- make_all_segment_contrasts(g)
            pvs = sapply(vlist, function(v){
                return(poly.pval2(y=y, poly=poly, v=v, sigma=sigma)$pv)
            })
            names(pvs) = g$cp*g$cp.sign
            return(pvs)
        },mc.cores=4)
        pvs = unlist(pvs)
        pvs = unique(pvs) ## Sometimes ties occur numerically. Not sure why.
    }

    null.pvals = sapply(c(10,20,50,100), mysim, lev=0)
    nonnull.pvals = sapply(c(10,20,50,100), mysim, lev=5)

    ## Do KS tests
    sapply(null.pvals, function(pvs){
        expect_true(ks.test(unlist(pvs), punif)$p.value > 0.05)
    })
    sapply(nonnull.pvals, function(pvs){
        expect_true(ks.test(unlist(pvs), punif)$p.value < 0.01)
    })

})


test_that("WBS polyhedron is exact",{

    ## Simulation settings
    n = 10
    lev = 0
    mn = c(rep(0,n/2), rep(lev,n/2))
    sigma = 1
    numIntervals = 30

    nsim=1000
    a = mclapply(1:nsim,function(isim){

        ## Generate original data.
        y.orig = mn + rnorm(n, 0, sigma)
        obj.orig = intervals(numIntervals=numIntervals, n=n)

        ## Fit WBS
        numSteps = 3
        g.orig = wildBinSeg_fixedSteps(y.orig, numIntervals=numIntervals, numSteps=numSteps, intervals=obj.orig)
        poly.orig = polyhedra(obj=g.orig$gamma, u=g.orig$u)

        ## Add noise and see if exact correspondence
        sigma.add = .2
        y.new = y.orig + rnorm(n,0,sigma.add)
        g.new = wildBinSeg_fixedSteps(y.new, numIntervals=numIntervals, numSteps=numSteps, intervals=obj.orig)

        ## Helper function to see if results match up, at every step.
        wbs_results_match <- function(obj.orig, obj.new){
            return(all.equal(obj.orig$results[,c("max.s", "max.b", "max.e", "max.sign")],
                             obj.new$results[,c("max.s", "max.b", "max.e", "max.sign")]))
        }
        ## See if the results are exactly the same.
        if(wbs_results_match(g.orig, g.new)==TRUE){
            expect_true(contained(poly.orig, y.new))
        }
    }, mc.cores=4)
})

test_that("Randomized wbs inference has uniform null p-values, and right-skewed nonnull p-values",{

    mysim <- function(n, lev, test.type=c("none", "ks.unif","ks.nonunif"),mc.cores=6){
        test.type = match.arg(test.type)
        cat("Sample size=", n, fill=TRUE)

        ## Simulation settings
        mn = c(rep(0,n/2), rep(lev,n/2))
        sigma = 1
        numIntervals = n
        numSteps=1
        numIS = 100

        ## Get randomized p-values
        nsim = 1000
        pvs = mclapply(1:nsim, function(isim){
            printprogress(isim,nsim)
            set.seed(isim)
            y = mn + rnorm(n, 0, sigma)
            g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps)
            poly = polyhedra(obj=g$gamma, u=g$u)
            vlist <- make_all_segment_contrasts(g)
            v = vlist[[1]]
            pv.rand = suppressWarnings(randomize_wbsfs(v=v, winning.wbs.obj=g,
                                                        sigma=sigma,
                                                        numIS = numIS))

        }, mc.cores=mc.cores)
        pvs = unlist(pvs)
        if(test.type=="ks.unif") expect_true(ks.test(unlist(pvs), punif)$p.value > 0.05)
        if(test.type=="ks.nonunif") expect_true(ks.test(unlist(pvs), punif)$p.value > 0.05)
        return(pvs)
    }
    null.pvals.by.samplesize = lapply(c(10,20,50,100), mysim, lev=0, test.type="ks.unif")
    nonnull.pvals.by.samplesize = lapply(c(10,20,50,100), mysim, lev=3, test.type="ks.nonunif")
})


test_that("Randomized wbs inference has uniform nulls in further steps (2+).",{

    mysim <- function(n, lev, test.type=c("none", "ks.unif","ks.nonunif"),mc.cores=6){
        test.type = match.arg(test.type)
        cat("Sample size=", n, fill=TRUE)

        ## Simulation settings
        mn = c(rep(0,n/2), rep(lev,n/2))
        sigma = 1
        numIntervals = n
        numSteps=4
        numIS = 100

        ## Get randomized p-values
        nsim = 100
        pvs = mclapply(1:nsim, function(isim){
            printprogress(isim,nsim)
            set.seed(isim)
            y = mn + rnorm(n, 0, sigma)
            g = wildBinSeg_fixedSteps(y, numIntervals=numIntervals, numSteps=numSteps)
            poly = polyhedra(obj=g$gamma, u=g$u)
            vlist <- make_all_segment_contrasts(g)
            ## v = vlist[[1]]
            pv.rands = sapply(vlist, function(v){
                suppressWarnings(randomize_wbsfs(v=v, winning.wbs.obj=g,
                                                 sigma=sigma,
                                                 numIS = numIS))
            }
            return(pv.rands)
        }, mc.cores=mc.cores)
        pvs = unlist(pvs)
        if(test.type=="ks.unif") expect_true(ks.test(unlist(pvs), punif)$p.value > 0.05)
        if(test.type=="ks.nonunif") expect_true(ks.test(unlist(pvs), punif)$p.value > 0.05)
        return(pvs)
    }
    null.pvals.by.samplesize = lapply(c(10,20,50,100), mysim, lev=0, test.type="ks.unif")

})



test_that("WBS output matched wbs::wbs()",{
print('not written yet')
})
