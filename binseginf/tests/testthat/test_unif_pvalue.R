context("Test whether we get null uniform p-values")

test_that("Null p-values are uniform under several scenarios.",{


    dosim <- function(lev, nsim, meanfun=onejump, mc.cores=1, numSteps=1,
                      ic.stop=FALSE,  alg, type){
        onesim <- function(isim){

            ## Base settings
            sigma.add = (if(type=="plain") 0 else 0.2)
            n = 10
            numIntervals = 20
            maxSteps = 10

            ## Generate mean and data
            mn = meanfun(lev=lev, n=n)
            y = mn + rnorm(n)

            ## Fit algorithm
            obj = alg(y=y, numSteps=numSteps, sigma.add=sigma.add, numIntervals=numIntervals)
            ## if(ic.stop){if(obj$ic_flag!="normal") return()}

            ## Conduct inference
            obj = addpv(obj, sigma=1, sigma.add=sigma.add, numIntervals=numIntervals, type=type, mn=mn)

            ## Store results
            results = list()
            results$pvs = obj$pvs
            results$zero = (obj$means==0)
            return(results)
        }
        results.list = Mclapply(nsim, onesim, mc.cores, Sys.time(), beep=TRUE)
        return(results.list)
    }

    testsim <- function(results.list){
        null.pvals = sapply(results.list, function(a){ a$pvs[a$zero] })
        expect_uniform(null.pvals)
        return(null.pvals)
    }

    ## (1) Two combinations: (lev=0 + numSteps=1), (lev=0 + numSteps=2)
    ## (2) With and without additive noise.
    ## (3) With and without IC stopping. (not coded for WBSFS, CBSFS, BSFS yet)

    mc.cores = 8
    algs = list(bsfs, fl, cbsfs, wbsfs)
    for(ic.stop in c(TRUE, FALSE)){
        for(i.alg in 1:4){
            alg = algs[[i.alg]]
            a = dosim(lev=0, nsim=5000, numSteps=1, alg=alg, ic.stop=ic.stop, mc.cores=mc.cores, type="plain")
            a = dosim(lev=0, nsim=5000, numSteps=1, alg=alg, ic.stop=ic.stop, mc.cores=mc.cores, type=(if(i.alg==4)"rand" else "addnoise"))
            a = dosim(lev=2, nsim=5000, numSteps=2, alg=alg, ic.stop=ic.stop, mc.cores=mc.cores, type="plain")
            a = dosim(lev=2, nsim=5000, numSteps=2, alg=alg, ic.stop=ic.stop, mc.cores=mc.cores, type=(if(i.alg==4)"rand" else "addnoise"))
        }
    }


})







## test_that("Separate test for whether FL null p-values are all uniform", {

##     ## Testing FL in flat mean
##     onesim <- function(isim){
##         n = 10
##         mn = rep(0, n)
##         y = mn + rnorm(n, 0, 1)
##         obj = fl(y=y, numSteps=1)
##         obj = addpv_fl(obj, sigma=1, type="plain", mn=mn)
##         return(obj$pvs)
##     }
##     pvslist = Mclapply(nsim=1000, onesim, 4, Sys.time())
##     qqunif(unlist(pvslist))
##     expect_uniform(unlist(pvslist))

##     ## Testing FL in nonflat mean
##     onesim <- function(isim){
##         n = 10
##         mn = c(rep(0, n/2), rep(2, n/2))
##         y = mn + rnorm(n, 0, 1)
##         obj = fl(y=y, numSteps=2)
##         obj = addpv_fl(obj, sigma=1, type="plain", mn=mn)
##         return(obj$pvs[which(obj$means==0)])
##     }
##     pvslist = Mclapply(nsim=4000, onesim, 4, Sys.time())
##     qqunif(unlist(pvslist))
##     expect_uniform(unlist(pvslist))

## })


## test_that("Uniform null FL p-values without ic stopping.",{

##     dosim <- function(lev, nsim, n=200, meanfun=onejump, mc.cores=1,
##                       numSteps=4, maxSteps=n/3,
##                       ic.stop=FALSE){
##         onesim <- function(isim){
##             mn = meanfun(lev=lev, n=n)
##             y = mn + rnorm(n)
##             results = list()
##             obj = fl(y=y, maxSteps=maxSteps, numSteps=numSteps, sigma.add=0.2,
##                      ic.stop=ic.stop)
##             if(ic.stop){if(obj$ic_flag!="normal") return()}
##             obj = addpv_fl(obj, sigma=1, sigma.add=0.2, type="addnoise", mn=mn)
##             results$bsfs_addnoise = obj$pvs
##             results$bsfs_addnoise_zero = (obj$means==0)
##             return(results)
##         }
##         results.list = Mclapply(nsim, onesim, mc.cores, Sys.time(), beep=TRUE)
##         return(results.list)
##     }

##     ## continue here!!
##     a = dosim(lev=0, nsim=600, n=10, numSteps=1, mc.cores=4)
##     a = dosim(lev=0, nsim=10000, maxSteps=10, n=20, numSteps=1, mc.cores=8)
##     save(a, file=file.path("../output", "fix-fl.rdata"))
##     load(file=file.path("~/desktop/tempoutput/fix-fl.rdata"))
##     aa = a[sapply(a,length)!=0]
##     obj = a[[1]]
##     aa = sapply(a, function(obj){obj[[1]][which(obj[[2]])]})
##     qqunif(aa)

## })

## test_that("Uniform null BS p-values without ic stopping.",{

##     dosim <- function(lev, nsim, n=200, meanfun=onejump, mc.cores=1,
##                       numSteps=4, maxSteps=n/3,
##                       ic.stop=FALSE, sigma.add=0.2){
##         onesim <- function(isim){
##             mn = meanfun(lev=lev, n=n)
##             y = mn + rnorm(n)
##             results = list()
##             obj = bsfs(y=y, numSteps=numSteps, sigma.add=sigma.add)
##             if(ic.stop){if(obj$ic_flag!="normal") return()}
##             obj = addpv(obj, sigma=1, sigma.add=sigma.add, type="addnoise", mn=mn)
##             results$pvs = obj$pvs
##             results$zero = (obj$means==0)
##             return(results)
##         }
##         results.list = Mclapply(nsim, onesim, mc.cores, Sys.time(), beep=TRUE)
##         return(results.list)
##     }

##     ## Flat
##     results = dosim(lev=0, nsim=2000, n=10, numSteps=1, mc.cores=8)
##     null.pvals = sapply(results, function(a){ a$pvs[a$zero] })
##     expect_uniform(null.pvals)


##     ## Nonflat
##     results = dosim(lev=2, nsim=2000, n=10, numSteps=2, mc.cores=8)
##     null.pvals = sapply(results, function(a){ a$pvs[a$zero] })
##     expect_uniform(null.pvals)
    
## })
