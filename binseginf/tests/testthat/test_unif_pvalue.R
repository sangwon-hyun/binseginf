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
        results.list = mclapply(1:nsim, onesim, mc.cores=mc.cores)
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
