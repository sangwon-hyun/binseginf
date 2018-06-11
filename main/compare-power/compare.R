## Synopsis: Simulation code to compare each method's power. To be run from compare-run.R
dosim <- function(lev, ichunk, nsim, n=200, meanfun=fourjump, mc.cores=1, numSteps=4, filename=NULL){
    cat("lev=", lev, " and ichunk", ichunk, fill=TRUE)
    outputdir = "../output"
    onesim <- function(isim){

        ## Generate data
        ## n = 200
        mn = meanfun(lev=lev, n=n)
        y = mn + rnorm(n, 0, 1)
        results = list()


        ## ## One attempt at /unifying/ the code
        ## obj = alg(y, numSteps=numSteps, numIntervals=numIntervals, sigma.add=sigma.add)
        ## obj1 = addpv(obj, sigma=1, mn=mn)
        ## if(!change) ##don't run obj2
        ## obj2 = addpv(obj, sigma=1, mn=mn, with.decluttering)
        ## results$bsfs_plain= obj$pvs
        ## results$bsfs_plain_zero = (obj$means==0)
        ## }, error=function(err){ print('error occurred during plain bsfs')})

        ## Plain BS inference
        tryCatch({
        obj = bsfs(y, numSteps=numSteps)
        obj = addpv(obj, sigma=1, mn=mn)
        results$bsfs_plain= obj$pvs
        results$bsfs_plain_zero = (obj$means==0)
        }, error=function(err){ print('error occurred during plain bsfs')})
        
        ## Noisy BS inference
        tryCatch({
        obj = bsfs(y, numSteps=numSteps, sigma.add=0.2)
        obj = addpv(obj, sigma=1, sigma.add=0.2, type="addnoise", mn=mn)
        results$bsfs_addnoise = obj$pvs
        results$bsfs_addnoise_zero = (obj$means==0)
        }, error=function(err){ print('error occurred during noisy bsfs')})
        
        ## Plain WBS inference
        tryCatch({
        obj = wbsfs(y, numSteps=numSteps, numIntervals=length(y))
        obj = addpv(obj, sigma=1, type="plain", mn=mn)
        results$wbsfs_plain = obj$pvs
        results$wbsfs_plain_zero = (obj$means==0)
        }, error=function(err){ print('error occurred during plain wbsfs')})
        
        ## Marginalized WBS inference
        tryCatch({
        obj = wbsfs(y, numSteps=numSteps, numIntervals=length(y))
        obj = addpv(obj, sigma=1, type="rand", mn=mn)
        results$wbsfs_marg = obj$pvs
        results$wbsfs_marg_zero = (obj$means==0)
        }, error=function(err){ print('error occurred during marginalized wbsfs')})
        
        ## Plain CBS inference
        tryCatch({
        obj = cbsfs(y, numSteps=numSteps/2)
        obj = addpv(obj, sigma=1, type="plain", mn=mn)
        results$cbsfs_plain = obj$pvs
        results$cbsfs_plain_zero = (obj$means==0)
        }, error=function(err){ print('error occurred during plain cbsfs')})
        
        ## Noisy CBS inference
        tryCatch({
        obj = cbsfs(y, numSteps=numSteps/2, sigma.add=0.2)
        obj = addpv(obj, sigma=1, type="addnoise", mn=mn, sigma.add=0.2)
        results$cbsfs_addnoise = obj$pvs
        results$cbsfs_addnoise_zero = (obj$means==0)
        }, error=function(err){ print('error occurred during noisy cbsfs')})
        
        ## Plain FL inference 
        tryCatch({
        obj = fl(y, numSteps=numSteps)
        obj = addpv_fl(obj, sigma=1, type="plain", mn=mn)
        results$fl_plain = obj$pvs
        results$fl_plain_zero = (obj$means==0)
        }, error=function(err){ print('error occurred during plain fl')})
        
        ## Noisy FL inference
        tryCatch({
        obj = fl(y, numSteps=numSteps, sigma.add=0.2) 
        obj = addpv_fl(obj, sigma=1, sigma.add=0.1, type="addnoise", mn=mn)
        results$fl_addnoise = obj$pvs
        results$fl_addnoise_zero = (obj$means==0)
        }, error=function(err){ print('error occurred during noisy fl')})

        ## ## IC-stopped FL inference (with decluttering? not written yet)
        ## tryCatch({
        ## obj = fl(y, numSteps=numSteps, sigma.add=0.2, ic.stop=TRUE) 
        ## obj = addpv_fl(obj, sigma=1, sigma.add=0.1, type="addnoise", mn=mn)
        ## results$fl_ic_addnoise = obj$pvs
        ## results$fl_ic_addnoise_zero = (obj$means==0)
        ## print('fl')
        ## }, error=function(err){ print('error occurred during noisy fl')})
        ##
        
        return(results)
    }

    ## Run the actual simulations
    start.time = Sys.time()
    results.list = Mclapply(nsim, onesim, mc.cores, Sys.time())

    ## Save or return
    ## if(!is.null(filename)) filename = paste0("compare-power-fourjump-lev-",
    ##                                          myfractions(lev), ".Rdata")

    if(!is.null(filename)) filename = paste0("compare-power-fourjump-lev-",
                                             myfractions(lev), "-ichunk-", ichunk, ".Rdata")
    save(results.list, file=file.path(outputdir, filename))
   ## return(results.list)
}
