## Synopsis: Simulation code to compare each method's power. To be run from compare-run.R

dosim <- function(lev, nsim, mc.cores=1){
    cat("lev=", lev, fill=TRUE)
    outputdir = "../output"
    onesim <- function(isim){

        ## Generate data
        n = 20
        mn = fourjump(lev=lev, n=n)
        y = mn + rnorm(n, 0, 1)
        results = list()


        ## Plain BS inference
        tryCatch({
        obj = bsfs(y, numSteps=4)
        obj = addpv(obj, sigma=1, mn=mn)
        results$bsfs_plain= obj$pvs
        results$bsfs_plain_zero = (obj$means==0)
        }, error=function(err){ print('error occurred during plain bsfs')})
        
        ## Noisy BS inference
        tryCatch({
        obj = bsfs(y, numSteps=4, sigma.add=0.2)
        obj = addpv(obj, sigma=1, sigma.add=0.2, type="addnoise", mn=mn)
        results$bsfs_addnoise = obj$pvs
        results$bsfs_addnoise_zero = (obj$means==0)
        }, error=function(err){ print('error occurred during noisy bsfs')})
        
        ## Plain WBS inference
        tryCatch({
        obj = wbsfs(y, numSteps=4, numIntervals=length(y))
        obj = addpv(obj, sigma=1, type="plain", mn=mn)
        results$wbsfs_plain = obj$pvs
        results$wbsfs_plain_zero = (obj$means==0)
        }, error=function(err){ print('error occurred during plain wbsfs')})
        
        ## Marginalized WBS inference
        tryCatch({
        obj = wbsfs(y, numSteps=4, numIntervals=length(y))
        obj = addpv(obj, sigma=1, type="rand", mn=mn)
        results$wbsfs_marg = obj$pvs
        results$wbsfs_marg_zero = (obj$means==0)
        }, error=function(err){ print('error occurred during marginalized wbsfs')})
        
        ## Plain CBS inference
        tryCatch({
        obj = cbsfs(y, numSteps=2)
        obj = addpv(obj, sigma=1, type="plain", mn=mn)
        results$cbsfs_plain = obj$pvs
        results$cbsfs_plain_zero = (obj$means==0)
        }, error=function(err){ print('error occurred during plain cbsfs')})
        
        ## Noisy CBS inference
        tryCatch({
        obj = cbsfs(y, numSteps=2, sigma.add=0.2)
        obj = addpv(obj, sigma=1, type="addnoise", mn=mn, sigma.add=0.2)
        results$cbsfs_addnoise = obj$pvs
        results$cbsfs_addnoise_zero = (obj$means==0)
        }, error=function(err){ print('error occurred during noisy cbsfs')})
        
        ## Plain FL inference 
        tryCatch({
        obj = fl(y, numSteps=4)
        obj = addpv_fl(obj, sigma=1, type="plain", mn=mn)
        results$fl_plain = obj$pvs
        results$fl_plain_zero = (obj$means==0)
        }, error=function(err){ print('error occurred during plain fl')})
        
        ## Noisy FL inference
        tryCatch({
        obj = fl(y, numSteps=1, sigma.add=0.2)
        obj = addpv_fl(obj, sigma=1, sigma.add=0.1, type="addnoise", mn=mn)
        results$fl_addnoise = obj$pvs
        results$fl_addnoise_zero = (obj$means==0)
        }, error=function(err){ print('error occurred during noisy fl')})

        ## IC-stopped FL inference (not written yet)
        
        return(results)
    }

    ## Run the actual simulations
    start.time = Sys.time()
    results.list = mclapply(1:nsim, function(isim){
        printprogress(isim, nsim, start.time=start.time)
        onesim(isim)
    }, mc.cores=mc.cores, mc.preschedule=FALSE)

    filename = paste0("compare-power-lev", lev, ".Rdata")
    save(results.list, file=file.path(outputdir, filename))
    ## return(results.list)
}

