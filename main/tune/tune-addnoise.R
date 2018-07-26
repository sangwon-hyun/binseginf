## Synopsis: See how detection/power depend on amount of additive
## noise. Ultimate goal is to show a phase transition type phenomenon.

dosim <- function(sigma.add, nsim, mc.cores=1){
    cat("sigma.add=", sigma.add, fill=TRUE)
    outputdir = "../output"
    onesim <- function(isim){

        ## Generate data
        n = 200
        mn = fourjump(lev=1, n=n)
        y = mn + rnorm(n, 0, 1)

        ## Conduct noisy BS inference
        tryCatch({
        obj = bsfs(y, numSteps=4, sigma.add=sigma.add)
        obj = addpv(obj, sigma=1, mn=mn)
        pvs = obj$pvs
        }, error=function(err){ print('error occurred during noisy bsfs')})

        return(pvs)
    }

    ## Run the actual simulations
    start.time = Sys.time()
    results.list = mclapply(1:nsim, function(isim){
        printprogress(isim, nsim, start.time=start.time)
        onesim(isim)
    }, mc.cores=mc.cores, mc.preschedule=FALSE)

    filename = paste0("tune-addnoise-sigma-add-", round(sigma.add,3), ".Rdata")
    save(results.list, file=file.path(outputdir, filename))
}
