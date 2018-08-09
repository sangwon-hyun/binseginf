## Synopsis: Simulation code to compare each method's power. To be run from compare-run.R

dosim <- function(lev, ichunk, nsim, n=200, meanfun=fourjump, mc.cores=1,
                  numSteps=4, filename=NULL, sigma.add=0.2, type,
                  max.numSteps = 10,
                  allsteps = 2:max.numSteps,
                  allsteps.cbs = 1:(max.numSteps/2)
                  ){

    assert_that(all(type %in% c("bsfs","nbsfs","wbsfs","mwbsfs","cbsfs","ncbsfs","fl","nfl")))

    cat("lev=", lev, " and ichunk", ichunk, fill=TRUE)
    outputdir = "../output"
    onesim <- function(isim){

        ## Generate data
        mn = meanfun(lev=lev, n=n)
        y = mn + rnorm(n, 0, 1)
        results = list()

        ## Global settings
        ## max.numSteps = 10
        ## allsteps = 2:max.numSteps
        ## ## allsteps = max.numSteps = 10 ## temporary
        ## allsteps.cbs = 1:(max.numSteps/2)

        ## Plain BS inference
        if(any(type=="bsfs")){tryCatch({
            ## Collect largest algorithm information
            obj = bsfs(y, numSteps=max.numSteps)
            poly.max = polyhedra(obj, numSteps=max.numSteps, record.nrows=TRUE)
            ## Collect each steps' inferences
            res = inf.by.step(obj, allsteps, poly.max, mn)
            results$bsfs = res$pvs.by.step
            results$bsfs_zero = res$zeros.by.step
        })}

        ## Noisy BS inference (under construction)
        if(any(type=="nbsfs")){tryCatch({
            nbsfs = nbsfs_zero = vector("list", length(allsteps))
            names(nbsfs) = names(nbsfs_zero) = paste0("step-",allsteps)
            for(ii in 1:length(allsteps)){
                numSteps = allsteps[ii]
                ## These two lines change, and nothing else
                obj = bsfs(y, numSteps=numSteps, sigma.add=sigma.add)
                obj = addpv(obj, sigma=1, sigma.add=sigma.add, type="addnoise", mn=mn)
                nbsfs[[ii]] = obj$pvs
                nbsfs_zero[[ii]] = (obj$means==0)
            }
            results$nbsfs = nbsfs
            results$nbsfs_zero = nbsfs_zero
        }, error=function(err){ print('error occurred during noisy bsfs')})
        }
        
        ## Plain WBS inference
        if(any(type=="wbsfs")){tryCatch({
            ## Collect largest algorithm information
            obj = wbsfs(y, numSteps=max.numSteps, numIntervals=length(y))
            poly.max = polyhedra(obj, numSteps=max.numSteps, record.nrows=TRUE)
            ## Collect each steps' inferences
            res = inf.by.step(obj, allsteps, poly.max, mn)
            results$wbsfs = res$pvs.by.step
            results$wbsfs_zero = res$zeros.by.step
        }, error=function(err){ print('error occurred during plain wbsfs')})}

        ## ## Marginalized WBS inference
        if(any(type=="mwbsfs")){tryCatch({

            mwbsfs = mwbsfs_zero = vector("list", length(allsteps))
            names(mwbsfs) = names(mwbsfs_zero) = paste0("step-",allsteps)
            for(ii in 1:length(allsteps)){
                numSteps = allsteps[ii]
                ## These two lines change, and nothing else
                obj = wbsfs(y, numSteps=numSteps, numIntervals=length(y))
                obj = addpv(obj, sigma=1, type="rand", mn=mn)
                mwbsfs[[ii]] = obj$pvs
                mwbsfs_zero[[ii]] = (obj$means==0)
            }
            results$mwbsfs = mwbsfs
            results$mwbsfs_zero = mwbsfs_zero

        }, error=function(err){ print('error occurred during marginalized wbsfs')})}
        
        ## Plain CBS inference
        if(any(type=="cbsfs")){tryCatch({
            ## Collect largest algorithm information
            obj = cbsfs(y, numSteps=max.numSteps/2)
            poly.max = polyhedra(obj, numSteps=max.numSteps/2, record.nrows=TRUE)
            ## Collect each steps' inferences
            res = inf.by.step(obj, allsteps.cbs, poly.max, mn)
            results$cbsfs = res$pvs.by.step
            results$cbsfs_zero = res$zeros.by.step
        }, error=function(err){ print('error occurred during plain cbsfs')})}
        
        ## ## Noisy CBS inference
        if(any(type=="ncbsfs")){tryCatch({
            ## Number of steps are really tricky
            ncbsfs = ncbsfs_zero = vector("list", length(allsteps))
            names(ncbsfs) = names(ncbsfs_zero) = paste0("step-",allsteps)
            for(ii in 1:length(allsteps.cbs)){
                numSteps = allsteps.cbs[ii]
                ## These two lines change, and nothing else
                obj = cbsfs(y, numSteps=numSteps, sigma.add=sigma.add)
                obj = addpv(obj, sigma=1, sigma.add=sigma.add, type="addnoise", mn=mn)
                ncbsfs[[ii]] = obj$pvs
                ncbsfs_zero[[ii]] = (obj$means==0)
            }
            results$ncbsfs = ncbsfs
            results$ncbsfs_zero = ncbsfs_zero


        }, error=function(err){ print('error occurred during noisy cbsfs')})}
        
        ## Plain FL inference 
        if(any(type=="fl")){tryCatch({
            ## Collect largest algorithm information
            obj = fl(y, numSteps=max.numSteps)
            poly.max = polyhedra(obj, numSteps=max.numSteps)
            ## Collect each steps' inferences
            res = inf.by.step(obj, allsteps, poly.max, mn)
            results$fl = res$pvs.by.step
            results$fl_zero = res$zeros.by.step
        }, error=function(err){ print('error occurred during plain fl')})}
        
        ## ## Noisy FL inference
        if(any(type=="nfl")){tryCatch({
            nfl = nfl_zero = vector("list", length(allsteps))
            names(nfl) = names(nfl_zero) = paste0("step-",allsteps)
            for(ii in 1:length(allsteps)){
                print(ii)
                numSteps = allsteps[ii]
                ## These two lines change, and nothing else
                obj = fl(y, numSteps=numSteps, sigma.add=sigma.add)
                obj = addpv(obj, sigma=1, sigma.add=sigma.add, type="addnoise", mn=mn)
                nfl[[ii]] = obj$pvs
                nfl_zero[[ii]] = (obj$means==0)
            }
            results$nfl = nfl
            results$nfl_zero = nfl_zero

        }, error=function(err){ print(paste0('error occurred during noisy fl isim=', isim)) })}
        
        return(results)
        }

    ## Run the actual simulations
    start.time = Sys.time()
    results.list = Mclapply(nsim, onesim, mc.cores, Sys.time())

    ## Save
    if(is.null(filename)){ filename = paste0("compare-power-multistep-", type, "-lev-",
                                             myfractions(lev), "-ichunk-", ichunk, ".Rdata")}
    save(results.list, file=file.path(outputdir, filename))
}



##' Helper to get a list of p-values at each step |numSteps|.
get.plain.pv.by.step <- function(obj, numSteps, poly.max, sigma=1){
    poly.this.step = snapshot(poly.max, numSteps)
    if(class(obj)=="cbsfs") numSteps = numSteps * 2
    vlist = make_all_segment_contrasts(obj, numSteps=numSteps)
    return(poly_pval2_from_vlist(obj$y, poly.this.step, vlist, sigma))
}

##' Collect information for all steps in |allsteps|.
inf.by.step <- function(obj, allsteps, poly.max, mn){

    pvs.by.step = lapply(allsteps, function(numSteps){
        a = get.plain.pv.by.step(obj, numSteps, poly.max, 1)
    })

    if(class(obj)=="cbsfs") allsteps = allsteps * 2
    zeros.by.step = sapply(allsteps, function(numSteps){
        vlist = make_all_segment_contrasts(obj, numSteps=numSteps)
        sapply(vlist, function(v) sum(v*mn)==0 )     })
    names(pvs.by.step) = names(zeros.by.step) = paste0("step-",allsteps)

    return(list(pvs.by.step=pvs.by.step, zeros.by.step=zeros.by.step))
}
