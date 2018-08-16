## Synopsis: Simulation code to compare each method's power. To be run from compare-run.R

dosim <- function(lev, ichunk, nsim, n=200, meanfun=fourjump, mc.cores=1,
                  numSteps=4, filename=NULL, sigma = 1, sigma.add=0.2, type,
                  outputdir = "../output"){

    assert_that(all(type %in% c("bsfs","nbsfs", "mbsfs", "wbsfs","mwbsfs",
                                "cbsfs","ncbsfs", "mbsfs", "fl","nfl", "mfl")))

    cat("lev=", lev, " and ichunk", ichunk, fill=TRUE)
    
    onesim <- function(isim){
        ## Generate data
        ## if(isim==5)browser()
        mn = meanfun(lev=lev, n=n)
        set.seed(isim)
        y = mn + rnorm(n, 0, 1)
        
        y.addnoise = rnorm(n, 0, sigma.add)
        results = list()

        ## Global settings
        ## max.numSteps = 10
        ## allsteps = allsteps.plain = 2:max.numSteps
        ## allsteps.cbs = 1:(max.numSteps/2) ## CBS should take half as many steps!
        ## allsteps.marg = 4 
        ## allsteps.cbs.marg = 2


        max.numSteps = 4
        allsteps = allsteps.plain = 4##2:max.numSteps

        ## Plain BS inference
        if(any(type=="bsfs")){tryCatch({
            ## Collect largest algorithm information
            obj = bsfs(y, numSteps=max.numSteps)
            poly.max = polyhedra(obj, numSteps=max.numSteps, record.nrows=TRUE)
            ## Collect each steps' inferences
            res = plain_inf_multistep(obj, allsteps, poly.max, mn, sigma)
            results$bsfs = res$pvs.by.step
            results$bsfs_zero = res$zeros.by.step
        })}

        ## Plain noisy BS inference (nonmarginalized)
        if(any(type=="nbsfs")){tryCatch({
            ## Collect largest algorithm information
            obj = bsfs(y=y, y.addnoise=y.addnoise, numSteps=max.numSteps, sigma.add=sigma.add)
            poly.max = polyhedra(obj, numSteps=max.numSteps, record.nrows=TRUE)
            ## Collect each steps' inferences
            res = plain_inf_multistep(obj, allsteps, poly.max, mn, sigma,
                                      shift=y.addnoise)
            results$nbsfs = res$pvs.by.step
            results$nbsfs_zero = res$zeros.by.step
            }, error=function(err){ print('error occurred during noisy bsfs')})
        }

        ## Marginalized noisy BS inference
        if(any(type=="mbsfs")){tryCatch({
            mbsfs = mbsfs_zero = vector("list", length(allsteps.marg))
            names(mbsfs) = names(mbsfs_zero) = paste0("step-",allsteps.marg)
            for(ii in 1:length(allsteps.marg)){
                numSteps = allsteps.marg[ii]
                obj = bsfs(y, numSteps=numSteps, sigma.add=sigma.add,
                           y.addnoise=y.addnoise)
                obj = addpv(obj, sigma=1, sigma.add=sigma.add, type="addnoise", mn=mn)
                mbsfs[[ii]] = obj$pvs
                mbsfs_zero[[ii]] = (obj$means==0)
            }
            results$mbsfs = mbsfs
            results$mbsfs_zero = mbsfs_zero
        }, error=function(err){ print('error occurred during noisy mbsfs')})
        }
       
 
        ## Plain (noisy) WBS inference
        if(any(type=="wbsfs")){tryCatch({
            ## Collect largest algorithm information
            obj = wbsfs(y, numSteps=max.numSteps, numIntervals=length(y))
            poly.max = polyhedra(obj, numSteps=max.numSteps, record.nrows=TRUE)
            ## Collect each steps' inferences
            res = plain_inf_multistep(obj, allsteps, poly.max, mn, sigma)
            results$wbsfs = res$pvs.by.step
            results$wbsfs_zero = res$zeros.by.step
        }, error=function(err){ print('error occurred during plain wbsfs')})}

        ## Marginalized WBS inference
        if(any(type=="mwbsfs")){tryCatch({
            mwbsfs = mwbsfs_zero = vector("list", length(allsteps.marg))
            names(mwbsfs) = names(mwbsfs_zero) = paste0("step-",allsteps.marg)
            for(ii in 1:length(allsteps.marg)){
                numSteps = allsteps.marg[ii]
                ## These two lines change, and nothing else
                obj = wbsfs(y, numSteps=numSteps, numIntervals=length(y))
                obj = addpv(obj, sigma=1, type="rand", mn=mn)
                mwbsfs[[ii]] = obj$pvs
                mwbsfs_zero[[ii]] = (obj$means==0)
            }
            results$mwbsfs = mwbsfs
            results$mwbsfs_zero = mwbsfs_zero

        }, error=function(err){ print('error occurred during marginalized wbsfs')})}
        
        ## Plain CBS inference (Non-marginalized)
        if(any(type=="cbsfs")){tryCatch({
            ## Collect largest algorithm information
            obj = cbsfs(y, numSteps=max.numSteps/2)
            poly.max = polyhedra(obj, numSteps=max.numSteps/2, record.nrows=TRUE)
            ## Collect each steps' inferences
            res = plain_inf_multistep(obj, allsteps.cbs, poly.max, mn, sigma)
            results$cbsfs = res$pvs.by.step
            results$cbsfs_zero = res$zeros.by.step
        }, error=function(err){ print('error occurred during plain cbsfs')})}

        ## Noisy CBS inference
        if(any(type=="ncbsfs")){tryCatch({
            ## Collect largest algorithm information
            obj = cbsfs(y=y, y.addnoise=y.addnoise, numSteps=max.numSteps/2,
                        sigma.add=sigma.add)

            poly.max = polyhedra(obj, numSteps=max.numSteps/2, record.nrows=TRUE)
            ## allsteps.cbs = process(allsteps.cbs)
            ## if(max.numSteps == length(obj$cp)) allsteps.cbs

            ## Collect each steps' inferences
            res = plain_inf_multistep(obj, allsteps.cbs, poly.max, mn, sigma,
                                      shift=y.addnoise)
            results$ncbsfs = res$pvs.by.step
            results$ncbsfs_zero = res$zeros.by.step
        }, error=function(err){ print('error occurred during ncbsfs')})}
        
        ## Marginalized noisy CBS inference
        if(any(type=="mcbsfs")){tryCatch({
            ## Number of steps are really tricky
            mcbsfs = mcbsfs_zero = vector("list", length(allsteps.cbs.marg))
            names(mcbsfs) = names(mcbsfs_zero) = paste0("step-",allsteps.cbs.marg)
            for(ii in 1:length(allsteps.cbs.marg)){
                numSteps = allsteps.cbs.marg[ii]
                ## These two lines change, and nothing else
                obj = cbsfs(y, numSteps=numSteps, sigma.add=sigma.add)
                obj = addpv(obj, sigma=1, sigma.add=sigma.add, type="addnoise", mn=mn)
                mcbsfs[[ii]] = obj$pvs
                mcbsfs_zero[[ii]] = (obj$means==0)
            }
            results$mcbsfs = mcbsfs
            results$mcbsfs_zero = mcbsfs_zero
        }, error=function(err){ print('error occurred during mcbsfs')})}
        
        ## Plain FL inference 
        if(any(type=="fl")){tryCatch({
            ## Collect largest algorithm information
            obj = fl(y, numSteps=max.numSteps)
            poly.max = polyhedra(obj, numSteps=max.numSteps)
            ## Collect each steps' inferences
            res = plain_inf_multistep(obj, allsteps, poly.max, mn, sigma)
            results$fl = res$pvs.by.step
            results$fl_zero = res$zeros.by.step
        }, error=function(err){ print('error occurred during plain fl')})}
        
        ## Noisy FL inference (Non-marginalized)
        if(any(type=="nfl")){tryCatch({
            ## Collect largest algorithm information
            obj = fl(y=y, y.addnoise=y.addnoise, numSteps=max.numSteps, sigma.add=sigma.add)
            poly.max = polyhedra.path(obj, numSteps=max.numSteps)
            ## Collect each steps' inferences
            res = plain_inf_multistep(obj, allsteps, poly.max, mn, sigma,
                                      shift=y.addnoise)
            results$nfl = res$pvs.by.step
            results$nfl_zero = res$zeros.by.step
        }, error=function(err){ print('error occurred during noisy fl')})}


        ## Marginalized noisy FL inference
        if(any(type=="mfl")){tryCatch({
            mfl = mfl_zero = vector("list", length(allsteps))
            names(mfl) = names(mfl_zero) = paste0("step-",allsteps)
            for(ii in 1:length(allsteps)){
                numSteps = allsteps[ii]
                ## These two lines change, and nothing else
                obj = fl(y, numSteps=numSteps, sigma.add=sigma.add)
                obj = addpv(obj, sigma=1, sigma.add=sigma.add, type="addnoise", mn=mn)
                mfl[[ii]] = obj$pvs
                mfl_zero[[ii]] = (obj$means==0)
            }
            results$mfl = mfl
            results$mfl_zero = mfl_zero

        }, error=function(err){ print(paste0('error occurred during noisy fl isim=', isim)) })}
        
        return(results)
        }

    ## Run the actual simulations.
    start.time = Sys.time()
    results.list = Mclapply(nsim, onesim, mc.cores, Sys.time())

    ## Save results
    if(is.null(filename)){ filename = paste0("compare-power-multistep-lev-",
                                             myfractions(lev), "-ichunk-", ichunk, ".Rdata")}
    save(results.list, file=file.path(outputdir, filename))
}
        


##' Helper for dosim(), in power comparison simulation. Collect information
##' (plain saturated p-values, zero-ness of mean) for all steps in |allsteps|.
##' @param allsteps All algorithm steps to test.
plain_inf_multistep <- function(obj, allsteps, poly.max, mn, sigma, shift=NULL){

    pvs.by.step = lapply(allsteps, function(numSteps){
        poly_pval2_from_vlist(y=obj$y.orig,
                              poly=snapshot(poly.max, numSteps),
                              vlist=make_all_segment_contrasts(obj, numSteps=numSteps),
                              sigma=sigma,
                              shift=shift)
    })

    zeros.by.step = lapply(allsteps, function(numSteps){
        vlist = make_all_segment_contrasts(obj, numSteps=numSteps)
        sapply(vlist, function(v) sum(v*mn)==0 )     })
    names(pvs.by.step) = names(zeros.by.step) = paste0("step-", allsteps)

    return(list(pvs.by.step=pvs.by.step, zeros.by.step=zeros.by.step))
}
