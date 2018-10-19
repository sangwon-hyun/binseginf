## Synopsis: Simulation code to compare each method's power. To be run from compare-run.R
dosim <- function(lev, ichunk, nsim, n=200, meanfun=fourjump, mc.cores=1,
                  numSteps=4, filename=NULL, sigma = 1, sigma.add=0.5, type,
                  outputdir = "../output", locs=1:n,

                  ## Additional settings
                  max.numSteps = 10,
                  allsteps=2:max.numSteps,
                  allsteps.cbs=1:(max.numSteps/2), ## CBS should take half as many steps!
                  allsteps.marg=4,
                  allsteps.cbs.marg=(allsteps.marg/2),## CBS should take half as many steps!

                  ## For decluttering
                  how.close=2,
                  seed.set=FALSE,

                  ##Temporary, for MFL
                  max.numIS.fl=2000,
                  min.num.things.fl=30,

                  ## multicore option
                  mc.preschedule=TRUE
                  ){

    assert_that(all(type %in% c("bsfs","nbsfs", "mbsfs", "wbsfs","mwbsfs",
                                "cbsfs","ncbsfs","mcbsfs", "fl","nfl", "mfl",
                                "dfl", "dbsfs", ## Temporarily added as separate
                                ## things, but should actually be
                                ## absorbed into plain inference
                                ## eventually.
                                "mdfl", "mdbsfs")),
                msg="|type| error")

    cat("lev=", lev, " and ichunk", ichunk, fill=TRUE)
    
    onesim <- function(isim){
        
        ## Generate data
        if(seed.set) set.seed(isim)
        mn = meanfun(lev=lev, n=n)
        y = mn + rnorm(n, 0, sigma)
        
        y.addnoise = rnorm(n, 0, sigma.add)
        results = list()


        ## max.numSteps = 4
        ## allsteps = allsteps.plain = 4

        ## Plain BS inference
        if(any(type=="bsfs")){tryCatch({
            obj = bsfs(y, numSteps=max.numSteps)
            poly.max = polyhedra(obj, numSteps=max.numSteps, record.nrows=TRUE)
            res = plain_inf_multistep(obj, allsteps, poly.max, mn, sigma, locs=locs)
            results$bsfs = res$pvs.by.step
            results$bsfs_zero = res$zeros.by.step
            results$bsfs_cps = obj$cp * obj$cp.sign 
        })}

        ## Plain noisy BS inference (nonmarginalized)
        if(any(type=="nbsfs")){tryCatch({
            obj = bsfs(y=y, y.addnoise=y.addnoise, numSteps=max.numSteps, sigma.add=sigma.add)
            poly.max = polyhedra(obj, numSteps=max.numSteps, record.nrows=TRUE)
            res = plain_inf_multistep(obj, allsteps, poly.max, mn, sigma,
                                      shift=y.addnoise, locs=locs)
            results$nbsfs = res$pvs.by.step
            results$nbsfs_zero = res$zeros.by.step
            results$nbsfs_cps = obj$cp * obj$cp.sign 
            }, error=function(err){ print('error occurred during noisy bsfs')})
        }

        ## "D"ecluttered plain BS inference
        if(any(type=="dbsfs")){tryCatch({
            obj = bsfs(y, numSteps=max.numSteps)
            ## poly.max = polyhedra(obj, numSteps=max.numSteps)
            poly.max = polyhedra(obj, numSteps=max.numSteps, record.nrows=TRUE)
            res = plain_inf_multistep(obj, allsteps, poly.max, mn, sigma,
                                      locs=locs, declutter=TRUE,
                                      how.close=how.close)
            results$dbsfs = res$pvs.by.step
            results$dbsfs_zero = res$zeros.by.step
            results$dbsfs_cps = obj$cp * obj$cp.sign 

        }, error=function(err){ print(paste0('error occurred during decluttered bsfs isim=', isim)) })}

        if(any(type=="mdbsfs")){tryCatch({
            mdbsfs = mdbsfs_zero = vector("list", length(allsteps.marg))
            names(mdbsfs) = names(mdbsfs_zero) = paste0("step-",allsteps.marg)
            for(ii in 1:length(allsteps.marg)){
                numSteps = allsteps.marg[ii]
                obj = bsfs(y, numSteps=numSteps, sigma.add=sigma.add,
                           y.addnoise=y.addnoise)
                obj = addpv(obj, sigma=sigma, sigma.add=sigma.add, type="addnoise", mn=mn,
                            locs=locs, declutter=TRUE)
                mdbsfs[[ii]] = obj$pvs
                mdbsfs_zero[[ii]] = (obj$means==0)
            }
            results$mdbsfs = mdbsfs
            results$mdbsfs_zero = mdbsfs_zero
            results$mdbsfs_cps = obj$cp * obj$cp.sign 
        }, error=function(err){ print('error occurred during noisy mdbsfs')})
        }
       
 


        ## Marginalized noisy BS inference
        if(any(type=="mbsfs")){tryCatch({
            mbsfs = mbsfs_zero = vector("list", length(allsteps.marg))
            names(mbsfs) = names(mbsfs_zero) = paste0("step-",allsteps.marg)
            for(ii in 1:length(allsteps.marg)){
                numSteps = allsteps.marg[ii]
                obj = bsfs(y, numSteps=numSteps, sigma.add=sigma.add,
                           y.addnoise=y.addnoise)
                obj = addpv(obj, sigma=sigma, sigma.add=sigma.add, type="addnoise", mn=mn,
                            locs=locs)
                mbsfs[[ii]] = obj$pvs
                mbsfs_zero[[ii]] = (obj$means==0)
            }
            results$mbsfs = mbsfs
            results$mbsfs_zero = mbsfs_zero
            results$mbsfs_cps = obj$cp * obj$cp.sign 
        }, error=function(err){ print('error occurred during noisy mbsfs')})
        }
       
 
        ## Plain (noisy) WBS inference
        if(any(type=="wbsfs")){tryCatch({
            obj = wbsfs(y, numSteps=max.numSteps, numIntervals=length(y))
            obj$y.orig = y
            poly.max = polyhedra(obj, numSteps=max.numSteps, record.nrows=TRUE)
            res = plain_inf_multistep(obj, allsteps, poly.max, mn, sigma, locs=locs)
            results$wbsfs = res$pvs.by.step
            results$wbsfs_zero = res$zeros.by.step
            results$wbsfs_cps = obj$cp * obj$cp.sign 
        }, error=function(err){ print('error occurred during plain wbsfs')})}

        ## Marginalized WBS inference
        if(any(type=="mwbsfs")){tryCatch({
            mwbsfs = mwbsfs_zero = vector("list", length(allsteps.marg))
            names(mwbsfs) = names(mwbsfs_zero) = paste0("step-",allsteps.marg)
            for(ii in 1:length(allsteps.marg)){
                numSteps = allsteps.marg[ii]
                obj = wbsfs(y, numSteps=numSteps, numIntervals=length(y))
                obj = addpv(obj, sigma=sigma, type="rand", mn=mn, locs=locs)
                mwbsfs[[ii]] = obj$pvs
                mwbsfs_zero[[ii]] = (obj$means==0)
            }
            results$mwbsfs = mwbsfs
            results$mwbsfs_zero = mwbsfs_zero
            results$mwbsfs_cps = obj$cp * obj$cp.sign 

        }, error=function(err){ print('error occurred during marginalized wbsfs')})}
        
        ## Plain CBS inference (Non-marginalized)
        if(any(type=="cbsfs")){tryCatch({
            obj = cbsfs(y, numSteps=max(allsteps.cbs))
            poly.max = polyhedra(obj, numSteps=max(allsteps.cbs), record.nrows=TRUE)
            res = plain_inf_multistep(obj, allsteps.cbs, poly.max, mn, sigma, locs=locs)
            results$cbsfs = res$pvs.by.step
            results$cbsfs_zero = res$zeros.by.step
            results$cbsfs_cps = obj$cp * obj$cp.sign 
        }, error=function(err){ print('error occurred during plain cbsfs')})}

        ## Noisy CBS inference
        if(any(type=="ncbsfs")){tryCatch({
            obj = cbsfs(y=y, y.addnoise=y.addnoise, numSteps=max(allsteps.cbs),
                        sigma.add=sigma.add)

            poly.max = polyhedra(obj, numSteps=max(allsteps.cbs), record.nrows=TRUE)
            res = plain_inf_multistep(obj, allsteps.cbs, poly.max, mn, sigma,
                                      shift=y.addnoise, locs=locs)
            results$ncbsfs = res$pvs.by.step
            results$ncbsfs_zero = res$zeros.by.step
            results$ncbsfs_cps = obj$cp * obj$cp.sign 
        }, error=function(err){ print('error occurred during ncbsfs')})}
        
        ## Marginalized noisy CBS inference
        if(any(type=="mcbsfs")){tryCatch({
            mcbsfs = mcbsfs_zero = vector("list", length(allsteps.cbs.marg))
            names(mcbsfs) = names(mcbsfs_zero) = paste0("step-",allsteps.cbs.marg)
            for(ii in 1:length(allsteps.cbs.marg)){
                numSteps = allsteps.cbs.marg[ii]
                obj = cbsfs(y, numSteps=numSteps, sigma.add=sigma.add)
                obj = addpv(obj, sigma=sigma, sigma.add=sigma.add, type="addnoise", mn=mn,
                            locs=locs)
                mcbsfs[[ii]] = obj$pvs
                mcbsfs_zero[[ii]] = (obj$means==0)
            }
            results$mcbsfs = mcbsfs
            results$mcbsfs_zero = mcbsfs_zero
            results$mcbsfs_cps = obj$cp * obj$cp.sign 
        }, error=function(err){ print('error occurred during mcbsfs')})}
        
        ## Plain FL inference 
        if(any(type=="fl")){tryCatch({
            obj = fl(y, numSteps=max.numSteps)
            poly.max = polyhedra(obj, numSteps=max.numSteps)
            res = plain_inf_multistep(obj, allsteps, poly.max, mn, sigma, locs=locs)
            results$fl = res$pvs.by.step
            results$fl_zero = res$zeros.by.step
            results$fl_cps = obj$cp * obj$cp.sign 
        }, error=function(err){ print('error occurred during plain fl')})}

        ## ## Plain FL inference with improvements: under construction!!
        ## if(any(type=="fl")){tryCatch({
        ##     obj = fl(y, numSteps=max.numSteps)
        ##     poly.max = polyhedra(obj, numSteps=max.numSteps)
        ##     res = plain_inf_multistep(obj, allsteps, poly.max, mn, sigma, locs=locs)
        ##     results$fl = res$pvs.by.step
        ##     results$fl_zero = res$zeros.by.step
        ##     results$fl_cps = obj$cp * obj$cp.sign 
        ## }, error=function(err){ print('error occurred during plain fl')})}
        
        ## Noisy FL inference (Non-marginalized)
        if(any(type=="nfl")){tryCatch({
            obj = fl(y=y, y.addnoise=y.addnoise, numSteps=max.numSteps, sigma.add=sigma.add)
            poly.max = polyhedra(obj, numSteps=max.numSteps)
            res = plain_inf_multistep(obj, allsteps, poly.max, mn, sigma,
                                      shift=y.addnoise, locs=locs)
            results$nfl = res$pvs.by.step
            results$nfl_zero = res$zeros.by.step
            results$nfl_cps = obj$cp * obj$cp.sign 
        }, error=function(err){ print('error occurred during noisy fl')})}

        ## Marginalized noisy FL inference
        if(any(type=="mfl")){tryCatch({
            mfl = mfl_zero = vector("list", length(allsteps.marg))
            names(mfl) = names(mfl_zero) = paste0("step-",allsteps.marg)
            for(ii in 1:length(allsteps.marg)){
                numSteps = allsteps.marg[ii]
                obj = fl(y, numSteps=numSteps, sigma.add=sigma.add)
                obj = addpv(obj, sigma=sigma, sigma.add=sigma.add, type="addnoise", mn=mn,
                            locs=locs, max.numIS=max.numIS.fl,
                            min.num.things=min.num.things.fl)
                mfl[[ii]] = obj$pvs
                mfl_zero[[ii]] = (obj$means==0)
            }
            results$mfl = mfl
            results$mfl_zero = mfl_zero
            results$mfl_cps = obj$cp * obj$cp.sign 

        }, error=function(err){ print(paste0('error occurred during noisy fl isim=', isim)) })}

        ## Marginalized noisy decluttered FL inference
        if(any(type=="mdfl")){tryCatch({
            mdfl = mdfl_zero = vector("list", length(allsteps.marg))
            names(mdfl) = names(mdfl_zero) = paste0("step-",allsteps.marg)
            for(ii in 1:length(allsteps.marg)){
                numSteps = allsteps.marg[ii]
                obj = fl(y, numSteps=numSteps, sigma.add=sigma.add)
                obj = addpv(obj, sigma=sigma, sigma.add=sigma.add, type="addnoise", mn=mn,
                            locs=locs, max.numIS=max.numIS.fl, declutter=TRUE,
                            min.num.things=min.num.things.fl)
                mdfl[[ii]] = obj$pvs
                mdfl_zero[[ii]] = (obj$means==0)
            }
            results$mdfl = mdfl
            results$mdfl_zero = mdfl_zero
            results$mdfl_cps = obj$cp * obj$cp.sign 

        }, error=function(err){ print(paste0('error occurred during noisy fl isim=', isim)) })}

        ## "D"ecluttered plain FL inference
        if(any(type=="dfl")){tryCatch({
            obj = fl(y, numSteps=max.numSteps)
            poly.max = polyhedra(obj, numSteps=max.numSteps)
            res = plain_inf_multistep(obj, allsteps, poly.max, mn, sigma,
                                      locs=locs, declutter=TRUE,
                                      how.close=how.close)
            results$dfl = res$pvs.by.step
            results$dfl_zero = res$zeros.by.step
            results$dfl_cps = obj$cp * obj$cp.sign 

        }, error=function(err){ print(paste0('error occurred during decluttered fl isim=', isim)) })}
        
        return(results)
        }

    ## Run the actual simulations.
    start.time = Sys.time()
    results.list = Mclapply(nsim, onesim, mc.cores, Sys.time(), mc.preschedule=mc.preschedule)

    ## Save results
    ## if(is.null(filename)){ filename = paste0("compare-power-lev-",
    ##                                          myfractions(lev), "-ichunk-",
    ##                                          ichunk, "-multistep-plain.Rdata")}
    save(results.list, file=file.path(outputdir, filename))
}
        

convert_pv_to_two_sided_test <- function(pv){
    return(2 * min(pv, 1-pv))
}

##' Helper for dosim(), in power comparison simulation. Collect information
##' (plain saturated p-values, zero-ness of mean) for all steps in |allsteps|.
##' @param allsteps All algorithm steps to test.
plain_inf_multistep <- function(obj, allsteps, poly.max, mn, sigma, shift=NULL,
                                locs=1:length(mn), declutter=FALSE, how.close=2){
    pvs.by.step = zeros.by.step = list()
    for(istep in 1:length(allsteps)){
        numSteps = allsteps[istep]
 
         if(declutter){
             out = declutter(obj$cp[1:numSteps],
                             obj$cp.sign[1:numSteps],
                             how.close=how.close)
             cp.clean = unlist(out$cp)
             cp.sign.clean = unlist(out$cp.sign)
 
             ## Handle the two-sided test cases.
             which.two.sided = which(is.na(cp.sign.clean))
             cp.sign.clean[which.two.sided] = 1
 
             ## Make contrasts
             vlist = make_all_segment_contrasts_from_cp(cp = cp.clean,
                                                        cp.sign=cp.sign.clean,
                                                        n=length(obj$y))
         } else {
             vlist = make_all_segment_contrasts(obj, numSteps=numSteps)
         }
         vlist = filter_vlist(vlist, locs=locs)
 
        ## Store inference results
        pvs = poly_pval2_from_vlist(y=obj$y.orig,
                                    poly=snapshot(poly.max, numSteps),
                                    vlist=vlist, sigma=sigma, shift=shift)
        if(declutter){
            if(length(which.two.sided)>0){
                pvs[which.two.sided] = sapply(pvs[which.two.sided], convert_pv_to_two_sided_test)
            }
        }
        pvs.by.step[[istep]] = pvs
        zeros.by.step[[istep]] = sapply(vlist, function(v) sum(v*mn)==0)
    }

    names(pvs.by.step) = names(zeros.by.step) = paste0("step-", allsteps)
    return(list(pvs.by.step=pvs.by.step, zeros.by.step=zeros.by.step))
}

