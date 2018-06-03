##' Conducts post-parseselection inference procedure for additive noise binary
##' segmentation inference -- segment test inference for all detected
##' changepoints -- for a given n-lengthed data vector y, with or without IC
##' stopping.
##' @param added.noise is a manually inputted additive noise vector. This must
##'     be generated from i.i.d. Gaussian noise with \code{sigma} standard
##'     deviation. A rough check is in place, but really just trusting the user
##'     at this point.
##' @param y data vector
##' @param consec The number of rises you'd like as a stopping rule, in terms of
##'     BIC.
##' @param max.numSteps The (maximum) number of algorithm steps you'd like to
##'     do, prior to BIC. If you say \code{icstop=FALSE}, then this is the
##'     actual, fixed number of algorithm steps used.
##' @param icstop Whether to use IC stopping. Defaults to \code{TRUE}.
##' @param verbose Whether to print progress of
##' @param mc.cores Whether to run importance sampler on multiple cores. Not
##'     recommended for normal user.
##' @param start.time Optional argument, for timing things when
##'     \code{verbose=TRUE}
##' @param postprocess \code{TRUE} if you'd like to centroid cluster the
##'     detected changepoints.
##' @param how.close How close would you like points to be in each cluster, for
##'     the centroid clustering, to be.
##' @param locs Which locations to test. Defaults to \code{1:length(y)}.
##' @param return An object of class \code{bsFs} with p-values.
inference_bsFs <- function(y=y, max.numSteps=20, consec=2, sigma, icstop=TRUE,
                           postprocess=TRUE, locs=1:length(y), numIS=100,
                           sigma.add = 0.2, bits=50,
                           inference.type=c("rows", "pre-multiply"),
                           numIntervals=length(y),
                           max.numIS=2000, verbose=FALSE, min.num.things=10,
                           added.noise=NULL,
                           mc.cores=1,
                           start.time=NULL,
                           how.close=5,
                           mn = NULL,
                           sim.options = list(retain.only.null.cases=FALSE),
                           bootstrap.inds=NULL){
    ## Basic checks
    inference.type = match.arg(inference.type)
    if(!is.null(added.noise)){
        if( abs(sd(added.noise) - sigma.add) > sigma.add/2){
            stop("Your added noise doesn't match sigma.add well.")
        }
    }

    ## Fit model and get IC information
    n = length(y)
    new.noise = rnorm(length(y),0,sigma.add)
    y.fudged = y + new.noise
    h.fudged = binSeg_fixedSteps(y.fudged, numSteps=max.numSteps)
    h.fudged$y.orig = y
    if(icstop){
        ic_obj = get_ic(h.fudged$cp, h.fudged$y, consec=consec,
                        sigma=sigma+sigma.add, type="bic")
        stoptime = ic_obj$stoptime
        if(ic_obj$flag!="normal"){
            ## print(ic_obj$flag)
            warning(paste0("IC stopping resulted in: ", ic_obj$flag))
            return(NA)
        }
    } else {
        stoptime = max.numSteps
        ic_obj=NULL
    }

    ## Collect stopped model and postprocess
    cp = h.fudged$cp[1:stoptime]
    cp.sign = h.fudged$cp.sign[1:stoptime]
    if(postprocess){
        cpobj = declutter_new(cp, cp.sign, how.close=how.close)
        cp = abs(cpobj)
        cp.sign = sign(cpobj)
    }

    ## It is important to form the contrasts /after/ decluttering
    ## Retain only the changepoints we want results from:
    vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n)
    vlist <- filter_vlist(vlist, locs)
    if(length(vlist)==0) return(list(pvs=c(), null.true=c()))

    ## (kind of) temporary addition
    if(sim.options$retain.only.null.cases){
        truth = sapply(vlist, function(myv){
            sum(myv*mn)
        })
        tol=1E-20
        retain = which(abs(truth)>tol)
        vlist = vlist[retain]
    }
    if(length(vlist)==0) return(NULL)

    ## Do noise-added inference
    parts.list = list()
    pvs = rep(NA, length(vlist))
    names(pvs) = names(vlist)
    for(iv in 1:length(vlist)){
        v = vlist[[iv]]
        if(verbose) printprogress(iv, length(vlist), "segment tests", fill=TRUE)
        result = randomize_addnoise(y=y, v=v, sigma=sigma, numIS=numIS,
                                    sigma.add=sigma.add, orig.fudged.obj = h.fudged,
                                    numSteps = stoptime+consec,
                                    ic.poly = ic_obj$poly, bits=bits,
                                    inference.type=inference.type,
                                    max.numIS=max.numIS, verbose=verbose,
                                    mc.cores= mc.cores,
                                    start.time=start.time,
                                    min.num.things=min.num.things)
        if(verbose) cat(fill=TRUE)
        parts.list[[iv]] = result$parts.so.far
        pvs[iv] = result$pv
    }

    ## Add some information to the original fudged object:
    h.fudged$pvs = pvs
    h.fudged$anyleft = TRUE
    h.fudged$vlist = vlist
    h.fudged$parts.list = parts.list
    h.fudged$stoptime = stoptime
    if(!is.null(mn)){
        h.fudged$resids = y - mn
        h.fudged$mn = mn
    }
    h.fudged$bootstrap.inds = bootstrap.inds
    return(h.fudged)
}



##' Does a single randomized wbs (rwbs) inference for a given y
##' @param numIntervals Defaults to \code{length(y)}, or if |intervals| is
##'     provided, the length of that.
##' @param intervals An object of class |intervals|.
##' @param postprocess \code{TRUE} if you'd like to centroid cluster the
##'     detected changepoints.
##' @param howclose How close would you like points to be in each cluster, for
##'     the centroid clustering, to be.
##' @return List of simulation results (p-values, locations tested, etc.)
inference_wbs <- function(y=y, max.numSteps=20, numIntervals=length(y),
                          intervals=NULL, consec=2,
                          sigma, postprocess=TRUE, how.close = 5,
                          better.segment=FALSE,
                          locs=1:length(y), numIS=100,
                          inference.type=inference.type,
                          bits=1000,
                          max.numIS=2000,
                          verbose=FALSE, mc.cores=1,
                          min.num.things=30){

    ## Basic checks
    if(!is.null(intervals)){
        if(is.null(intervals)){stop("Provide either |numIntervals| or |intervals|.")}
        numIntervals = intervals$numIntervals
    }

    n=length(y)

    ## Fit initial WBS for a generous number of steps
    if(is.null(intervals) & !is.null(numIntervals)){
        intervals = intervals(numIntervals=numIntervals, n=n)
    }
    g = wildBinSeg_fixedSteps(y, intervals=intervals, numSteps=max.numSteps,
                              inference.type='none')
    cumsum.y = cumsum(y)

    ## Collect the IC information and polyhedron
    ic_obj = get_ic(g$cp, g$y, consec=consec, sigma=sigma, type="bic")
    ic_poly = ic_obj$poly
    stoptime  = ic_obj$stoptime

    if(verbose) cat("stoptime is", stoptime, fill=TRUE)

    ## Check for flag
    if(ic_obj$flag!="normal" ){
        return(NULL)
    }

    ## Extract changepoints from stopped model and declutter
    cp = g$cp[1:stoptime]
    cp.sign = g$cp.sign[1:stoptime]

    if(postprocess){
        cpobj = declutter_new(cp, cp.sign, how.close=how.close)
        cp = abs(cpobj)
        cp.sign = sign(cpobj)
    }

    ## Form contrasts
    if(better.segment){
        vlist <- make_all_segment_contrasts_from_wbs(wbs_obj=g, cps=cp)
    } else {
        vlist <- make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=length(y))
    }

    ## Retain only the changepoints we want results from:
    retain = which((abs(as.numeric(names(vlist))) %in% locs))
    if(length(retain)==0) return(list(pvs=c(), null.true=c()))
    vlist = vlist[retain]

    ## Calculate the p-values
    pvs = sapply(1:length(vlist), function(iv){
        printprogress(iv, length(vlist), type = "tests")
        v = vlist[[iv]]
        cumsum.v = cumsum(v)
        pv = suppressWarnings(randomize_wbsfs(v=v, winning.wbs.obj=g,
                                              sigma=sigma,
                                              numIS=numIS,
                                              inference.type=inference.type,
                                              cumsum.y=cumsum.y,
                                              cumsum.v=cumsum.v,
                                              stop.time=stoptime+consec,
                                              ic.poly=ic_poly,
                                              bits=bits,
                                              max.numIS=max.numIS,
                                              mc.cores=mc.cores, min.num.things=min.num.things))
        return(pv)})
    names(pvs) = names(vlist)
    return(list(pvs=pvs, locs.all=cp*cp.sign, locs.retained=as.numeric(names(pvs)),
                vlist=vlist) )
}
