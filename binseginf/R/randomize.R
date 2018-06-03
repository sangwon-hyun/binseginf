##' Conduct noise-added saturated model inference, for fused lasso or binary
##' segmentation (But also applies to sequential changepoint methods that
##' creates a valid polyhedron and has $cp and $cp.sign, like fused lasso from
##' the genlassoinf::dualpathSvd2()). Implements an importance sampling scheme.
##' @param y Data
##' @param sigma Gaussian noise standard deviation in data around mean.
##' @param sigma.add Amount of Gaussian noise to add. NOTE if equal to zero,
##'     then reverts back to a non-marginalized noise addition.
##' @param v Test contrast vector. This is assumed to be fixed in the model
##'     selection event.
##' @param orig.fudged.poly Original polyhedron object (of class
##'     \code{polyhedra}) from adding a single noise.
##' @param numSteps Number of steps originally taken. Defaults to \code{NA}.
##' @param numIS Number of importance sampling replicates originally
##'     desired. The importance sampling continues until there are
##'     \code{min.num.things} number of valid importance sampling draws.
##' @param max.numIS Maximum number of importance sampling replicates to take.
##' @param inference.type If equal to "rows", then the calculation of TG
##'     statistics is done without modification. If equal to "pre-multiply",
##'     then internally Gy = \eqn{\Gamma y} and Gv = \eqn{\Gamma v}. This serves
##'     two purposes: first, if each polyhedra is too big in the first place,
##'     this helps. Secondly, it speeds up calculation considerably
##'     (poly_pval_from_inner_products() does this work, given Gv and Gy.)
##' @param ic.poly Polyhedron of the (2-rise) stopping event.
##' @param mc.cores Number of cores to use for importance sampling. Defaults to
##'     \code{1}.
##' @param start.time Start time of running time; usually output from
##'     \code{Sys.time()}.
##' @param verbose If \code{TRUE}, the importance sampling progress
##' @export
randomize_addnoise <- function(y, sigma, sigma.add, v, orig.fudged.poly=NULL,
                               numSteps=NA, numIS, bits=1000,
                               orig.fudged.obj = NULL, ic.poly=NULL,
                               max.numIS=2000,
                               inference.type = c("rows", "pre-multiply"),
                               verbose=FALSE,
                               min.num.things=10,
                               mc.cores=1,
                               start.time=NULL
                               ){

    ## New: Get many fudged TG statistics.
    inference.type = match.arg(inference.type)
    if(sigma.add==0) numIS=max.numIS=1

    ## If applicable, append IC poly to the original poly
    poly = orig.fudged.poly
    if(!is.null(ic.poly)){
        poly$gamma = rbind(poly$gamma, ic.poly$gamma) ## I don't like this rbind..
        poly$u = c(poly$u, ic.poly$u)
    }

    ## Helper function
    one_IS_addnoise = function(isim, numIS.cumulative){
        if(verbose) {printprogress(isim+numIS.cumulative, numIS+numIS.cumulative,
                                  "importance sampling replicate",
                                  start.time = start.time)}


        new.noise = rnorm(length(y),0,sigma.add)

        if(inference.type=="rows"){
            obj.new = partition_TG(y=y, poly=poly, shift=new.noise,
                                   v=v, sigma=sqrt(sigma^2), bits=bits)
        } else if (inference.type=="pre-multiply"){
            premult = polyhedra.bsFs(orig.fudged.obj,
                                     inference.type="pre-multiply",
                                     new.noise=new.noise, v=v,
                                     numSteps=numSteps, y=y)
            ## Append IC stopping to Gy, Gv, Gw
            if(!is.null(ic.poly)){
                ic.Gy = as.numeric(ic.poly$gamma%*%y)
                ic.Gv = as.numeric(ic.poly$gamma%*%v)
                ic.Gw = as.numeric(ic.poly$gamma%*%new.noise)
                premult$Gy = c(premult$Gy, ic.Gy)
                premult$Gv = c(premult$Gv, ic.Gv)
                premult$Gw = c(premult$Gw, ic.Gw)
                premult$u = c(premult$u, ic.poly$u)
            }
            obj.new = poly_pval_from_inner_products(Gy=premult$Gy,
                                                    Gv=premult$Gv, v=v, y=y,
                                                    sigma=sigma,
                                                    u=premult$u - premult$Gw,
                                                    bits=bits)
            pv = obj.new$pv
            if(is.nan(obj.new$pv)) obj.new$pv=0 ## temporary fix ## Need to deal
                                                ## with this.
        } else {
            stop("inference type not recognized!")
        }

        ## Handle boundary cases
        pv.new = obj.new$pv
        weight.new = obj.new$denom
        if(is.nan(pv.new)) return(c(0,0)) ## Actually not calculable
        if(pv.new > 1 | pv.new < 0)  browser() ## Not sure why this would happen, but anyway!
        if(weight.new < 0 | weight.new > 1){
            weight.new=0 ## Nomass problem is to be caught here.
        }
        info = cbind(pv=pv.new, weight=weight.new, vlo=obj.new$vlo,
                     vty=obj.new$vy, vup=obj.new$vup, sigma=sigma)
        return(info)
    }

    ## Importance-sample until you have some amount of variation.
    parts.so.far = cbind(c(Inf, Inf, Inf, Inf, Inf, Inf))[,-1, drop=FALSE]
    rownames(parts.so.far) = c("pv", "weight", "vlo", "vty", "vup", "sigma")
    numIS.cumulative = 0
    done = FALSE
    while(!done){
        parts = mcmapply(one_IS_addnoise, 1:numIS, numIS.cumulative,
                         mc.cores=mc.cores)

        ## Combine the new parts with the prexisting.
        parts.so.far = cbind(parts.so.far, parts)

        ## Handling the problem of p-value being NaN/0/1
        things = sum((parts.so.far["weight",]>0)&
                     (parts.so.far["pv",] != 1) &
                     (parts.so.far["pv",] != 0))

        enough.things = (things >= min.num.things)
        numIS.cumulative = numIS.cumulative + numIS
        reached.limit = numIS.cumulative > max.numIS
        if(reached.limit | enough.things | sigma.add == 0){ done = TRUE }
    }

    ## Calculate randomized TG statistic.
    pv = sum(unlist(Map('*', parts.so.far["pv",], parts.so.far["weight",])))/
        sum(unlist(parts.so.far["weight",]))

    return(list(things=things, min.num.things=min.num.things,
                numIS.cumulative=numIS.cumulative, parts.so.far=parts.so.far,
                pv=pv, sigma=sigma, v=v, parts=parts))
}


##' Conduct randomized inference for WBS. Implements an importance sampling
##' scheme.
##' @param y Data
##' @param sigma Gaussian noise standard deviation in data around mean.
##' @param winning.wbs.obj Winning WBS object of class \code{wbsFs}.
##' @param v Test contrast vector. This is assumed to be fixed in the model
##'     selection event.
##' @param orig.fudged.poly Original polyhedron object (of class
##'     \code{polyhedra}) from adding a single noise.
##' @param numSteps Number of steps originally taken. Defaults to \code{NA}.
##' @param numIS Number of importance sampling replicates originally
##'     desired. The importance sampling continues until there are
##'     \code{min.num.things} number of valid importance sampling draws.
##' @param max.numIS Maximum number of importance sampling replicates to take.
##' @param inference.type If equal to "rows", then the calculation of TG
##'     statistics is done without modification. If equal to "pre-multiply",
##'     then internally Gy = \eqn{\Gamma y} and Gv = \eqn{\Gamma v}. This serves
##'     two purposes: first, if each polyhedra is too big in the first place,
##'     this helps. Secondly, it speeds up calculation considerably
##'     (poly_pval_from_inner_products() does this work, given Gv and Gy.)
##' @param ic.poly Polyhedron of the (2-rise) stopping event.
##' @param mc.cores Number of cores to use for importance sampling. Defaults to
##'     \code{1}.
##' @param start.time Start time of running time; usually output from
##'     \code{Sys.time()}.
##' @param bits Number of precision bits to use for the calculation of the
##'     Gaussian probabilities (regarding Vup and Vlo and vty).
##' @param verbose If \code{TRUE}, the importance sampling progress
##' @export 
randomize_wbsfs <- function(v, winning.wbs.obj, numIS = 100, sigma,
                            inference.type=c("rows", "pre-multiply"),
                            cumsum.y=NULL,cumsum.v=NULL, stop.time=min(winning.wbs.obj$numSteps,
                                                                       length(winning.wbs.obj$cp)),
                            ic.poly=NULL, bits=50, max.numIS=2000,
                            min.num.things=30, verbose=FALSE,
                            mc.cores=1,
                            warn=FALSE,
                            start.time=NULL){

    numIntervals = winning.wbs.obj$numIntervals
    numSteps = winning.wbs.obj$numSteps

    ## Basic checks
    if(inference.type=="pre-multiply" & (is.null(cumsum.y) | is.null(cumsum.v)) ){
        stop("Provide cumulative sums of y and v, if you want to use the pre-multiply option.")
    }

    ## Helper function (bundler) for a single importance sampling replicate
    one_IS_wbs = function(isim, numIS.cumulative){
        if(verbose) printprogress(isim+numIS.cumulative, numIS+numIS.cumulative,
                                  "importance sampling replicate",
                                  start.time = start.time)


        rerun_wbs(v=v, winning.wbs.obj=winning.wbs.obj,
                  numIntervals=numIntervals,
                  numSteps=winning.wbs.obj$numSteps,
                  sigma=sigma,
                  inference.type=inference.type,
                  cumsum.y=cumsum.y,
                  cumsum.v=cumsum.v,
                  stop.time=stop.time,
                  ic.poly=ic.poly,
                  bits=bits,
                  warn=warn)
    }

    ## Actual importance sampling is run here
    done = FALSE
    parts.so.far = cbind(c(Inf,Inf))[,-1,drop=FALSE]
    rownames(parts.so.far) = c("pv", "weight")
    numIS.cumulative=0
    while(!done){
        numIS.cumulative = numIS.cumulative + numIS

        ## Collect parts and combine with preexisting
        parts = mcmapply(one_IS_wbs, 1:numIS , numIS.cumulative, mc.cores=mc.cores)
        parts.so.far = cbind(parts.so.far, parts)

        ## Handling issue of p-value being NaN/0/1
        things = sum(parts.so.far["weight",] > 0)
        enough.things = (things > min.num.things)
        reached.limit = (numIS.cumulative > max.numIS)
        if( reached.limit | enough.things){ done = TRUE }
    }

    ## Calculate p-value
    pv = sum(unlist(Map('*', parts.so.far["pv",], parts.so.far["weight",])))/
        sum(unlist(parts.so.far["weight",]))

    ## Return information from this simulation.
    return(list(things=things, min.num.things=min.num.things, numIS.cumulative=numIS.cumulative,
                parts.so.far=parts.so.far, pv=pv, sigma=sigma, v=v))
}

##' Helper for WBSFT randomization. Rerun WBS to get _new_, single set of denom
##' and numers from a new TG statistic, which is calculated from the _single_
##' new set of halfspaces that characterize the maximization of the original
##' |winning.wbs.obj| but among /different/ set of |numIntervals|-|numSteps|
##' spaces.
##' @param winning.wbs.obj Original contrast. We call it winning because we will
##'     extract only the winning locations and /those/ winners' enclosing
##'     intervals.
##' @param v test contrast
##' @param data frame with two columns ("pv" and "weight") and single row.
rerun_wbs <- function(winning.wbs.obj, v, numIntervals, numSteps, sigma,
                      cumsum.y=NULL,cumsum.v=NULL, inference.type, stop.time=numSteps,
                      ic.poly=NULL, bits=50, warn=FALSE){

    ## Basic checks
    assert_that(is_valid.wbsFs(winning.wbs.obj))

    ## New intervals added onto old winning intervals
    n = length(v)
    stopifnot(n==length(winning.wbs.obj$y))
    winning_se = rbind(winning.wbs.obj$results[1:stop.time, c("max.s", "max.e")])
    colnames(winning_se) = c("s", "e")
    intervals.new = intervals(numIntervals=numIntervals-stop.time, n=n, existing=winning_se)
    intervals.new = add2(intervals=intervals.new,
                         winning.wbs.obj=winning.wbs.obj,
                         stop.time=stop.time)

    ## Create new halfspaces (through |mimic| option)
    if(inference.type=="rows"){
        g.new = wildBinSeg_fixedSteps(y=winning.wbs.obj$y, numSteps= numSteps,
                                      intervals= intervals.new, mimic=TRUE,
                                      wbs.obj=winning.wbs.obj,
                                      inference.type=inference.type)
        poly.new = polyhedra(obj=g.new$gamma, u=g.new$u)

        ## Partition TG to denom and numer
        pvobj = partition_TG(y=winning.wbs.obj$y, poly.new, v=v, sigma=sigma,
                             correct.ends=TRUE, warn=warn)
        pv = pvobj$pv
        if(is.nan(pv)) pv=0 ## temporary fix
        weight = pvobj$denom


    } else {

        ## new way using new function
        g.new = wildBinSeg_fixedSteps(y=winning.wbs.obj$y, numSteps= numSteps,
                                      intervals= intervals.new, mimic=TRUE,
                                      wbs.obj=winning.wbs.obj,
                                      cumsum.y=cumsum.y,
                                      cumsum.v=cumsum.v,
                                      inference.type="pre-multiply",
                                      stop.time=stop.time,
                                      ic.poly=ic.poly,
                                      v=v)

        ## Calculate TG denom and numer directly
        pvobj = poly_pval_from_inner_products(Gy=g.new$Gy, Gv=g.new$Gv, v=v, y=g.new$y,
                                              sigma=sigma, u=g.new$u, bits=bits, warn=warn)
        pv = pvobj$pv
        if(is.nan(pv)) pv=0 ## temporary fix
        weight = pvobj$denom
    }

    ## Special handling so that, if Vup<Vlo, then the weight, which is the prob
    ## along the line trapped in the polyhedron, is manually assigned zero.
    if(weight<0 | weight>1) weight = 0

    info = data.frame(pv=pv,weight=weight)
    return(info)
}

