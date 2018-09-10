##' Collecting the basic information to be shared across importance sampling
##' replicates.
get_basics <- function(y, poly, v, sigma, bits){

    ## Form everything
    stopifnot(sum(v*v)==1)
    vy = sum(v*y)
    vv = sum(v^2)
    sd = sigma * sqrt(vv)
    Gv = poly$gamma %*% v
    Gv[which(abs(Gv)<1E-15)] = 0
    rho = Gv / vv
    Gy = poly$gamma %*% y

    return(list(vv=vv, sd=sd, Gy=Gy, rho=rho,
                rho.times.vy = rho*vy, u = poly$u, vy=vy))
}


##' Reusing the basic information that is shared across all importance sampling
##' replicates.
partition_TG2 <- function(obj.basic, Gw, bits=2000){

    ## Load obj.basic that contains shared information
    new.u = obj.basic$u - Gw
    vec = (new.u - obj.basic$Gy + obj.basic$rho.times.vy) / obj.basic$rho
    vlo = suppressWarnings(max(vec[obj.basic$rho>0]))
    vup = suppressWarnings(min(vec[obj.basic$rho<0]))
    vy = max(min(obj.basic$vy, vup), vlo) ## This is required.

    ## Myabe there is a copying cost?

    ## Calculate a,b,z for TG = (F(b)-F(z))/(F(b)-F(a))
    z = Rmpfr::mpfr(vy/obj.basic$sd, precBits=bits) ## Even this is unnecessary, I think.
    a = Rmpfr::mpfr(vlo/obj.basic$sd, precBits=bits)
    b = Rmpfr::mpfr(vup/obj.basic$sd, precBits=bits)

    numer = as.numeric(Rmpfr::pnorm(b)-Rmpfr::pnorm(z))
    denom = as.numeric(Rmpfr::pnorm(b)-Rmpfr::pnorm(a))
    pv = as.numeric(numer/denom)
    return(list(denom=denom, numer=numer, pv=pv, vlo=vlo, vy=vy, vup=vup))
}

reuse_basics = partition_TG2 


partition_TG3 <- function(u, Gy, rho.times.vy, rho, vy, Gw, sd, bits=2000){

    ## Load obj.basic that contains shared information
    new.u = u - Gw
    vec = (new.u - Gy + rho.times.vy) / rho
    vlo = suppressWarnings(max(vec[rho>0]))
    vup = suppressWarnings(min(vec[rho<0]))
    vy = max(min(vy, vup), vlo) ## This is required.

    ## Myabe there is a copying cost?

    ## Calculate a,b,z for TG = (F(b)-F(z))/(F(b)-F(a))
    z = Rmpfr::mpfr(vy/sd, precBits=bits) ## Even this is unnecessary, I think.
    a = Rmpfr::mpfr(vlo/sd, precBits=bits)
    b = Rmpfr::mpfr(vup/sd, precBits=bits)

    numer = as.numeric(Rmpfr::pnorm(b)-Rmpfr::pnorm(z))
    denom = as.numeric(Rmpfr::pnorm(b)-Rmpfr::pnorm(a))
    pv = as.numeric(numer/denom)
    return(list(denom=denom, numer=numer, pv=pv, vlo=vlo, vy=vy, vup=vup))
}



## partition_TG <- function(y, poly, v, sigma, nullcontrast=0, bits=50, reduce,
##                          correct.ends=FALSE, shift=NULL, ic.poly=NULL, warn=TRUE, obj.basic){

##     ## This is the part that is made in advance
##     ## vy = sum(v*y)
##     vy = obj.basic$vy
##     ## vv = sum(v^2)
##     vv = obj.basic$v
##     sd = sigma*sqrt(vv)
##     poly$u = poly$u - poly$gamma%*%shift
##     pvobj <- poly.pval2(y, poly, v, sigma, correct.ends=correct.ends)
##     vup = pvobj$vup
##     vlo = pvobj$vlo
##     vy = max(min(vy, vup),vlo)

##     ## Calculate a,b,z for TG = (F(b)-F(z))/(F(b)-F(a))
##     z = Rmpfr::mpfr(vy/sd, precBits=bits)
##     a = Rmpfr::mpfr(vlo/sd, precBits=bits)
##     b = Rmpfr::mpfr(vup/sd, precBits=bits)
##     if(!(a<=z &  z<=b) & warn){
##         warning("F(vlo)<vy<F(vup) was violated, in partition_TG()!")
##     }

##     ## Separately store and return num&denom of TG
##     numer = as.numeric(Rmpfr::pnorm(b)-Rmpfr::pnorm(z))
##     denom = as.numeric(Rmpfr::pnorm(b)-Rmpfr::pnorm(a))

##     ## Form p-value as well.
##     pv = as.numeric(numer/denom)

##     return(list(denom=denom, numer=numer, pv=pv, vlo=vlo, vy=vy, vup=vup))
## }



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
randomize_addnoise2 <- function(y, sigma, sigma.add, v, orig.fudged.poly=NULL,
                                numSteps=NA, numIS=100, bits=1000,
                                orig.fudged.obj = NULL, ic.poly=NULL,
                                max.numIS=2000,
                                inference.type = c("rows", "pre-multiply"),
                                verbose=FALSE,
                                min.num.things=10,
                                mc.cores=1,
                                start.time=NULL,
                                warn=FALSE,
                                stable.thresh=1E-3
                                ){

    ## New: Get many fudged TG statistics.
    inference.type = match.arg(inference.type)
    if(sigma.add==0) numIS=max.numIS=1

    if(inference.type!="rows"){ stop()}

    ## If applicable, append IC poly to the original poly
    poly = orig.fudged.poly

    ## Form shared information
    ## obj.basic = get_basics(y=y, poly=poly, v=v, sigma=sigma)

    ## Form everything
    stopifnot(sum(v*v)==1)
    vy = sum(v*y)
    vv = sum(v^2)
    sd = sigma * sqrt(vv)
    Gv = poly$gamma %*% v
    Gv[which(abs(Gv)<1E-15)] = 0
    rho = Gv / vv
    Gy = poly$gamma %*% y
    u = poly$u
    rho.times.vy = rho*vy


    ## Helper function
    one_IS_addnoise2 = function(isim, numIS.cumulative){
        if(verbose) {printprogress(isim+numIS.cumulative, numIS+numIS.cumulative,
                                  paste0("importance sampling replicate (with ", things, " valid draws so far):"),
                                  start.time = start.time)}

        new.noise = rnorm(length(y),0,sigma.add)

        ## obj.new = partition_TG(y=y, poly=poly, shift=new.noise,
        ##                            v=v, sigma=sqrt(sigma^2), bits=bits)
        Gw = poly$gamma %*% new.noise
        bits = 2000
        ## obj.new = reuse_basics(obj.basic, Gw=Gw, bits=bits)
        obj.new = partition_TG3(u, Gy, rho.times.vy, rho, vy, Gw, sd, bits=bits)

        ## Handle boundary cases
        pv.new = obj.new$pv
        weight.new = obj.new$denom
        if(is.nan(pv.new)){ ## This happens when the TG components are
                            ## numerically incalculable because vlo and vty are
                            ## too high.
            emptyrow = rbind(c(0,0,NA,NA,NA,NA)) ## Setting both pv and
                                                      ## weight to zero so it
                                                      ## doesn't count.
            colnames(emptyrow)= c("pv", "weight", "vlo", "vty", "vup", "sigma")
            return(emptyrow)
        }
        if(pv.new > 1 | pv.new < 0)  browser() ## Not sure why this would happen, but anyway!
        if(weight.new < 0 | weight.new > 1){
            weight.new = 0 ## Nomass problem is to be caught here.  But nomass
            ## problem probably doesn't happen for additive noise.
        }
        info = cbind(pv=pv.new, weight=weight.new, vlo=obj.new$vlo,
                     vty=obj.new$vy, vup=obj.new$vup, sigma=sigma)
        return(info)
    }

    ## one_IS_addnoise = function(isim, numIS.cumulative, investigate){
    ##     if(verbose) {printprogress(isim+numIS.cumulative, numIS+numIS.cumulative,
    ##                               paste0("importance sampling replicate (with ", things, " valid draws so far):"),
    ##                               start.time = start.time)}

    ##     new.noise = rnorm(length(y),0,sigma.add)

    ##     obj.new = partition_TG(y=y, poly=poly, shift=new.noise,
    ##                            v=v, sigma=sqrt(sigma^2), bits=bits, obj.basic=obj.basic)

    ##     ## Handle boundary cases
    ##     pv.new = obj.new$pv
    ##     weight.new = obj.new$denom
    ##     if(is.nan(pv.new)){ ## This happens when the TG components are
    ##                         ## numerically incalculable because vlo and vty are
    ##                         ## too high.
    ##         emptyrow = rbind(c(0,0,NA,NA,NA,NA)) ## Setting both pv and
    ##                                                   ## weight to zero so it
    ##                                                   ## doesn't count.
    ##         colnames(emptyrow)= c("pv", "weight", "vlo", "vty", "vup", "sigma")
    ##         return(emptyrow)
    ##     }
    ##     if(pv.new > 1 | pv.new < 0)  browser() ## Not sure why this would happen, but anyway!
    ##     if(weight.new < 0 | weight.new > 1){
    ##         weight.new = 0 ## Nomass problem is to be caught here.  But nomass
    ##         ## problem probably doesn't happen for additive noise.
    ##     }
    ##     info = cbind(pv=pv.new, weight=weight.new, vlo=obj.new$vlo,
    ##                  vty=obj.new$vy, vup=obj.new$vup, sigma=sigma)
    ##     return(info)
    ## }



    ## Importance-sample until you have some amount of variation.
    parts.so.far = cbind(c(Inf, Inf, Inf, Inf, Inf, Inf))[,-1, drop=FALSE]
    rownames(parts.so.far) = c("pv", "weight", "vlo", "vty", "vup", "sigma")
    numIS.cumulative = things = 0
    pv.so.far = c(-1)
    pv_from_parts <- function(parts.so.far){
        sum(unlist(Map('*', parts.so.far["pv",], parts.so.far["weight",])))/
                                                sum(unlist(parts.so.far["weight",]))}
    done = FALSE
    while(!done){
        parts = mcmapply(one_IS_addnoise2, 1:numIS,
                         numIS.cumulative=numIS.cumulative , mc.cores=mc.cores)
        ## onesim <- function(isim){one_IS_addnoise(isim, numIS.cumulative, obj.basic)}
        ## parts = mcmapply(onesim, 1:numIS, mc.cores=mc.cores)

        ## Combine the new parts with the prexisting.
        parts.so.far = cbind(parts.so.far, parts)

        ## Handling the problem of p-value being NaN/0/1
        things = sum((parts.so.far["weight",] > 0) &
                     (parts.so.far["pv",] != 1)) ## &
                     ## (parts.so.far["pv",] != 0)) ## Temporarily excluded b/c it seems extraneous.
        pv.latest = pv_from_parts(parts.so.far)

        enough.things = (things >= min.num.things)
        numIS.cumulative = numIS.cumulative + numIS
        reached.limit = numIS.cumulative > max.numIS
        stable.enough = (stable(pv.latest, pv.so.far, stable.thresh) & things > 10) 

        ## Not used now: Maybe small p-values should just be accepted, for practical reasons
        ## pvtol = 1E-15 pv.is.really.small = (latest.pv < pvtol)


        if(reached.limit | enough.things | sigma.add == 0){ done = TRUE }
        pv.so.far = c(pv.so.far, pv.latest)
    }

    ## Calculate randomized TG statistic.
    pv = pv_from_parts(parts.so.far)

    return(list(things=things, min.num.things=min.num.things,
                numIS.cumulative=numIS.cumulative, parts.so.far=parts.so.far,
                pv=pv, sigma=sigma, v=v, parts=parts))
}

