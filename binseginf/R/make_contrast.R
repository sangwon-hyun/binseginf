##' Helper function for making plain segment contrasts from a bsfs/wbsfs/cbsfs
##' object or path object (from genlassoinf package). Essentially, the
##' requirements are that \code{obj} has a nonempty component called \code{cp}
##' and \code{cp.sign} and \code{y}.
##'
##' @param obj Result from running one of: \code{bsfs(), bsft(), wbsfs(),
##'     wbsft(), cbsfs()}.
##'
##' @export
make_all_segment_contrasts <- function(obj, numSteps=obj$numSteps){

    ## Basic checks
    if(length(obj$cp)==0) stop("No detected changepoints!")
    if(all(is.na(obj$cp)))stop("No detected changepoints!")
    assert_that(!is.null(obj$y))

    cp = obj$cp
    cp.sign = obj$cp.sign

    ## Handle CBS specially because it detects 2 at a time, sometimes with
    ## overlaps
    if("cbsfs" %in% class(obj) ){
        numSteps = numSteps * 2
        cp = obj$cp.all
        cp.sign = sign(obj$cp.sign.all)
    }

    cp = cp[1:numSteps]
    cp.sign = cp.sign[1:numSteps]
    if(any(is.na(cp))){
        na.cp = which(is.na(cp))
        cp = cp[-na.cp]
        cp.sign = cp.sign[-na.cp]
    }

    return(make_all_segment_contrasts_from_cp(cp, cp.sign, length(obj$y)))
}


##' Take |cp| and |cp.sign| and make all segment contrast vectors.
##' @param cp An integer vector of changepoints.
##' @param cp.sign An integer vector of changepoint signs.
##' @param n Data length.
##' @return A list of test contrast vectors, whose names are the signed
##'     changepoint locations.
##' @export
make_all_segment_contrasts_from_cp <- function(cp, cp.sign, n, scaletype = c("segmentmean", "unitnorm")){

    ## Basic checks
    scaletype = match.arg(scaletype)

    ## Augment the changepoint set for convenience
    ord = order(cp)
    cp_aug = c(0,cp[ord], n)
    sn_aug = c(NA,cp.sign[ord],NA)
    dlist = list()

    ## Make each contrast
    for(ii in (2:(length(cp_aug)-1))){
        d = rep(0,n)
        ind1 = (cp_aug[ii-1]+1):cp_aug[ii] ## 1 to 3, 4 to 9
        ind2 = (cp_aug[ii]+1):cp_aug[ii+1]
        d[ind1] = -1/length(ind1)
        d[ind2] = 1/length(ind2)
        if(scaletype == "unitnorm"){
            d = d/sqrt(sum(d*d))
        } else if (scaletype == "segmentmean"){
        } else {
            stop("scaletype not written yet!")
        }
        dlist[[ii-1]] = d * sn_aug[ii]
    }
    names(dlist) = (cp * cp.sign)[ord]
    return(dlist)
}


##' Makes segment contrasts with endpoints as the winning intervals' start/ends
##' @param wbsfs_obj An object of class \code{wbsfs}.
##' @param cps if not \code{NULL}, only form contrasts from changepoints in this
##'     set.
##' @return A list of test contrast vectors, whose names are the signed
##'     changepoint locations.
##' @export
make_all_segment_contrasts_from_wbs <- function(wbsfs_obj, cps=NULL, scaletype = c("segmentmean", "unitnorm")){

    n = length(wbsfs_obj$y)
    cp = wbsfs_obj$cp
    cp.sign = wbsfs_obj$cp.sign
    winning.intervals = lapply(1:nrow(wbsfs_obj$results), function(myrow)wbsfs_obj$results[myrow, c("max.s", "max.e")])

    scaletype = match.arg(scaletype)

    ## Make each contrast
    dlist = list()
    for(ii in 1:length(cp)){
        d = rep(0,n)
        ind1 = (winning.intervals[[ii]]["max.s"]+1):cp[ii]
        ind2 = (cp[ii]+1):winning.intervals[[ii]]["max.e"]
        d[ind1] = -1/length(ind1)
        d[ind2] = 1/length(ind2)
        dlist[[ii]] = d * cp.sign[ii]

        if(scaletype == "unitnorm"){
            d = d/sqrt(sum(d*d))
        } else if (scaletype == "segmentmean"){
            ## Do nothing
        } else {
            stop("scaletype not written yet!")
        }
        assert_that(length(d)==n)
    }
    names(dlist) = (cp * cp.sign)

    ## Only return the requested ones
    if(!is.null(cps)){
        dlist = dlist[which(abs(as.numeric(names(dlist))) %in% cps)]
    }
    return(dlist)
}
