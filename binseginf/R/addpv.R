## Synopsis: the functions defined here aim to streamline the process of
## producing p-values. The general structure is \code{obj<-bsFs(y, ...)}, then
## \code{obj<-test(obj)}, which appends (adds) the inference results to the same
## object.

##' Function generic for \code{test.OOO()} functions.
##' @export
addpv <- function(obj,...) UseMethod("addpv")

##' Appends the inference results to an object of class |bsFs|.
##' @param obj object of type bsFs
##' @param loc only test locations in \code{loc}.
##' @param type One of \code{ c("plain", "addnoise")}. If equal to
##'     \code{"addnoise"}, then \code{sigma.add} needs to be provided.
##' @param sigma Noise level (standard deviation) of data.
##' @param sigma.add Additive noise. Defaults to NULL, in which case no additive
##'     noise randomization inference is done.
##' @param mn original mean vector.
##' @export
addpv.bsfs <- function(obj, loc=NULL, type=c("plain", "addnoise"), sigma,
                       sigma.add=NULL, declutter=FALSE, mn=NULL, min.num.things=30){

    ## Basic checks
    assert_that(class(obj)=="bsfs")
    assert_that(is.null(obj$pvs))
    type = match.arg(type)
    if(type=="addnoise"){
        assert_that(obj$noisy)
        assert_that(!is.null(obj$sigma.add))
    }

    ## Form the test contrasts
    vlist <- make_all_segment_contrasts(obj)
    vlist <- filter_vlist(vlist, loc)

    ## Obtain p-values
    if(type=="plain"){
        poly.nonfudged = polyhedra(obj)
        pvs = sapply(vlist, function(v){
            pv = poly.pval2(y=obj$y, poly=poly.nonfudged, v=v, sigma=sigma, bits=5000)$pv
        })

    } else if (type=="addnoise") {
        poly.fudged = polyhedra(obj)
        pvs = sapply(vlist, function(v){
            pv = randomize_addnoise(y=obj$y, v=v, sigma=sigma, numIS=10,
                                    sigma.add=sigma.add,
                                    orig.fudged.poly=poly.fudged, bits= 5000,
                                    max.numIS=2000,
                                    min.num.things=min.num.things)$pv})
    } else {
        stop("|type| argument is wrong!")
    }

    ## Add to object
    obj$pvs = pvs
    obj$vlist = vlist
    if(!is.null(mn)){obj$means = sapply(vlist, function(v){ sum(v*mn) })}
    return(obj)
}


##' Appends the inference results to an object of class |bsFs|.
##' @param obj object of type bsFs
##' @param loc only test locations in \code{loc}.
##' @param type One of \code{ c("plain", "addnoise")}. If equal to
##'     \code{"addnoise"}, then \code{sigma.add} needs to be provided.
##' @param sigma Noise level (standard deviation) of data.
##' @param mn original mean vector.
##' @export
addpv.wbsfs <- function(obj, loc=NULL, type=c("plain", "rand"), sigma,
                        declutter=FALSE, mn=NULL, min.num.things=30){

    ## Basic checks
    assert_that(class(obj)=="wbsfs")
    assert_that(is.null(obj$pvs))
    type = match.arg(type)

    ## Form the test contrasts
    vlist <- make_all_segment_contrasts(obj)
    vlist <- filter_vlist(vlist, loc)

    ## Obtain p-value
    if(type=="plain"){
        pvs = sapply(vlist, function(v){
            pv = poly.pval2(y=obj$y, poly=polyhedra(obj=obj$gamma, u=obj$u), v=v, sigma=sigma, bits=5000)$pv
        })
    } else if (type=="rand") {

        vlist <- filter_vlist(make_all_segment_contrasts(obj), loc)

        ## Get the p-values
        pvs = sapply(vlist, function(v){
            pv = randomize_wbsfs(v=v, winning.wbs.obj=obj,
                                 sigma=sigma, numIS=10,
                                 inference.type="pre-multiply",
                                 cumsum.y=cumsum(obj$y),
                                 cumsum.v=cumsum(v), bits=2000,
                                 max.numIS=2000,
                                 min.num.things=min.num.things)$pv
        })

    } else {
        stop("|type| argument is wrong!")
    }

    ## Add to object
    obj$pvs = pvs
    obj$vlist = vlist
    if(!is.null(mn)){obj$means = sapply(vlist, function(v){ sum(v*mn) })}
    return(obj)
}


##' Appends the inference results to an object of class |bsFs|.
##' @param obj object of type bsFs
##' @param loc only test locations in \code{loc}.
##' @param type One of \code{ c("plain", "addnoise")}. If equal to
##'     \code{"addnoise"}, then \code{sigma.add} needs to be provided.
##' @param sigma Noise level (standard deviation) of data.
##' @param sigma.add Additive noise. Defaults to NULL, in which case no additive
##'     noise randomization inference is done.
##' @param mn original mean vector.
##' @export
addpv.cbsfs <- function(obj, loc=NULL, type=c("plain", "addnoise"), sigma,
                       sigma.add=NULL, declutter=FALSE, mn=NULL, min.num.things=30){

    ## Basic checks
    assert_that(class(obj)=="cbsfs")
    assert_that(is.null(obj$pvs))
    type = match.arg(type)
    if(type=="addnoise"){
        assert_that(obj$noisy)
        assert_that(!is.null(obj$sigma.add))
    }

    ## Form the test contrasts
    vlist <- make_all_segment_contrasts(obj)
    vlist <- filter_vlist(vlist, loc)

    ## Obtain p-values
    if(type=="plain"){
        poly.nonfudged = polyhedra(obj)
        pvs = sapply(vlist, function(v){
            pv = poly.pval2(y=obj$y, poly=poly.nonfudged, v=v, sigma=sigma, bits=5000)$pv
        })

    } else if (type=="addnoise") {
        poly.fudged = polyhedra(obj)
        pvs = sapply(vlist, function(v){
            pv = randomize_addnoise(y=obj$y, v=v, sigma=sigma, numIS=10,
                                    sigma.add=obj$sigma.add,
                                    orig.fudged.poly=poly.fudged, bits= 5000,
                                    max.numIS=2000,
                                    min.num.things=min.num.things)$pv})
    } else {
        stop("|type| argument is wrong!")
    }

    ## Add to object
    obj$pvs = pvs
    obj$vlist = vlist
    if(!is.null(mn)){obj$means = sapply(vlist, function(v){ sum(v*mn) })}
    return(obj)
}


##' Appends the inference results to an object created as a result of
##' \code{genlassoinf::dualpathSvd2()}. It is important that, when
##' \code{type=="addnoise"}, that the object is a result of running fused lasso
##' on noise-added response. This is flagged by obj$noisy
##' @param obj object of |path| type.
##' @param loc only test locations in \code{loc}.
##' @param type One of \code{ c("plain", "addnoise")}. If equal to
##'     \code{"addnoise"}, then \code{sigma.add} needs to be provided.
##' @param sigma Noise level (standard deviation) of data.
##' @param sigma.add Additive noise. Defaults to NULL, in which case no additive
##'     noise randomization inference is done.
##' @param mn original mean vector.
##' @export
addpv_fl <- function(obj, loc=NULL, type=c("plain", "addnoise"), sigma,
                       sigma.add=NULL, declutter=FALSE, mn=NULL){

    ## Basic checks
    assert_that(is.null(obj$pvs))
    type = match.arg(type)
    if(type=="addnoise"){
        assert_that(!is.null(obj$noisy))
        assert_that(obj$noisy)
        assert_that(!is.null(obj$sigma.add))
    }


    ## Get randomized p-value
    vlist <- make_all_segment_contrasts(obj)
    vlist <- filter_vlist(vlist, loc)

    ## Obtain p-values
    if(type=="plain"){
        poly.nonfudged = polyhedra_fl(obj)
        pvs = sapply(vlist, function(v){
            pv = poly.pval2(y=obj$y, poly=poly.nonfudged, v=v, sigma=sigma, bits=5000)$pv
        })

    } else if (type=="addnoise") {

        poly.fudged = polyhedra_fl(obj)
        pvs = sapply(vlist, function(v){
            pv = randomize_addnoise(y=obj$y, v=v, sigma=sigma, numIS=10,
                                    sigma.add=sigma.add,
                                    orig.fudged.poly=poly.fudged, bits= 5000,
                                    max.numIS=2000, min.num.things=30)$pv})
    } else {
        stop("|type| argument is wrong!")
    }

    ## Add to object
    obj$pvs = pvs
    obj$vlist = vlist
    if(!is.null(mn)){obj$means = sapply(vlist, function(v){ sum(v*mn) })}
    return(obj)

}

##' Helper to harvest polyhedra from FL object.
polyhedra_fl <- function(obj){
    Gobj = genlassoinf::getGammat.naive(obj=obj, y=obj$y,
                                        condition.step=obj$maxsteps)
    poly = polyhedra(obj=Gobj$G, u=Gobj$u)
    return(poly)
}
