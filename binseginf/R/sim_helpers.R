##' Function to generate one-jump mean
##' @export
onejump <- function(lev, n){
    c(rep(0,n/2),rep(lev,n/2))
}

##' Functions to generate two-jump mean
##' @export
twojump <- function(lev,n){c(rep(0,n/3),rep(lev,n/3), rep(0,n/3))}


##' Functions to generate two-jump mean
##' @export
twonarrowjump <- function(lev,n){stopifnot(n%%20==0);c(rep(0,n/2-n/20),rep(lev,n/10), rep(0,n/2-n/20))}

##' Function to generate four-jump mean
##' @export
fourjump <- function(lev,n){
    c(rep(0,n/5), rep(lev,n/5), rep(0,n/5), rep(-2*lev, n/5), rep(0,n/5) )
}

##' Function to generate four-jump mean /with same size jumps/.
##' @export
fourjump_samesize <- function(lev,n){
    c(rep(0,n/5), rep(lev,n/5), rep(0,n/5), rep(-1*lev, n/5), rep(0,n/5) )
}


##' Function to generate four-jump mean with narrow spikes
##' @export
fourjump_spiky <- function(n, lev){
    mn = rep(0, n)
    mn[n/3+seq(from=-n/100, to = n/100)] = lev
    mn[2*n/3+seq(from=-n/100, to = n/100)] = -2*lev
    return(mn)
}

##' Function to generate four-jump mean with one spike and one plateau
##' @export
fourjump_hybrid <- function(n,lev){
    mn = rep(0, n)
    mn[n/3+seq(from=-n/100, to = n/100)] = lev
    mn[2*n/3+seq(from=-n/10, to = n/10)] = -2*lev
    return(mn)
}



##' Simple function to do all binary segmentation segment test after sample
##' splitting.
##' @param y data vector.
samplesplit <- function(y, mn, numSteps, only.test.nulls){
    assert_that(length(y)%%2==0)
    ind = seq(from=1, to=length(y),by=2)
    y1 = y[ind]
    mn1 =  mn[ind] ## This is actually tricky business, because what might be
                  ## exactly zero in ind1 might not actually be zero in ind2
    y2 = y[-ind]
    obj = bsfs(y1, numSteps)
    vlist = make_all_segment_contrasts(obj)
    vlist = filter_vlist(vlist, mn=mn1, only.test.nulls=only.test.nulls)
    ## print(mn1)
    ## print(sapply(vlist, function(v){sum(v*mn1)}))
    pvs = sapply(vlist, function(v){  pv = ztest(y2, v) })
    return(pvs)
}
