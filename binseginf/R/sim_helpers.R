##' Function to generate one-jump mean
##' @export
onejump <- function(lev, n){
    c(rep(0,n/2),rep(lev,n/2))
}

##' Functions to generate two-jump mean
##' @export
twojump <- function(lev,n){c(rep(0,n/3),rep(lev,n/3), rep(0,n/3))}


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

