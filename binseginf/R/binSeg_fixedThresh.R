##' Main function for binary segmentation for fixed threshold. This is actually
##' a wrapper for binary segmentation with fixed threshold. It creates an
##' environment and creates the variables there, then runs
##' binseg.by.thresh.inner() all in this environment, and returns the relevant
##' guy
##' @param s Starting index, in the vector-valued data. Must be an integer
##'     larger than or equal to 1, and strictly smaller than \code{e}.
##' @param e Ending index, in the vector-valued data. Must be an integer smaller
##'     than or equal to \code{n}, and strictly larger than \code{s}.
##' @param thresh Threshold for \deqn{\tilde{X}}. This serves as a stopping rule
##'     for the recursion.
##' @param y The original data.
##' @param verbose Set to true if you'd like to see algorithm progression.
##' @param return.env Set to true if you'd like to get the environment
##'     containing the algorithm output, instead of a list.
##' @import pryr
binSeg_fixedThresh = function(y, thresh, s=1, e=length(y), verbose=FALSE, return.env=FALSE, numIntervals=NULL){

    n = length(y)

    ## Create new environment |env|
    env = new.env()
    env$infotable = data.frame(j=1,k=1,s=1,b=1,e=1,pass=TRUE,dir=+1,len=0)[-1,]
    env$y = y

    ## Run binary segmentation on |env|
    binseg.by.thresh.inner(y, thresh, s, e, 1, 1, verbose, env=env)
    rownames(env$infotable) = c()

    ## Gather output from |env| and return it.
    obj = structure(list(infotable = env$infotable,
                         y = y,
                         thresh = thresh,
                         cp = env$infotable[which(env$infotable[,"pass"]),"b"],
                         cp.sign = env$infotable[which(env$infotable[,"pass"]),"dir"]),
                    class = "bsFt")

    if(return.env){ return(env) } else{ return(obj) }

}

##' Checks of object is of class "bsFt", produced by binSeg_fixedThresh().
##' @param obj bsFt object
is_valid.bsFt <- function(obj){
    if(!(all(names(obj) %in% c("infotable",
                               "y",
                               "thresh",
                               "cp",
                               "cp.sign")))) stop("obj must contain certain elements!")
    TRUE
}



##' Inner function for binary segmentation with fixed threshold. The wrapper
##' \code{binseg.by.thresh()} is intended to be used by user. Note, when
##' \code{thresh} is set to zero, then this can be used to collect the
##' unbalanced haar wavelet basis.
##' @param s Starting index, in the vector-valued data. Must be an integer
##'     larger than or equal to 1, and strictly smaller than \code{e}.
##' @param e Ending index, in the vector-valued data. Must be an integer smaller
##'     than or equal to \code{n}, and strictly larger than \code{s}.
##' @param j The depth of the recursion on hand.
##' @param k The indexing of the node location, from left to right, in the
##'     \emph{complete} binary tree.
##' @param thresh Threshold for \deqn{\tilde{X}}. This serves as a stopping rule
##'     for the recursion.
##' @param y The original data.
##' @param n The length of the data \code{y}.

binseg.by.thresh.inner <- function(y, thresh, s=1, e=length(y), j=1, k=1, verbose=F, env=NULL){

    n = length(y)

    ## If segment is 2-lengthed, terminate
    if(e-s<1){
        newrow = data.frame(j=j,k=k,s=s,b=NA,e=e,pass=NA,dir=NA,len=e-s)
        env$infotable = rbind(env$infotable, newrow)

    ## Otherwise, calculate CUSUMs
    } else {
        all.bs = (s:(e-1))
        all.cusums = sapply(all.bs, function(b) cusum(s=s, b=b, e=e, y=y))
        names(all.cusums) = all.bs
        all.cusums = c(rep(NA,s-1), all.cusums, rep(NA,n-s))
        sn = sign(all.cusums)

        ## Obtain breakpoint and its sign
        b = which.max(abs(all.cusums))
        z = sn[b]
        if(verbose) cat("the (potential) changepoint", b, "is being considered, between s and e:",s,e, fill=T)

        ## Check threshold exceedance, then store
        if(abs(all.cusums[b]) < thresh){

            newrow = data.frame(j=j,k=k,s=s,b=b,e=e,pass=FALSE,dir=z,len=e-s)
            env$infotable = rbind(env$infotable, newrow)
            return(env)
        } else {
            if(verbose) cat("the biggest cusum was", all.cusums[b],
                            "which passed the threshold:", thresh,
                            "from between s and e:",s,e, fill=T)

            newrow = data.frame(j=j,k=k,s=s,b=b,e=e,pass=TRUE,dir=z,len=e-s)
            env$infotable = rbind(env$infotable, newrow)
        }

        ## Recurse
        binseg.by.thresh.inner(y, thresh, s, b, j+1, 2*k-1, verbose, env=env)
        binseg.by.thresh.inner(y, thresh, b+1, e, j+1, 2*k, verbose, env=env)
    }
}



##' Print function for bsFt class
print.bsFt <- function(obj){
    cat("Changepoint set is ", obj$cp*obj$cp.sign, "produced from threshold =", obj$thresh, fill=TRUE)
}
