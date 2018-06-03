##' Wild binary segmentation, with fixed threshold
##' @param y numeric vector to contain data
##' @param numSteps numeric of number of stepsemas
##' @param numIntervals number of random intervals to sample
##' @param tol tolerance to handle divide by zero instances
##'
##' @return either an environment address which contains output
##'     (\code{env$intervals} and \code{env$tree}), or a list of wild binary
##'     segmentation output.
##' @export
wildBinSeg_fixedThresh <- function(y, thresh, numIntervals = NULL,
                return.env=FALSE, seed=NULL, verbose=FALSE, intervals = NULL, augment=FALSE){

    ## Basic checks
    if(any(duplicated(y))) stop("y must contain all unique values")
    if(is.null(numIntervals) & is.null(intervals)){
        stop("Provide input for generating intervals, or the intervals themselves!")}
    if(!is.null(intervals)) stopifnot(.is_valid_intervals(intervals))

    ## Generate the random intervals
    if(!is.null(seed)) set.seed(seed)

    if(is.null(intervals)){
        intervals = generate_intervals(length(y), numIntervals)
    }

    ## De-duplicate the intervals
    intervals = .deduplicate_intervals(length(y), intervals)
    starts = intervals$starts
    ends = intervals$ends

    ## Create new environment |env|
    env = new.env()
    env$tree = env$signs = NULL
    env$intervals = intervals ## Carry through the intervals!

    ## Run WBS
    .wildBinSeg_fixedThresh_inner(y, thresh, 1, length(y), 0, 1, verbose, env=env, augment=augment)

    ## Clean the result and return
    .clean_env(env)
    obj = structure(list(tree = env$tree,
                         y = y,
                         thresh = thresh,
                         signs=env$signs,
                         intervals = env$intervals,
                         cp = extract_cp_from_tree(env$tree, "cp"),
                         cp.sign = extract_cp_from_tree(env$tree,
                                                         "sign"),
                         fixedInterval = !is.null(intervals),
                         augment = augment),
                    class = "wbsFt")

    if(return.env){ return(env) } else{ return(obj) }
}

##' Inner function for wbs-ft
.wildBinSeg_fixedThresh_inner <- function(y, thresh, s, e, j, k, verbose=FALSE, env=NULL, augment){

    if(verbose) cat('j,k are', j,k,fill=TRUE)
    if(verbose) cat('s,e are', s,e,fill=TRUE)

    ## If segment is 2-lengthed, terminate
    if(e-s<1){
        return()

    ## Otherwise, calculate CUSUMs
    } else {

        ## Extract the qualifying intervals at this branch
        m = which(.get_which_qualify(s,e,env$intervals))
        if(augment) m = c(m,0)
        if(length(m)==0) return()

        ## Form a matrix of results and save them
        semat = .make_semat(m, s=s, e=e, env$intervals, y, thresh)
        env$tree = addd(env$tree, j, k, semat)

        ## Add the cusum sign information to the |env|ironment
        newsigns = lapply(m, function(single.m){
            .make_signs(single.m, s, e, env$intervals, y)})
        adddd(newsigns, m, env)

        ## If threshold is /not/ exceeded, then terminate
        if(all(semat[,"passthreshold"]==F)){
            return()
        ## If threshold is exceeded, then recurse.
        } else {
            ## Extract changepoin
            passed = which(semat[,"maxhere"] & semat[,"passthreshold"])
            b = semat[passed,"b"]

            ## Recurse
            .wildBinSeg_inner(y, thresh, s,   b,  j+1,  2*k-1, verbose=verbose, env=env, augment=augment)
            .wildBinSeg_inner(y, thresh, b+1, e, j+1, 2*k, verbose=verbose, env=env, augment=augment)
        }
    }
}

#' is_valid for wbs-ft
#' @param obj wbsft object
#' @return TRUE if valid
#' @export
is_valid.wbsFt <- function(obj){
    if(!all(names(obj) %in%
            c("tree",
              "y",
              "thresh",
              "signs",
              "intervals",
              "cp",
              "cp.sign",
              "augment",
              "fixedInterval"))) stop("obj must contain certain elements!")
  TRUE
}


## Print wild binary segmentation
print.wbsFt <- function(obj){
    stopifnot(is_valid.wbsFt(obj))
    print("The changepoints of this object are: ")
    print(obj$cp * obj$cp.sign)
    print("The tree structure is")
    print(obj$tree)
    print("produced from a threshold of")
    print(obj$thresh)
}

## #' is_valid for semat; not using this because it destroys the data frame, makes it into a general structure, which screws code.
## #'
## #' @param semat semat object
## #'
## #' @return TRUE if valid
## #' @export
## is_valid.semat <- function(semat){
##     if(!all(colnames(semat) %in% c("m", "b", "maxcusum", "maxhere", "maxhere",
##                                    "passthreshold"))) stop("semat must be a matrix that contains certain elements!")
##   TRUE
## }

## Collect list of (starts,ends) that qualify at this branch for comparison
.get_which_qualify <- function(s,e, intervals){
    return(sapply(intervals$se, function(se){ return(s<=se[1] && se[2]<=e)}))
}
