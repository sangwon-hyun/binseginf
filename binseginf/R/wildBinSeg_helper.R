##TODO: Most of these are obsolete. Get rid of the ones that are not used.

##' Adds a matrix containing information about wbs selection even, as a last
##' element to an existing list of such matrices. the name of this new element
##' is of the form "j,k". this is in lieu of |cplist| which was used to store
##' information for sbs-fixed-threshold events. also note: more than one
##' maximizer per branch is /fine/!
##' @param mylist a list of matrices
addd = function(mylist,j,k,semat){

    ## Basic checks: list is empty? elements are all correctly formated matrices?
    ## if(length(mylist)==0) stop("you are trying to add to an empty list!")
    ## mylist=null
    ## if(length(mylist)==0) mylist= list(1)
    ## if(!all(sapply(1:length(mylist), function(ii) colnames((mylist)[[ii]])[1]=="m"))) stop("mylist is empty!")

    ## Augment the list with a new element and name it appropriately
    mylist = c(mylist, list(semat))
    names(mylist)[length(mylist)] = paste(j, k, sep = ",")

    ## return the list
    return(mylist)
}

##' Make an addition to the sign list. Specifically, adds elements
##' @param newsigns the signs to add.
##' @param M the indices in env$intervals
##' @param env environment where signs are embedded in
adddd = function(newsigns,M,env){

    ## Basic checks
    stopifnot(class(newsigns) == "list")

    ## Append the new signs
    Map(function(newsign, m){
        env$signs = c(env$signs,list(newsign))
        names(env$signs)[length(env$signs)] = m
        return()
    },newsigns, M)


    ## No return
    return()
}


##' Helper function for a /single/ start and end.
.maximize <- function(s, e, y, getb=TRUE){

    ## collect all cusums and signs
    all.bs = (s:(e-1))
    all.cusums = sapply(all.bs, function(b){cusum(s=s, b=b, e=e, y=y)})
    names(all.cusums) = all.bs
    all.cusums = c(rep(NA,s-1), all.cusums, rep(NA,length(y)-s))
    sn = sign(all.cusums)

    ## obtain cusum-maximizing breakpoint and its sign
    b = which.max(abs(all.cusums))
    z = sn[b]

    ## return maximizer and sign
    if(getb){ return(b) } else {  return(all.cusums[b]) }

}



##' Matrix of selection information at this branch
##' @param m indices of intervals that are to be considered (i.e. are contained
##'     in the relevant interval (s:e) at runtime of wbs()
##' @param s start location of the current binseg call.
##' @param e end location of the current binseg call.
##' @param intervals set of random intervals, drawn between 1 and 60
.make_semat = function(m, s, e, intervals, y, thresh){

    ## Make bare matrix
    mymat = matrix(NA, ncol=7, nrow = length(m), dimnames=NULL)
    mymat = data.frame(mymat)
    colnames(mymat) = c("m", "s", "b", "e", "maxcusum", "maxhere",
                        "passthreshold")

    ## Handle the case of m containing zeros.
    which.zero = which(m==0)
    if(length(which.zero)>0){
        m.without.zeros = m[-which.zero]
        se.for.each.m = c(intervals$se[m.without.zeros], list(c(s,e)))
    } else {
        se.for.each.m = intervals$se[m]
    }

    ## Fill in information about selection
    mymat[,"m"] = m
    mymat[,"s"] = sapply(se.for.each.m, function(vec)vec[1])
    mymat[,"b"] = sapply(se.for.each.m,
                       function(se){
                           .maximize(se[1], se[2], y, TRUE) })
    mymat[,"e"] = sapply(se.for.each.m, function(vec)vec[2])
    mymat[,"maxcusum"] = sapply(se.for.each.m,
                       function(se){ maximize(se[1], se[2], y, FALSE) })

    ## Mark which guy is the maximizing interval
    mymat[,"maxhere"] = rep(FALSE,length(m))
    mymat[which.max(abs(mymat[,"maxcusum"])),"maxhere"] = TRUE

    ## Check if that guy's max cusum passed the threshold
    mymat[,"passthreshold"] = rep(FALSE,length(m))
    max.ind = which.max(abs(mymat[,"maxcusum"]))
    passed.thresh = (abs(mymat[max.ind,"maxcusum"]) > thresh)
    mymat[which.max(abs(mymat[,"maxcusum"])),"passthreshold"] = passed.thresh


    return(mymat)
}


##' Matrix of selection information at this branch
##' @param m index of intervals that are to be considered (i.e. are contained
##'     in the relevant interval (s:e) at runtime of wbs()
##' @param s start location of the current binseg call.
##' @param e end location of the current binseg call.
##' @param intervals set of random intervals, drawn between 1 and 60
.make_signs <- function(m, s, e, intervals, y, thresh){

    ## Basic checks
    stopifnot(length(m)==1)

    ## Calculate e/hings
    if(m==0){
      se = c(s,e)
    } else {
      se = intervals$se[[m]]
    }
    s=se[1]; e=se[2];
    all.bs = (s:(e-1))
    all.cusums = sapply(all.bs, function(b){cusum(s=s, b=b, e=e, y=y)})
    names(all.cusums) = all.bs
    sn = sign(all.cusums)


    ## Make bare matrix
    mymat = matrix(NA, ncol=6, nrow = e-s, dimnames=NULL)
    mymat = data.frame(mymat)
    colnames(mymat) = c("s", "b", "e", "sn", "maxhere", "cusums")

    ## Fill in information about selection
    mymat[,"s"] = rep(s, e-s)
    mymat[,"b"] = all.bs
    mymat[,"e"] = rep(e, e-s)
    mymat[,"sn"] = sn
    mymat[, "maxhere"] = FALSE
    mymat[which.max(abs(all.cusums)), "maxhere"] = TRUE
    mymat[,"cusums"] = all.cusums

    return(mymat)
}



##' Function to generate random intervals for wild binary segmentation.
##' @param n length of data.
##' @param numIntervals Number of intervals to draw.
##' @param seed seed number for random interval generation; defaults to NULL.
##' @param start.end.list Manual list of starts and ends. Literally an R list
##' with two equal length vectors, each named |start| and |end|
##' @return List containing starts and ends and intervals etc.
##' @export
generate_intervals <- function(n, numIntervals, seed=NULL, start.end.list = NULL){

    ## Basic checks
    stopifnot(n>1)

    if(!is.null(seed)) set.seed(seed)
    if(!is.null(start.end.list)){
        stopifnot(names(start.end.list) %in% c("start","end"))
        stopifnot(length(start.end.list[["start"]]) == length(start.end.list[["end"]]))
    }

    ## Use inputted start.end.list, or draw intervals
    if(!is.null(start.end.list)){
       starts = start.end.list[[1]]
       ends = start.end.list[[2]]
       reverses = (starts >= ends)
       duplicates = (starts == ends)
       numIntervals = length(starts)
    } else {
        done.drawing = FALSE
        while(!done.drawing){
            starts = sample(1:n, numIntervals*3, replace=TRUE)
            ends = sample(1:n, numIntervals*3, replace=TRUE)
            reverses = (starts >= ends)
            duplicates = (starts == ends)
            done.drawing = (3*numIntervals-sum(duplicates)> numIntervals)
            if(numIntervals==0) done.drawing = TRUE
        }
    }

    ## Function to make interval given start and end indices, not necessarily
    ## start<end, with an option to reverse them, or rarely, return NULL when
    ## start == end
    makeInterval = function(start, end, reverse, duplicate){
        if(duplicate) return(NULL)
        if(reverse){
            return(end:start)
        } else {
            return(start:end)
        }
    }

    ## Take intervals (identical intervals are eliminated! since they don't play
    ## a role anywhere further along the way, and the max-CUSUM comparisons made
    ## using the de-duplicated set of drawn intervals is still fair/same)
    intervals = Map(makeInterval, starts, ends, reverses, duplicates)
    intervals = intervals[!duplicates]
    intervals = intervals[1:numIntervals]

    ## Startpoints
    starts = sapply(intervals,function(se)se[1])
    ends = sapply(intervals,function(se)se[length(se)])

    ## return
    return(structure(list(starts = starts,
                          ends = ends,
                          intervals = intervals,
                          se = Map(c,starts,ends)),class="intervals"))
}



##' Helper function that takes in a tree (list of |semat|'s), and extracts the
##' changepoints
##' @param tree A list of |semat|'s. From the environment |env| you created with
##'     \code{wbs()}, simply use \code{env$tree}.
##' @param returntype One of \code{c("cp","sign")}, for whether to return the
##'     changepoint or the sign of teh changepoints.
extract_cp_from_tree = function(tree, returntype = c("cp","sign")){

    ## Extract changepoints and sign from the tree
    if(returntype == "cp"){
        all.cps = sapply(tree,
                         function(semat){
                             passed = which(semat[,"maxhere"] & semat[,"passthreshold"])
                             return(semat[passed,"b"])})
        return(unlist(all.cps))
    } else if (returntype == "sign"){
        all.signs = sapply(tree,
                           function(semat){
                               passed = which(semat[,"maxhere"] & semat[,"passthreshold"])
                               return(sign(semat[passed,"maxcusum"]))})
        return(unlist(all.signs))
    } else {
        stop("|returntype| argument should be one of \"cp\" or \"sign\"")
    }
}



##' Wrapper for \code{cusum()}, to produce /signed/ cusum (unsigned means that
##' the < contrast , y> is manually made to be positive, and signed means it retains its original
##' sign) contrast given \code{s,b,e} as start, break and end.
##' @param s starting index.
##' @param b breakpoint index.
##' @param e end index.
##' @param y data.
signed_contrast <- function(s,b,e,n=NULL,y){
            cusum(s,b,e,n,y, contrast.vec=TRUE, unsigned=FALSE)
}

##' Wrapper for \code{cusum()}, to produce /unsigned/ cusum contrast given
##' \code{s,b,e} as start, break and end. (Unsigned means that the < contrast ,
##' y> is manually made to be positive, and signed means it retains its original
##' sign).
##' @param s starting index.
##' @param b breakpoint index.
##' @param e end index.
##' @param y data.
unsigned_contrast <- function(s,b,e,n=NULL,y){
            cusum(s,b,e,n,y, contrast.vec=TRUE, unsigned=TRUE)
}


##' Helper function to /manually/ make contrasts.
##' @param test.bps Breakpoint location to test.
##' @param adj.bps Directly adjacent breakpoint locations.
##' @param sn Sign (direction) of proposed breakpoint at |test.bps|; use +1 or -1.
##' @param n length of data.
##' @examples
##' make_contrast(20,c(1,40),+1,60)
##' @export
make_contrast <- function(test.bp, adj.bps, sn, n){

    ## Basic checks
    stopifnot(all(c(test.bp, adj.bps) %in% 0:n))
    stopifnot(min(adj.bps)<=test.bp)
    stopifnot(max(adj.bps)>=test.bp)
    stopifnot(length(sn)==1)
    stopifnot(sn %in% c(-1,1))

    ## Make contrast
    d = rep(0,n)
    d[(min(adj.bps)):(test.bp)] = -1/(test.bp-min(adj.bps)+1)
    d[(test.bp+1):(max(adj.bps))] = +1/(max(adj.bps)-test.bp)
    return(sn*d)
}
make_segment_contrast = make_contrast


## Checking if intervals is correct.
.is_valid_intervals <- function(intervals){
   return(all(names(intervals) %in% c("starts","ends","intervals","se")))
}


## Deduplicating any intervals.
.deduplicate_intervals <- function(n, intervals){
    ## Basic checks
    stopifnot(.is_valid_intervals(intervals))
    if(all(sapply(intervals$se, is.null))) return(generate_intervals(n=n,numIntervals=0))

    ## Get unique guys, form new intervals
    unique.se = unique(intervals$se)
    unique.start.end.list = list(sapply(unique.se, function(se)se[1]),
                                 sapply(unique.se, function(se)se[2]))
    return(generate_intervals(n=n, start.end.list = unique.start.end.list))

}
