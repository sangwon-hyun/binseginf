##' Constructor for intervals object.
##' @param nrow Number of rows in the empty matrix
##' @param existing A 2-row numeric matrix containing start and end points.
##'     Column names should be names "s" and "e" respectively.
##' @param distance minimum distance between s and e, i.e. e is at least
##'     s+distance.
##' @param maxlength maximum distance between s and e, i.e. (e-s <= maxlength)
##' @return creates an all-NA matrix of dimension nrow x 3. The first two
##'     columns must be the numeric (no check yet), but the last column can be
##'     of any type you want. Initializes to numeric.
##' @export
intervals <- function(numIntervals, n, comprehensive=FALSE, existing=NULL, distance=0, maxlength=n-1) {

    ## Basic checks
    if(!is.null(existing)){
        assert_that(all(colnames(existing) %in% c("s", "e")),
                   msg = "|existing| must be a 2-column matrix with columns names s and e.")
    }

    ## If |comprehensive| == TRUE, draw /all/ possible intervals.
    if(comprehensive){
        all.se = t(combn(n,2)) ## This is computationally intensive and quite unnecessary.
        x.all = all.se[,1]
        y.all = all.se[,2]
        if(!is.null(existing))stop("Can't do both of: have comprehensive intervals, and exclude |existing| intervals!")
    } else {
        ## numIntervals = 10000
        enough.intervals = FALSE
        x.all = c()
        y.all = c()
        num.valid.intervals = 0
        fac = 5 ## This is just a way to sample a lot initially
        while(!enough.intervals){
            x = sample(n, size=numIntervals*fac, replace=TRUE)
            y = sample(n, size=numIntervals*fac, replace=TRUE)
            too.close = which(unlist(Map(function(a,b){b-a<=distance}, x, y)))
            too.long = which(unlist(Map(function(a,b){b-a>maxlength}, x, y)))
            eliminate = c(too.close, too.long)
            x.all = c(x.all, x[-eliminate])
            y.all = c(y.all, y[-eliminate])
            ## If |existing| matrix is supplied, then exclude these from consideration.
            if(!is.null(existing)){
                to.exclude = unlist(apply(existing,1, function(myrow){
                    which((x.all == myrow["s"]) & (y.all == myrow["e"]))
                }))
                x.all = x.all[-to.exclude]
                y.all = y.all[-to.exclude]
            }
            enough.intervals = (length(x.all) >= numIntervals)
        }
    }

    ## Visualizing the frequency; should be even barplots!
    ## xy = cbind(x.all, y.all)
    ## plot((xy), xlim=c(0,n), ylim=c(0,n))
    ## xylist = lapply(1:nrow(xy), function(irow)xy[irow,])
    ## allpairs = unique(xylist)
    ## freqs = rep(NA,length(allpairs))
    ## for(ipair in 1:length(allpairs)){
    ##     mypair = allpairs[[ipair]]
    ##     freqs[ipair] = sum(sapply(xylist, function(myxy)all(myxy==mypair)))
    ## }
    ## barplot(freqs)

    structure(list(starts=x.all,
                   ends=y.all,
                   cusummat=NA,
                   numIntervals=numIntervals), class="intervals")
}

##' Method to check if object is a valid object of the |interval| class.
is_valid.intervals <- function(obj){
    return(all(names(intervals) %in% c("starts", "ends", "cusummat", "numIntervals")))
}

addcusum <- function(intervals,...){ UseMethod("addcusum")}

##' Calculates all cusums
addcusum.intervals <- function(intervals, y){

    cumsums = cumsum(y)
    cusummat = matrix(NA, nrow=length(y)*intervals$numIntervals, ncol=4)
    colnames(cusummat) = c("s","b", "e", "cusum")
    tracker.i = 0
    for(ii in 1:intervals$numIntervals){
        s = intervals$starts[ii]
        e = intervals$ends[ii]
        cusums = getcusums2(s=s, e=e, cumsums)
        ncusum = length(cusums)
        irow = tracker.i + (1:ncusum)
        cusummat[irow, "cusum"] = cusums
        cusummat[irow, "s"] = rep(s, ncusum)
        cusummat[irow, "e"] = rep(e, ncusum)
        cusummat[irow, "b"] = c(s:(e-1))
        tracker.i = tracker.i + ncusum
    }
    intervals$cusummat = cusummat[1:tracker.i,]
    return(intervals)
}

restrict <- function(intervals,...){ UseMethod("restrict")}

##' Sees which guys in object correspond to cusums calculated between |s| and
##' |e|. In other words, returns which guys are /qualified/.
restrict.intervals <- function(intervals, s=1, e=n){
    match.i = (intervals$cusummat[,c("s")] >= s) & (intervals$cusummat[,c("e")] <= e)
    if(length(match.i) == 0) stop("No start and end points correspond to")
    return(which(match.i))
}


maximize <- function(intervals,...){ UseMethod("maximize")}
##' Get the cusum maximizer, among qualifying guys.
maximize.intervals <- function(intervals, qual.inds){
    max.i = qual.inds[which.max(abs(intervals$cusummat[qual.inds,"cusum"]))]
    max.sign = sign(intervals$cusummat[max.i,"cusum"])
    return(list(max.i=max.i,
                max.sign=max.sign,
                max.s=intervals$cusummat[max.i,"s"],
                max.b=intervals$cusummat[max.i,"b"],
               max.e=intervals$cusummat[max.i,"e"]))}


form_rows <- function(intervals,...){ UseMethod("form_rows")}
##' Forms polyhedron's rows, given a maximizing b and all the qualifying indices
##' @param y Defaults to NULL, but you can input y to assert that your contrasts
##'     should be correct.
form_rows.intervals <- function(intervals, max.s, max.b, max.e, max.sign, qual.inds, n, y=NULL){

    winning.contrast <- cusum(s=max.s, b=max.b, e=max.e, unsigned=TRUE, n=n, sign.of.vy=max.sign,
                              contrast.vec=TRUE)
    submat = intervals$cusummat[qual.inds,,drop=FALSE]

    rows1 <- t(apply(submat, 1, function(myrow){
        winning.contrast - cusum(s=myrow["s"], b=myrow["b"], e=myrow["e"],
              contrast.vec=TRUE, unsigned=FALSE, n=n)
    }))

    rows2 <- t(apply(submat, 1, function(myrow){
        winning.contrast + cusum(s=myrow["s"], b=myrow["b"], e=myrow["e"],
                                 contrast.vec=TRUE, unsigned=FALSE, n=n)
    }))

    rows3 = rbind(winning.contrast)

    all.rows = rbind(rows1, rows2, rows3)

    ## Optionally, assert that all halfspaces should actually contain y
    if(!is.null(y)) assert_that(all(all.rows%*%y>=0))
    return(all.rows)
}

## Forms information (Gy, Gv, and u) based on each winning event.
form_info <- function(intervals,...){ UseMethod("form_info")}
form_info.intervals <- function(intervals, max.s, max.b, max.e, max.sign,
                                qual.inds, cumsum.y, cumsum.v, Gy.ic, Gv.ic, u.ic){
    if(is.null(cumsum.y)|is.null(cumsum.v)) stop("cumsum.v and cumsum.y need to be provided!")

    ## Form winning Gy and Gv
    winning.Gy = cusum_fast(s=max.s, b=max.b, e=max.e, cumsums.aug=c(0,cumsum.y))
    winning.Gv = sign(winning.Gy) * cusum_fast(s=max.s, b=max.b, e=max.e,
                                               cumsums.aug=c(0,cumsum.v))
    winning.Gy = abs(winning.Gy)

    submat = submat.old = intervals$cusummat[qual.inds,1:3, drop=FALSE]
    submat[,2:3] = submat[,2:3] + 1 ## for adding 1 each to b and e

    ## Form G %*% y entries
    K = nrow(submat)
    CMy = c(0,cumsum.y)[submat]
    Gy1 = -sqrt( (submat[,3] - submat[,2])/((submat[,3] - submat[,1])*(submat[,2] - submat[,1]) )) *
        (CMy[(K+1):(2*K)] - CMy[(1:K)]) +
        sqrt( (submat[,2] - submat[,1] )/( (submat[,3]-submat[,1])* (submat[,3] - submat[,2] ) )) *
        (CMy[((2*K+1) :(3*K)) ] - CMy[ ( (K+1  ):(2*K) ) ]  )
    Gy = c(winning.Gy - Gy1, winning.Gy + Gy1, winning.Gy)

    ## Form G %*% v entries in the /exact same way/
    ## Gv1 <- apply(submat, 1, gv)
    CMv = c(0,cumsum.v)[submat]
    Gv1 = -sqrt( (submat[,3] - submat[,2])/((submat[,3] - submat[,1])*(submat[,2] - submat[,1]))) *
        (CMv[(K+1):(2*K)] - CMv[(1:K)]) +
        sqrt( (submat[,2] - submat[,1] )/( (submat[,3]-submat[,1])* (submat[,3] - submat[,2] ) )) *
        (CMv[((2*K+1) :(3*K)) ] - CMv[ ( (K+1  ):(2*K) ) ]  )
    Gv = c(winning.Gv - Gv1, winning.Gv + Gv1, winning.Gv)

    ## Also add things about IC if applicable.
    u = c(rep(0, length(Gv)), u.ic)
    Gy = c(Gy, Gy.ic)
    Gv = c(Gv, Gv.ic)
    stopifnot(length(Gy)==length(Gv) & length(Gy) == length(u))

    return(list(Gy=Gy, Gv=Gv, u=u))
}





##' Add a single interval to existing set of intervals
##' @param old.intervals an object of class |intervals|
##' @param new.s new start point
##' @param new.e new end point
##' @export
add.intervals <- function(intervals, new.s, new.e){
    ## Basic checks
    stopifnot(new.s < new.e)

    ## Append new interval information
    new.starts = c(intervals$starts, new.s)
    new.ends = c(intervals$ends, new.e)

    ## Return as |intervals| class object
    structure(list(starts=new.starts,
                   ends=new.ends),
              class="intervals")
}


add2 <- function(intervals,...){ UseMethod("add2")}

##' Add winning intervals from the original object, to the |intervals| object.
add2.intervals <- function(intervals, winning.wbs.obj, stop.time=NULL){

    ## Make ic-stopped results table
    if(is.null(stop.time)) stop.time = nrow(winning.wbs.obj$results)
    winning.results = winning.wbs.obj$results[1:stop.time,,drop=FALSE]

    ## Form new intervals
    intervals$starts = c(intervals$starts, winning.results[1:stop.time,"max.s"])
    intervals$ends = c(intervals$ends, winning.results[1:stop.time,"max.e"])
    intervals$numIntervals = intervals$numIntervals + stop.time
    return(intervals)
}


##' After making intervals, you can attempt to plot them.
##' @export
plot.intervals <- function(obj){

    ## Basic checks
    stopifnot(is_valid.intervals(obj))

    graphics::plot(NA,
                ylim = c(0, obj$numIntervals),
                xlim = c(0,max(obj$ends)),
                xlab = "intervals",
                ylab = "")
    for(ii in 1:numIntervals){
        graphics::lines(x=c(obj$starts[ii], obj$ends[ii]), y = c(ii,ii))
    }
}
