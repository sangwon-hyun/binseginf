##' Computes the CUSUM (cumulative sum) statistic.
##'
##' Note, we calculate this as the right-to-left difference, by default.
##' @param s starting index.
##' @param b breakpoint index.
##' @param e end index.
##' @param y data.
##' @param n length of data. Defaults to \code{length(y)}.
##' @param right.to.left Whether you want right-to-left difference in the cusum
##'     calculation. Defaults to TRUE.
##' @param contrast.vec If TRUE, then the contrast vector v for cusum=v'y is
##'     returned.
##' @param unsigned if TRUE, then returns the contrast vector that makes the
##'     contrast with |y| positive, or the absolute value of the contrast.
##' @export

cusum <- function(s,b,e,n=NULL, y=NULL, right.to.left = TRUE, contrast.vec = FALSE, unsigned = FALSE, sign.of.vy=NULL){

    ## Form temporary quantities
    nn = e-b+1
    n1 = b+1-s
    n2 = e-b
    ## browser()
    ## if(nn==1) stop("n cannot be 1!") ## Not sure why this was here.
    if(b>=e) stop("b must be strictly smaller than e!")
    if(s>b) stop("b must be larger than or equal to!")
    if(unsigned & is.null(y) & is.null(sign.of.vy)) stop("Can't produce an unsigned linear contrast or contrast vector without a |y| or |sign.of.vy|")
    ## if(!is.null(y)) stopifnot(length(y)==n)
    if(is.null(y)){ if(is.null(n)) stop("Must provide one of y or n.")}
    if(is.null(n)) n = length(y)

    ## Form contrast
    v = rep(0,n)
    v[s:b] = -1/n1
    v[(b+1):e]  = 1/n2
    v = v * sqrt(1/((1/n1)+(1/n2)))
    if(!right.to.left) v = -v

    ## Appropriately sign it
    if(unsigned){
        if(!is.null(y)){
            v = v * sign(sum(v*y))
        } else if (!is.null(sign.of.vy)){
            v = v * sign.of.vy
        } else {
            stop("Provide either |y| or |sign.of.vy| for an unsigned quantity!")
        }
    }

    ## Return the right thing
    if(contrast.vec){
        return(v)
    } else if(!is.null(y)) {
        return(sum(v*y))
    } else {
        stop("Check your options for cusum() again!")
    }
}



##' Function to plot qqlot of p-values. Use extra parameter
##' @param pp numeric vector of p-values.
##' @param main label to plot as main title.
##' @param cols colors if |pp| is a list of numeric vectors.
##' @export
qqunif <- function(pp, main=NULL, plot.it=TRUE, cols=NULL, type=c("p","l"),
                   legend.location="bottomright",
                   ...){

    type = match.arg(type)
    if(type=="l"){
        pch=NULL
        lty=1
    } else {
        pch=16
        lty=NULL
    }

    ## Internal helper
    myplotter <- function(xy,main,...){
            graphics::plot(xy, axes=FALSE, ylim=c(0,1), xlim=c(0,1), xlab="expected",ylab="observed",...)
            graphics::axis(2); graphics::axis(1)
            graphics::abline(0,1)
            if(!is.null(main)) graphics::title(main=main)
    }

    ## If single numeric vector, plot it.
    if(class(pp)!="list"){
        if(is.null(cols)) cols="black"
        xy <- stats::qqplot(x=pp,
                            y=seq(from=0,to=1,length=length(pp)), plot.it=FALSE)
        if(plot.it) myplotter(xy, main, type=type,col=cols,...)
        return(invisible(xy))

    ## Else, if a list of p-values is given, plot all of them
    } else {
        assertthat::assert_that(!is.null(cols))
        allpoints = lapply(pp, function(pvs){qqunif(pvs, plot.it=FALSE)})
        if(plot.it){
            myplotter(allpoints[[1]], main, col=cols[1], pch=16, type=type,...)
            if(length(allpoints)>1){
                for(ii in 2:length(allpoints)){
                    points(allpoints[[ii]], col = cols[ii], pch=16,type=type,...)
                }
            }
        }
        if(length(names(allpoints))>0){
            legend(legend.location,legend=names(pp),col=cols,lty = rep(lty, length(pp)),
                   pch=rep(pch,length(pp)),...)
        }
        return(invisible(allpoints))
    }
}

qqunif_line <- function(pp, main=NULL, cols=NULL,ltys=NULL,lwds=NULL,...){

    ## Internal helper
    myplotter <- function(xy,main,...){
            graphics::plot(xy, axes=FALSE, ylim=c(0,1), xlim=c(0,1),xlab="expected",ylab="observed",...)
            graphics::axis(2); graphics::axis(1)
            graphics::abline(0,1)
            if(!is.null(main)) graphics::title(main=main)
    }

    ## Plot the list of things
    assertthat::assert_that(!is.null(cols))
    if(is.null(ltys))ltys=rep(1,length(pp))
    if(is.null(lwds))lwds=rep(1,length(pp))
    allpoints = lapply(pp, function(pvs){qqunif(pvs, plot.it=FALSE)})
    ## lwds=rep(2,length(pp))
    myplotter(allpoints[[1]], main, col=cols[1], lty=ltys[1], lwd=lwds[1], type='l')
    if(length(allpoints)>1){
        for(ii in 2:length(allpoints)){
            lines(allpoints[[ii]], col = cols[ii], lty=ltys[ii], lwd=lwds[ii])
        }
    }
    if(length(names(allpoints))>0){
        legend("bottomright",legend=names(pp),col=cols,lty=ltys, lwd=lwds)
    }
    return(invisible(allpoints))
}



##' Function to /add/ qq plot points of p-values
##' @param pp numeric vector of p-values.
##' @param main label to plot as main title.
##' @param ... other parameters for \code{qqplot()}.
qqunif_add <- function(pp, main=NULL,...){
    xy <- stats::qqplot(x=pp, y=seq(from=0,to=1,length=length(pp)),plot.it=FALSE)
    graphics::points(xy,...)
}


##' Get all cusums, given start point \code{s} and end point \code{e}
##'
##' @param s Starting index, between \code{1} and \code{n}
##' @param e Ending index, between \code{1} and \code{n}
##' @param y \code{n}-lengthed data vector.
##' @param unsigned \code{TRUE} to return alsolute value.
##' @param simplereturn \code{TRUE} just to return the cusums.
##'
##' @return list of information about the cusum calculations and maximizers in
##'     this interval.
##' @export
getcusums <- function(s,e,y, unsigned=FALSE, simplereturn=FALSE){

    if(s<=0 | e<= 0) stop("must enter valid e,s >=1 ")

    ## Get all cusum
    cusums =  sapply(s:(e-1), function(b){ cusum(s=s,b=b,e=e,y=y, unsigned=unsigned) })
    if(simplereturn)return(cusums)

    names(cusums) = paste("b=",s:(e-1))
    contrasts =  t(sapply(s:(e-1), function(b){cusum(s=s,b=b,e=e,y=y, contrast.vec=TRUE, unsigned = unsigned)}))

    ## Get signs
    signs = sign(cusums)
    abs.cusums = signs*cusums

    return(list(bmax = which.max(abs.cusums)+s-1,
                bmax.cusums = which.max(abs.cusums),
                inds = s:(e-1),
                cusum = max(abs.cusums),
                allcusums = cusums,
                contrasts = contrasts,
                signs=signs))
}

##' Newer, 10 times faster function for getcusum().
getcusums_fast <- function(s, e, cumsums){
    bvec = (s:(e-1))
    n = e-s+1
    cumsums.aug = c(0,cumsums)
    return(-sqrt((e-bvec)/(n*(bvec-s+1)))*(cumsums.aug[bvec+1]-cumsums.aug[s-1+1]) +
        sqrt((bvec-s+1)/(n*(e-bvec)))*(cumsums.aug[e+1]-cumsums.aug[bvec+1]))
}
getcusums2 = getcusums_fast


##' Newer, 10 times faster function for cusum().
cusum_fast <- function(s,e,b,cumsums.aug,cumsums=NULL){
    ## cumsums.aug = c(0,cumsums)
    n = e-s+1
    return(-sqrt((e-b)/(n*(b-s+1)))*(cumsums.aug[b+1]-cumsums.aug[s]) +
        sqrt((b-s+1)/(n*(e-b)))*(cumsums.aug[e+1]-cumsums.aug[b+1]))
}
cusum2 = cusum_fast


##' Does a little more than getcusums2().
get_morethan_cusums2 <- function(s,e,cumsums){
    cusums <- getcusums2(s,e,cumsums)
    max.b.ind = which.max(abs(cusums))
    max.b = max.b.ind + s - 1
    max.z = sign(cusums[max.b.ind])
    ## max.z = sign(cusums[which.max(cusums)])
    return(list(max.b = max.b,
                max.z = max.z,
                max.cusum = cusums[max.b.ind]))
}



## Either return maximizing breakpoint, or the maximum cusum. If \code{m}
## includes zero, then that is handled to correspond to \code{c(s,e)}
## @param m set of indices that correspond to intervals. If equal to 0, coputes things correspond to \code{c(s,e)}.
## @param s start indices of the interval of interest.
## @param e end index of the interval of interest.
## @param interval set of intervals, produced by \code{generate_intervals()}./
.get_max_b <- function(m,s,e, intervals, y, type=c("cusums","max.b", "max.z")){
    type = match.arg(type)
    if(m!=0){
      s <- intervals$starts[[m]]
      e <- intervals$ends[[m]]
    }
    allcusums = getcusums(s=s,e=e,y=y,unsigned=TRUE)$allcusums
    max.b.before.adding.s <- which.max(allcusums)
    if(type=="cusums"){
        return(max(getcusums(s=s,e=e,y=y,unsigned=TRUE)$allcusums))
    } else if (type == "max.b"){
        max.b <- max.b.before.adding.s + (s - 1)
        return(max.b)
    } else if (type == "max.z"){
        allcusums = getcusums(s=s,e=e,y=y,unsigned=FALSE)$allcusums
        max.zs <- sign(allcusums)
        max.z <- max.zs[max.b.before.adding.s]
        return(max.z)
    } else {
        stop("type is not defined")
    }
}


##' Function to crete 1d fused lasso regression matrix.
##' @param m data length
##' @return 1d Fused lasso regression matrix
dual1d_Dmat = function(m){
  D = matrix(0, nrow = m-1, ncol = m)
  for(ii in 1:(m-1)){
    D[ii,ii] = -1
    D[ii,ii+1] = 1
  }
  return(D)
}


##' Creates bootstrap sample of numeric vector \code{vec}, optionally with a
##' seed.
##' @param vec Numeric vector
##' @param seed Seed for random number generation
##' @return resampled \code{vec}.
##' @export
bootstrap_sample <- function(vec, size=length(vec), seed=NULL){
    if(!is.null(seed)) set.seed(seed)
    return(vec[bootstrap_ind(n=length(vec), size=size)])
        ## sample.int(length(vec),replace=TRUE)]
}

##' Creates a set of bootstrap indices.
bootstrap_ind <- function(n, size=n, seed=NULL){
    if(!is.null(seed)) set.seed(seed)
    return(sample.int(n, size=size, replace=TRUE))
}


##' Scale \code{resid} to have empirical std of \code{std}.
##' @param resid numeric vector (meant to be residuals from changepoint model).
##' @param std desired empirical standard deviation
##' @return Properly scaled \code{resid}.
##' @export
scale_resid <- function(resid, std){
    return(resid*(1/stats::sd(resid))*std)
}



##' Only compare things on nodes whose start and end are >1
prune_of_1_length_segments <- function(Tcurr,Scurr,Ecurr){
  Tcurr.copy = Tcurr
  long.enough <- sapply(Tcurr, function(t){
    if(is.null(t)) return(TRUE)
    s = extract(Scurr,t[1],t[2])
    e = extract(Ecurr,t[1],t[2])
    return(e-s>=1)
  })
  Tcurr.copy[!long.enough] = NULL
  return(Tcurr.copy)
}

##' Gets piecewise mean, given changepoint locations.
##' @param y Data.
##' @param cp Changepoints.
##' @return Piecewise constant vector of data whose entries are sample means of
##'     segments.
##' @export
piecewise_mean <- function(y, cp){
  stopifnot(all(1<=cp & cp<=length(y)))
  ## stopifnot(all.equal(sort(cp),cp))
  cp = sort(cp)
  ## segments = lapply(1:(length(cp)+1), function(ii){ v=c(0,cp,length(y));(v[ii]+1):(v[ii+1]) })
  segments = make_segment_inds(cp, length(y))
  segment.means = sapply(segments, function(mysegment){mean(y[mysegment])})
  cleanmn = rep(NA,length(y))
  ## lapply(1:length(segments), function(ii){cleanmn[segments[[ii]]] <<- segment.means[ii]})
  for(ii in 1:length(segments)){
      cleanmn[segments[[ii]]] <- segment.means[ii]
  }
  return(cleanmn)
}



## filters NULL elements out of a list.
.filternull <- function(mylist){
    emptyguys = unlist(lapply(mylist, function(myobj) return(length(myobj)==0)))
    return(mylist[which(!emptyguys)])
}


##' Helper function to take all neighboring-to-each-other clusters,
##' And declutter them by removing all but (rounded up) centroids
##' @param coords.sign optional set of signs.
##' @export
declutter_new <- function(coords, coords.sign=NULL, how.close = 1){

    ## unsorted.coords = coords
    ## unsorted.coords.sign = coords.sign
    coords.order = order(coords)
    coords = coords[coords.order]
    if(is.null(coords.sign)) coords.sign=rep(1,length(coords))
    coords.sign = coords.sign[coords.order]

    ## Error checking
    if(length(coords)<=1){
      if(length(coords)==0) cat('\n',"attempting to declutter", length(coords), "coordinates",'\n')
      return(coords)
    }

    ## Get the clique memberships
    adjacent.diffs = abs(coords[1:(length(coords)-1)] - coords[2:length(coords)])
    cliq.num=1
    cliq.vec=rep(NA,length(coords))
    for(ii in 1:length(adjacent.diffs)){
      if(adjacent.diffs[ii] <= how.close){  ## used to be ==1
        cliq.vec[ii] = cliq.vec[ii+1] = cliq.num
      } else {
        cliq.vec[ii] = cliq.num
        cliq.num = cliq.num + 1
        if(ii==length(adjacent.diffs)){
            cliq.vec[ii+1] = cliq.num
        }
      }
    }

    ## Get the centroids of clique memberships
    ## if(any(table(cliq.vec)>=2)){
    ## browser()

    ## Get cliq centers
    cliq.centers = sapply(1:max(cliq.vec),function(cliq.num){
        floor(mean(which(cliq.vec==cliq.num)))
    })

    ## Get majority vote of signs in each clump
    cliq.signs = sapply(1:max(cliq.vec), function(cliq.num){
        this.clump.signs = coords.sign[which(cliq.vec==cliq.num)]
        sign(sum(this.clump.signs))
    })

    ## Handle when the clump sign is zero
    cliq.signs[which(cliq.signs==0)] = coords.sign[cliq.centers[which(cliq.signs==0)]]

    ## Return processed coords
    processed.coords = coords[cliq.centers] * cliq.signs
    return(processed.coords)
}


##' Return piecewise mean.
##' @param y data vector.
##' @param cp changepoint vector. Assumed to be sorted.
##' @export
get_piecewise_mean <- function(y, cp){
    if(all.equal(sort(cp), cp)!=TRUE) stop ("cp is not sorted!")
    ## y = rnorm(20)
    ## cp = c(5,10)
    aug.cp = c(0,cp,length(y))
    segment.inds = sapply(1:(length(cp)+1),
                          function(ii){ (aug.cp[ii]+1):aug.cp[ii+1]})
    mn = rep(NA,length(y))
    for(ind in segment.inds) mn[ind] <- mean(y[ind])
    return(mn)
}

##' Get changepoint locations from a piecewise constant mean.
##' @param mn data vector.
##' @export
get_cp_from_piecewise_mean <- function(mn, cp){
    diffs = mn[2:length(mn)] - mn[1:(length(mn)-1)]
    return(which(abs(diffs) > 1E-10))
}



##' Takes a named list of n-length contrast vectors, and filters them so that
##' only the contrasts that are desired i.e. contained in \code{locs}.
##'
##' @param locs Only test points in locs
##' @param only.test.nulls If TRUE, only test NULLs.
##' @export
filter_vlist <- function(vlist, locs=NULL, only.test.nulls=FALSE, mn=NULL){

    ## Filter by location
    if(!is.null(locs)){
        retain = which(abs(as.numeric(names(vlist))) %in% locs)
        if(length(retain)==0){
            return(list())
        }
        vlist = vlist[retain]
    }

    ## Only test the null contrasts
    if(only.test.nulls){
        assertthat::assert_that(!is.null(mn))
        means = sapply(vlist, function(v){ sum(v*mn) })
        tol = 1E-10
        which.null = which(abs(means)<tol)
        vlist = vlist[which.null]
    } 

    return(vlist)
}

##' Make n-length vector that induces piecewise mean that is cut at cp. Used for
##' simulation examples.
##' @param y Original data vector.
##' @param cp Changepoint location. Defaults to \code{c()}.
##' @return n-length vector with means
##' @export
make_pw_mean <- function(y, cp=c()){
   n = length(y)
    if(length(cp)==0){
        ord = c()
    } else {
        ord = order(cp)
    }
    cp_aug = c(0, cp[ord], n)
    d = rep(0,n)
    for(ii in 1:(length(cp)+1)){
        ind = (cp_aug[ii]+1):cp_aug[ii+1]
        d[ind] = mean(y[ind])
    }
    return(d)
}


#' Computes the Hausdorff Distance between Two Sets
#' 
#' The two sets, set1 and set2, must be numeric vectors. If one.sided is true,
#' the hausdorff returns max_(a in set1) min_(b in set2) |a-b|.
#' 
#' If either sets are empty, NA is returned.
#'
#' @param set1 a numeric vector
#' @param set2 a numeric vector
#' @param one.sided a logical to determine if a one-sided distance or 
#'   two-sided distance is computed.
#'
#' @return a numeric (length 1)
#' @export
hausdorff <- function(set1, set2, one.sided = F){
  if(!is.numeric(set1) | !is.numeric(set2)) stop(paste("set1 and set2 must be",
    "numerics"))
  
  #handle corner cases
  if(length(set1) == 0 | length(set2) == 0) return(NA)
  if(length(set2) == 1) set2 = c(set2, set2)
  if(length(set1) == 1) set1 = c(set1, set1)

  dist.mat = sapply(set1, function(i){abs(i-set2)})
  dist.vecx = apply(dist.mat, 2, min)
  
  if(!one.sided) dist.vecy = apply(dist.mat, 1, min) else dist.vecy = 0
  
  max(dist.vecx, dist.vecy)
}

#' Compute the Jump Locations
#' 
#' Here, an idx i is a changepoint if i+1 != i.
#'
#' @param obj a numeric vector
#' @param tol the numeric threshold to determine if a location is a changepoint
#' @param ... not used
#'
#' @return a set of numeric integers.
#' @export
jumps.numeric <- function(obj, tol = 1e-10, ...){
  if(!is.numeric(obj)) stop("obj must be numeric")
  if(length(obj) <= 1) return(numeric(0))
  
  dif = abs(diff(obj))
  idx = which(dif > tol)
  
  idx
}



##' Helper for tests; checks of \code{vec} is uniform.
##' @export
expect_uniform <- function(vec){
    expect_equal(ks.test(unlist(vec),punif)$p.value<0.05, FALSE)
}

##' Generate Laplace noie with sigma=1
##' @export
lapl <- function(n,samp=NULL){ rexp(n,rate=sqrt(2)) * sample(c(-1,1),n,replace=TRUE)}




##' Train a binseg changepoint model size using cross validation.
##' @export
cv.bsfs <- function(y, max.numSteps=30, numsplit=2){
    if(length(y)%%2 != 0) y = y[-length(y)] ## Just in case y is odd lengthed
    testerrors = matrix(nrow=max.numSteps, ncol=numsplit)
    testinds = lapply(1:numsplit, function(ii)seq(from=ii, to=length(y), by=numsplit))

    ## Cycle through each partition of the data and calculate test error
    for(jj in 1:numsplit){
        leaveout = testinds[[jj]]
        testData <- y[leaveout]
        trainData <- y[-leaveout]
        obj = bsfs(trainData, max.numSteps)
        for(numSteps in 1:max.numSteps){
            cp = obj$cp[1:numSteps]
            testcp = round(cp / length(trainData) * length(testData))
            testerrors[numSteps,jj] = mean((testData - piecewise_mean(trainData, testcp))^2)
        }
    }
    errors = apply(testerrors, 1, mean)
    return(list(k=which.min(errors), errors=testerrors))
}


##' Helper to make segment indices. 
##' @param cp Changepoints.
##' @param n Data length.
make_segment_inds <- function(cp,n){
    lapply(1:(length(cp)+1),
           function(ii){ v=c(0,cp,n);(v[ii]+1):(v[ii+1]) })
}



##' Linear prediction of vector |ptail| to the next value.
linpredict <- function(ptail,p){
    ind = 1:(length(ptail))
    g = lm(y~x,data=data.frame(y=ptail,x=ind))
    y = predict(g, newdata = data.frame(x=length(ptail)+1))
    return(y)
}


##' Helper function to take all neighboring-to-each-other clusters, and
##' declutter them by removing all but middle value of centroids. The sign is
##' only returned if the signs all agree.
##' @param coords (not necessarily sorted) Coordinates.
##' @param coords.sign Accompanying signs.
##' @param how.close how close you want the cluster members to be.
##' @return List containing |cp| and |cp.sign|, which are lists containing the
##'     de-cluttered changepoints and signs. Additionally, the list contains the
##'     raw, uncluttered versions |cp.raw| and |cp.sign.raw|.
##' @export
declutter <- function(coords, coords.sign, how.close = 1){

    ## Preprocess
    unsorted.coords = coords
    unsorted.coords.sign = coords.sign
    ord = order(unsorted.coords)
    coords = unsorted.coords[ord]
    coords.sign = unsorted.coords.sign[ord]

    ## Define a helper to identify centroid cluster.
    get_center <- function(members){
        mn = mean(members)
        stay = max(members[members<=mn])
        return(stay)
    }

    ## error checking
    if(length(coords)<=1){
      if(length(coords)==0) cat('\n',"attempting to declutter", length(coords), "coordinates",'\n')
      return(coords)
    }

    ## get the clique memberships
    adjacent.diffs = abs(coords[1:(length(coords)-1)] - coords[2:length(coords)])

    cliq.num = 1
    cliq.vec = rep(NA,length(coords))
    for(ii in 1:length(adjacent.diffs)){
      if(adjacent.diffs[ii] <= how.close){  ## used to be ==1
        cliq.vec[ii] = cliq.vec[ii+1] = cliq.num
      } else {
        cliq.num = cliq.num+1
      }
    }
    unique.cliq.nums = unique(cliq.vec[!is.na(cliq.vec)])

    ## Get membership
    members.list = list()
    if(length(unique.cliq.nums)!=0){
        for(ii in 1:length(unique.cliq.nums)){
            cliq.num = unique.cliq.nums[ii]
            members.list[[ii]] = which(cliq.vec == cliq.num)
        }
    }

    ## Add back all others
    all.other = 1:length(unsorted.coords)
    all.other= all.other[!(all.other%in%unlist(members.list))]
    members.list = c(members.list,
                     lapply(all.other, function(a)a))

    ## Determine one or two-sidedness
    signs.list = lapply(members.list, function(members)coords.sign[members])
    one.sided = which(sapply(signs.list, function(signs) all(signs==signs[1])))
    if(length(one.sided)==0){
        two.sided = (1:length(members.list))
    } else {
        two.sided = (1:length(members.list))[-one.sided]
    }
    clustered.signs.list = list()
    for(ii in one.sided){
        clustered.signs.list[[ii]] = signs.list[[ii]][1]
    }
    clustered.signs.list[two.sided] = NA

    ## Likewise, get centers
    clustered.members.list = unlist(lapply(members.list, get_center))
    clustered.coords.list = unlist(lapply(clustered.members.list, function(members){ coords[members]}))

    ## Also prepare uncluttered output
    unclustered.coords.list = lapply(members.list, function(members){coords[members]})
    unclustered.signs.list = signs.list

    output = list(cp = clustered.coords.list,
                  cp.sign = clustered.signs.list,
                  cp.raw = unclustered.coords.list,
                  cp.sign.raw = unclustered.signs.list)
    
    return(output)
}


##' A helper function to print the progress of a simulation.
printprogress <- function(isim, nsim, type="simulation", lapsetime=NULL,
                          lapsetimeunit="seconds", start.time=NULL,
                          fill=FALSE){

    ## If lapse time is present, then use it
    if(fill) cat(fill=TRUE)
    if(is.null(lapsetime) & is.null(start.time)){
            cat("\r", type, " ", isim, "out of", nsim)
    } else {
        if(!is.null(start.time)){
            lapsetime = round(difftime(Sys.time(), start.time,
                                       units = "secs"), 0)
            remainingtime = round(lapsetime * (nsim-isim)/isim,0)
            endtime = Sys.time() + remainingtime
        }
        cat("\r", type, " ", isim, "out of", nsim, "with lapsed time",
            lapsetime, lapsetimeunit, "and remaining time", remainingtime,
            lapsetimeunit, "and will finish at", strftime(endtime))
    }
    if(fill) cat(fill=TRUE)
}
