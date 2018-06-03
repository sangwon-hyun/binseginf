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




##' From polyhedron and data vector and contrast, gets the probability that vtY
##' is in the polyhedron, conditional on PvperpY. i.e. the probability that Y
##' stays in the polyhedron, fixing n-1 dimensions of it.
##' @param y data
##' @param poly polyhedra object produced form \code{polyhedra(wbs_object)}
##' @param sigma data noise (standard deviation)
##' @param nullcontrast the null value of \eqn{v^T\mu}, for \eqn{\mu = E(y)}.
##' @param v contrast vector
##'
##' @return list of two vectors: denominators and numerators, each named
##'     \code{denom} and \code{numer}.
partition_TG <- function(y, poly, v, sigma, nullcontrast=0, bits=50, reduce,correct.ends=FALSE, shift=NULL, ic.poly=NULL, warn=TRUE){

    ## Basic checks
    stopifnot(length(v)==length(y))

    vy = sum(v*y)
    vv = sum(v^2)
    sd = sigma*sqrt(vv)

    ## Shift polyhedron by a constant \R^n shift if needed
    if(!is.null(shift)){
        stopifnot(length(shift)==length(y))
        poly$u = poly$u - poly$gamma%*%shift
    }

    ## Add stopping component to the polyhedron at this point, if needed
    if(!is.null(ic.poly)){
        poly$gamma = rbind(poly$gamma, ic.poly$gamma)
        poly$u = c(poly$u, ic.poly$u)
    }

    ## Just in case |poly| doesn't contain |vup| and |vlo|, we manually form it.
    ## This is because in order to partition the TG statistic, we need to form
    ## these anyway.
    pvobj <- poly.pval2(y, poly, v, sigma, correct.ends=correct.ends)
    vup = pvobj$vup
    vlo = pvobj$vlo
    vy = max(min(vy, vup),vlo)

    ## Make it so that vlo<vup is ensured

    ## Calculate a,b,z for TG = (F(b)-F(z))/(F(b)-F(a))
    z = Rmpfr::mpfr(vy/sd, precBits=bits)
    a = Rmpfr::mpfr(vlo/sd, precBits=bits)
    b = Rmpfr::mpfr(vup/sd, precBits=bits)
    if(!(a<=z &  z<=b) & warn){
        warning("F(vlo)<vy<F(vup) was violated, in partition_TG()!")
    }

    ## Separately store and return num&denom of TG
    numer = as.numeric(Rmpfr::pnorm(b)-Rmpfr::pnorm(z))
    denom = as.numeric(Rmpfr::pnorm(b)-Rmpfr::pnorm(a))

    ## Form p-value as well.
    pv = as.numeric(numer/denom)
    ## if(!(0 <= pv & pv <= 1)) print("pv was not between 0 and 1, in partition_TG()!")

    return(list(denom=denom, numer=numer, pv=pv, vlo=vlo, vy=vy, vup=vup))
}




##' Function to plot qqlot of p-values. Use extra parameter
##' @param pp numeric vector of p-values.
##' @param main label to plot as main title.
##' @param cols colors if |pp| is a list of numeric vectors.
##' @export
qqunif <- function(pp, main=NULL, plot.it=TRUE, cols=NULL, type=c("p","l"),...){

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
            graphics::plot(xy, axes=FALSE, ylim=c(0,1), xlim=c(0,1),xlab="expected",ylab="observed",...)
            graphics::axis(2); graphics::axis(1)
            graphics::abline(0,1)
            if(!is.null(main)) graphics::title(main=main)
    }

    ## If single numeric vector, plot it.
    if(class(pp)!="list"){
            
        xy <- stats::qqplot(x=pp,
                     y=seq(from=0,to=1,length=length(pp)), plot.it=FALSE, pch=16, type=type)
        if(plot.it) myplotter(xy, main)
        return(invisible(xy))

    ## Else, if a list of p-values is given, plot all of them
    } else {
        assert_that(!is.null(cols))
        allpoints = lapply(pp, function(pvs){qqunif(pvs, plot.it=FALSE)})
        if(plot.it){
            myplotter(allpoints[[1]], main, col=cols[1],pch=16, type=type)
            if(length(allpoints)>1){
                for(ii in 2:length(allpoints)){
                    points(allpoints[[ii]], col = cols[ii], pch=16,type=type)
                }
            }
        }
        if(length(names(allpoints))>0){
            legend("bottomright",legend=names(pp),col=cols,lty = rep(lty, length(pp)),
                   pch=rep(pch,length(pp)))
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
    assert_that(!is.null(cols))
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

##' Gets piecewise mean, given segments
piecewise_mean <- function(y,cp){
  stopifnot(all(1<=cp & cp<=length(y)))
  ## stopifnot(all.equal(sort(cp),cp))
  cp = sort(cp)
  segments = lapply(1:(length(cp)+1), function(ii){ v=c(0,cp,length(y));(v[ii]+1):(v[ii+1]) })
  segment.means = sapply(segments, function(mysegment){mean(y[mysegment])})
  cleanmn = rep(NA,length(y))
  lapply(1:length(segments), function(ii){cleanmn[segments[[ii]]] <<- segment.means[ii]})
  return(cleanmn)
}


##' Helper function for making segment contrasts from a wildBinSeg object OR
##' bsFt object or path object (from genlassoinf package), or just cps
##' @param obj Result from running wbs()
##' @export
make_all_segment_contrasts <- function(obj){

    ## Basic checks
    if(length(obj$cp)==0) stop("No detected changepoints!")
    if(all(is.na(obj$cp)))stop("No detected changepoints!")
    assert_that(!is.null(obj$y))#, msg="in make_all_segment_contrasts(), you need an object that contains the data vector.")

    return(make_all_segment_contrasts_from_cp(obj$cp, obj$cp.sign, length(obj$y)))
}


##' Take |cp| and |cp.sign| and make all segment contrasts.
##' @param cp integer vector of changepoints.
##' @param cp.sign integer vector of changepoint signs.
##' @param n length of data.
##' @export
make_all_segment_contrasts_from_cp <- function(cp, cp.sign, n, scaletype = c("segmentmean", "unitnorm")){

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
##' @param cps if not null, these are the only ones we actually want.
make_all_segment_contrasts_from_wbs <- function(wbs_obj, cps=NULL, scaletype = c("segmentmean", "unitnorm")){

    n = length(wbs_obj$y)
    cp = wbs_obj$cp
    cp.sign = wbs_obj$cp.sign
    winning.intervals = lapply(1:nrow(wbs_obj$results), function(myrow)wbs_obj$results[myrow, c("max.s", "max.e")])

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
        if(length(d)!=n) browser()
    }
    names(dlist) = (cp * cp.sign)

    ## Only return the requested ones
    if(!is.null(cps)){
        dlist = dlist[which(abs(as.numeric(names(dlist))) %in% cps)]
    }
    return(dlist)
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
##' @import glmgen
##' @export
get_piecewise_mean <- function(y, cp){
  if(all.equal(sort(cp), cp)!=TRUE) stop ("cp is not sorted!")
  aug.cp = c(0,cp,length(y))
  segment.inds = sapply(1:(length(cp)+1),
                      function(ii){ (aug.cp[ii]+1):aug.cp[ii+1]})
  mn = rep(NA,length(y))
  for(ind in segment.inds) mn[ind] <- mean(y[ind])
  return(mn)
}


##' Takes a named list of n-length contrast vectors, and filters them so that
##' only the contrasts that are desired i.e. contained in \code{visc}.
##' @export
filter_vlist <- function(vlist, visc=NULL){
    if(!is.null(visc)){
        retain = which(abs(as.numeric(names(vlist))) %in% visc)
        if(length(retain)==0){
            return(data.frame(pvs=NA, locs=NA))
        }
        vlist = vlist[retain]
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
