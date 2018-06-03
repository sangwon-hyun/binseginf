##' Function that takes changepoint location vector \code{cp} (in the order that
##' they entered, from fixed step SBS of fixed step WBS), and calcualtes BIC of
##' each intermediate changepoint model (piecewise constant underlying mean
##' Gaussian model).
##' @param cp Vector of changepoints
##' @param y Data vector.
##' @param sigma Noise standard deviation.
##' @param verbose \code{TRUE} to print progress
##'
##' @return A list containing two elements: a numeric vector of BICs, and a list
##'     containing basis vectors \eqn{a} of the projection \eqn{P=aa^T}from one
##'     changepoint column basis space to the next.
get_ic <- function(cp, y, sigma, consec=2, maxsteps=length(cp), type="bic", verbose=FALSE){
    ## Basic checks
    if(type!="bic") stop("Only BIC is coded so far!")

    ## Collect things
    n = length(y)
    D = dual1d_Dmat(n)
    ic = pen = RSS = rep(NA, maxsteps)
    resid = list()

    ## Collect BIC at each step 0 ~ (maxsteps-1)
    allsteps = 0:(pmin(maxsteps,length(cp)) )
    for(ii in allsteps ){

        if(verbose)  cat('step', ii, '\n')

        ## Form proj null(D_{-B}) by orth proj onto row(D_{-B}) = col(t(D_{-B})) ~= tD
        if(ii==0){
            tDb = rep(1,n)
            curr.proj = .proj(tDb)
        } else {
            tDb = make.tDb(cps = cp[1:ii], n=n)
            curr.proj = .proj(tDb)
        }
            ## y.fitted = rep(mean(y),n)
            y.fitted = (curr.proj) %*% y

        ## Obtain RSS and penalty
        myRSS = sum( (y - y.fitted)^2 )
        ## mydf  = n-rr
        mydf = ii + 1
        if(ii==0) prev.df = mydf
        mypen = (sigma^2) * mydf * log(n)

        ## Obtain (2norm-scaled) residual projection vectors
        if(ii==0){
            myresid = rep(NA, n)
        } else {
            myresid = get_basis(cps_old=(if(ii==1) c() else cp[1:(ii-1)]),
                                cp_new=cp[ii], n=n)
        }

        ## Store BIC and re
        ic[ii+1] <- myRSS + mypen
        pen[ii+1] <- mypen
        RSS[ii+1] <- myRSS
        resid[[ii+1]] <- myresid

        ## Update things for next step.
        prev.proj = curr.proj
    }
    names(ic) = names(pen) = names(RSS) = names(resid) = allsteps

    ## Flag the result
    flag = ic_flag(ic, consec, cp, maxsteps)

    ## Assign appropriate stop time, based on the flag
    if(flag == "zero.stop"){
        stoptime = 0
    } else if (flag == "not.enough.steps") {
        stoptime = NA
    } else if (flag == "didnt.stop") {
        stoptime = NA
    } else if(flag == "normal"){
        stoptime = .whichrise(ic,consec) - 1
        stoptime = pmin(stoptime, length(y)-consec-1)
    }

    ## Get directions
    seqdirs = c(.getorder(ic))
    names(seqdirs) = allsteps

    if(flag=="normal"){
        ## Get order of ICs

        ## Make empty things before collecting halfspaces
        newrows = matrix(NA, nrow = 2*(stoptime+consec+1),
                         ncol = length(y))
        newu = rep(NA, 2*(stoptime+consec+1))
        irow = 0

        ## Collect halfspaces
        all.relevant.steps = (1:(stoptime + consec))-1
        for(ii in all.relevant.steps){

            residual = resid[[toString(ii+1)]]
            const = pen[toString(ii+1)] - pen[toString(ii)]
            pos.resid = sign(as.numeric(t(residual)%*%y))*residual/sqrt(sum(residual^2))

            if(seqdirs[toString(ii+1)] < 0){
                ## Add one row \sqrt{C} < z_a \times a^Ty
                newrows[irow+1,] = pos.resid #(sign(as.numeric(t(residual)%*%y)) * residual)/sqrt(sum(residual^2))
                newu[irow+1] = sqrt(const)
                irow = irow + 1
            } else {
                ## Add two rows -\sqrt{C} < z_a \times a^Ty < \sqrt{C}
                newrows[irow+1,] = pos.resid # (sign(as.numeric(t(residual)%*%y)) * residual) / sqrt(sum(residual^2))
                newrows[irow+2,] = -pos.resid #(-sign(as.numeric(t(residual)%*%y)) * residual)/sqrt(sum(residual^2))
                newu[irow+1] = -sqrt(const)
                newu[irow+2] = -sqrt(const)
                irow = irow + 2
            }
        }
        ## Form polyhedra
        poly = polyhedra(obj = trim(newrows), u = trim(newu))
    } else {
        poly = make_empty.polyhedra()
    }

    return(structure(list(ic=ic, consec=consec, resid=resid, pen=pen, RSS=RSS,
                          stoptime=stoptime,y=y, type=type, flag=flag, poly=poly, seqdirs=seqdirs), class="ic"))
}

ic_flag <- function(ic, consec=2, cp, maxsteps){

    ## Flag the result, case by case.
    if(length(ic) < consec+1){
        ##  warning("Not enough steps to do forward sequential BIC/AIC!")
        return("not.enough.steps")
    } else if (.whichrise(ic,consec) + consec - 1 > pmin(maxsteps, length(cp))){
        ##  warning("Didn't stop!")
        return("didnt.stop")
    } else if (.whichrise(ic,consec) - 1 == 0){
        return("zero.stop")
    } else {
        return("normal")
    }
}



##' Does two things: Flags the result of the IC path, and locates first point of
##' \code{consec} rises in IC path
##' @param ic numeric vector containing information criteria.
##' @return First point at which ic rises.
##' @export
.whichrise <- function(ic, consec = 2){

    ind = 1
    done = FALSE
    while(ind < (length(ic)-consec+1) ){
        ictol = 1E-10
        if( all(ic[(ind+1):(ind+consec)] > ic[(ind):(ind+consec-1)] + ictol )) break
        ind = ind+1
    }
    first.point.of.rise = pmin(pmax(ind,1),length(ic))
    return(first.point.of.rise=first.point.of.rise)
}


##' is_valid for ic objects
##' @return TRUE if valid
is_valid.ic <- function(obj){
    if(!(all(names(obj)%in% c("ic","consec","resid","pen","type","stoptime","RSS",
                              "y", "flag") ))){
        stop("obj needs to be produced from get_ic()")
    }
    TRUE
}



## Helper function to project on /column/ space of matrix.
.proj <- function(mymat){
    return(mymat %*% solve(t(mymat)%*%mymat, t(mymat)))
}

## Returns a sequence of +1 and -1 for sequential incr and decrements in a
## vector; assume first step always dips; _almost_ always true
## @param ic a numeric vector
.getorder <- function(ic){
    seqdir = c(NA, sign(ic[2:(length(ic))] - ic[1:(length(ic)-1)]))
    names(seqdir) = 0:(length(ic) - 1)
    return(seqdir)
}



##' Helper function to make basis vectors for the piecewise constant space
##' broken at changepoints at |cps|. |tDb| just stands for D^T[,-b].
make.tDb <- function(cps=c(), n){
    cps = sort(cps)
    if(length(cps)==0){
        return(cbind(rep(1/n,n)))
    } else {
        ncp = length(cps)
        tDb = matrix(0, nrow=n, ncol=ncp+1)
        cps.aug = c(0,cps,n)
        for(ii in 1:(ncp+1)){
            nonzero.inds = (cps.aug[ii]+1):(cps.aug[ii+1])
            tDb[nonzero.inds, ii] = 1/length(nonzero.inds)
        }
    }
    return(tDb)
}

##' Get basis vector of the residual subspace between the two spanned by
##' \code{c(cps_old)} and \code{c(cps_old, cp_new)}.
get_basis <- function(cps_old, cp_new, n){
    if(length(cps_old)!=0){
        sorted_cps_old = sort(cps_old)
    } else {
        sorted_cps_old=cps_old
    }
    sorted_cps_old = c(0,sorted_cps_old,n)
    imin = which.min(cp_new > sorted_cps_old)
    cp_right = sorted_cps_old[imin]
    cp_left = sorted_cps_old[imin-1]
    basisvec = make_all_segment_contrasts_from_cp(cp = c(cp_left, cp_new, cp_right),
                                                  cp.sign = rep(1,3),
                                                  n=n, scaletype = "unitnorm")[[toString(cp_new)]]
    return(basisvec)
}


