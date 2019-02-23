## Synposis: temporary location for helpers for the script more-power-helpers.R


##' Calculates full-model test contrasts from two-step binseg.
##' @param y data vector
##' @return list of vectors, whose names equal to changepoint locations.
two_step_contrasts <- function(y){
    obj = binseginf::bsfs(y, numSteps=2)
    vlist = lapply(obj$cp, function(cp){
        binseginf::make_all_segment_contrasts_from_cp(cp, +1, n)
    })
    names(vlist) = obj$cp
    return(vlist)
}


##' Calculates a p-value from a list of endpoints.
##' @param v contrast vector
##' @param vslist list of endpoints along v
##' @param y data vector
pv_from_endpoints <- function(v, vslist, y){

    vty = sum(v*y)
    all.up = which(sapply(vslist, function(pair){ all(pair>vty) }))
    all.in = which(sapply(vslist, function(pair){
        pair[1] < vty & vty < pair[2]}))
    if(length(all.in)==0){
        ## This needs to be maanged more carefully
        return(c())
    }
    all.below = which(sapply(vslist, function(pair){all(pair<vty)}))


    ## Calculate the masses
    getmass <- function(pair){
        pair = as.numeric(pair)
        getmass_from_pair(v, pair)}

    ## Get probability masses
    all.masses = lapply(vslist, getmass)
    in.mass = getmass(c(vty, vslist[[all.in]][2]))
    upper.masses = lapply(vslist[all.up], getmass)

    ## Calculate p-value
    pv = (sum(unlist(upper.masses)) + in.mass)/sum(unlist(all.masses))
    return(pv)
}

##' Calculates probability mass between a pair of endpoints, in distribution of
##' vty assuming y is from normal with zero mean and sd=1.
##' @param pair is the couplet of endpoints along contrast vector v.
##' @return mass between pair[1] and pair[2], in distribution of vty.
getmass_from_pair <- function(v, pair, sigma=1){
    stopifnot(all.equal(sum(v*v), 1))
    mass = pnorm(pair[2]) - pnorm(pair[1])
    return(mass)
}



##' Copied over from binseginf package, file name bsfs.R.
.cusum_contrast_full <- function(start, idx, end, n){
  res <- rep(0, n)
  res[start:end] <- .cusum_contrast(start, idx, end)

  res
}

##' Returns a gamma matrix for |b| winning in the endpoints |s| and |e|.
##' @param s start
##' @param b break
##' @param e end
##' @param sign.win Sign of winner
##' @return matrix
comparison <- function(s, b, e, n, sign.win){
    if(s+1==e) return(rbind(rep(Inf, n))[-1,])

    ## Make the matrix of competitors
    mat = matrix(ncol=3, nrow=e-s-1)
    mat[,1] = rep(s, e-s-1)
    breaks = seq(from=s, to=e-1)
    breaks = breaks[which(breaks!=b)]
    mat[,2] = breaks
    mat[,3] = rep(e, e-s-1)

    win.contrast <- .cusum_contrast_full(s, b, e, n)
    lose.contrast <- t(apply(mat, 1, function(x){
        .cusum_contrast_full(x[1], x[2], x[3], n)
    }))

    ## add inequalities to compare winning split to all other splits
    res <- .vector_matrix_signedDiff(win.contrast, lose.contrast, sign.win,
                                     rep(1, nrow(lose.contrast)))
    res2 <- .vector_matrix_signedDiff(win.contrast, lose.contrast, sign.win,
                                      -rep(1, nrow(lose.contrast)))
    newrows = rbind(res, res2)
    return(newrows)
}

##' Returns a gamma matrix for CUSUM winning among all the CUSUMs from the
##' endpoints |s| and |e|.
##' @param s start of competitor.
##' @param e end of competitor.
##' @param sbe.win start, break and end of winner.
##' @param n data length
##' @param sign.win Sign of winner
##' @return matrix
comparison_external <- function(s, e, sbe.win, n, sign.win){

    if(s==e) return(rbind(rep(Inf, n))[-1,])##return(NULL)

    ## Make the matrix of competitors
    mat = matrix(ncol=3, nrow=e-s)
    mat[,1] = rep(s, e-s)
    mat[,2] = c(s:(e-1))
    mat[,3] = rep(e, e-s)

    win.contrast <- .cusum_contrast_full(sbe.win[1], sbe.win[2], sbe.win[3], n)
    lose.contrast <- t(apply(mat, 1, function(x){
        .cusum_contrast_full(x[1], x[2], x[3], n)
    }))

    ## add inequalities to compare winning split to all other splits
    res <- .vector_matrix_signedDiff(win.contrast, lose.contrast, sign.win,
                                     rep(1, nrow(lose.contrast)))
    res2 <- .vector_matrix_signedDiff(win.contrast, lose.contrast, sign.win,
                                      -rep(1, nrow(lose.contrast)))
    newrows = rbind(res, res2)
    return(newrows)
}



##' Obtain binseg's first step's endpoints along v.
##' @param n data length.
##' @param j tested location.
##' @param j.win.sign sign of tested location.
##' @return A single gamma matrix
first_step <- function(n, j, j.win.sign){

    ## Collect all halfspaces
    newrows = comparison(1, j, n, n, j.win.sign) ## Test this function!!
    newrows.list = list(newrows)
    names(newrows.list) = paste0("second-step-", j)

    return(newrows)
}

##' Obtain binseg's second step's polyhedra along v. This is a union (list) of
##' gamma matrices.
##' @param n data length.
##' @param j tested location.
##' @param j.win.sign sign of tested location.
##' @return list of couplets of endpoints.
second_step <- function(n, j, j.win.sign){
    
    ## Helper for internally calculating second step rows given first step
    ## detection in first step |cand|.
    secondstep <- function(cand){
        ## if(cand==1) return(rbind(rep(Inf, n))[-1,])#return(NULL)
        if(cand <= j){
            newrows.left = comparison_external(1, cand, c(cand+1,j,n), n,
                                               j.win.sign) 
            newrows.right = comparison(cand+1, j, n, n, j.win.sign)
        } else {
            newrows.left = comparison(1, j, cand, n, j.win.sign)
            newrows.right = comparison_external(cand+1, n, c(1,j,cand), n,
                                                j.win.sign) 
        }
        all.newrows = rbind(newrows.left, newrows.right)
        if(ncol(all.newrows)!=n) browser()
        return(all.newrows)
    }

    ## Collect a union of polyhedra
    first.step.candidates = (1:(n-1))[-j]
    newrows.list = lapply(first.step.candidates, secondstep)

    ## Intersect each with the characterization of $ii$.
    newrows.list2 = list() ## Second steps' rows
    for( ii in 1:length(first.step.candidates )){
        cand = first.step.candidates[[ii]]
        secondstep.newrows = newrows.list[[ii]]
        cand.row.list = lapply(c(-1,1), function(cand.sign){
            firststep.newrows = first_step(n, cand, cand.sign)
            combined.newrows = rbind(firststep.newrows, secondstep.newrows)
            rownames(combined.newrows) = c(rep(1,nrow(firststep.newrows)),
                                           rep(2,nrow(secondstep.newrows)))
            return(combined.newrows)
        })
        names(cand.row.list) = c(-1,1)
        newrows.list2[[ii]] = cand.row.list
    }
    names(newrows.list2) = paste0("step2-cand", first.step.candidates)
    return(newrows.list2)
}



##' Helper to get vup and vlo.
##' @param y data vector.
##' @param v contrast vector.
##' @param newrows gamma matrix.
##' @param sigma noise level.
getvuplo <- function(y, v, newrows, sigma=1){
    poly = binseginf::polyhedra(newrows, u=rep(0, nrow(newrows)))
    obj = poly.pval2(y, poly, v, sigma, bits=5000)
    return(c(vlo=obj$vlo,vup=obj$vup))
}

##' Gets union of list of endpoint pairs.
##' @param all.pairs.mat 2-colum matrix of pairs of endpoints.
##' @return 2-row data frame containing interval endpoints.
get_union <- function(all.pairs.mat){
    
    ## all.pairs.mat = do.call(rbind, all.pairs.list)
    all.pairs.mat = all.pairs.mat[which(all.pairs.mat[,"vlo"] <= all.pairs.mat[,"vup"]),]
    idf = Intervals(all.pairs.mat)
    mat = as.data.frame(interval_union(idf))
    vslist = lapply(1:nrow(mat),function(irow){
        as.matrix(mat)[irow,]
    })
    return(vslist)
}

##' Makes spike (+-1 (default), or +-l at changepoint location) contrast.
##' @param n data length.
##' @param cp changepoint location.
##' @param cp.sign changepoint sign.
spike_contrast<- function(n, cp, cp.sign, l=1){
    v = rep(0,n)
    if(l==1){
        v[cp+1] = +1
        v[cp] = -1
    } else {
        ind1 = cp+(1:l)
        ind1 = ind1[which(ind1%in%(1:n))]
        ind2 = cp+1-(1:l)
        ind2 = ind2[which(ind2%in%(1:n))]
        v[ind1] = +1/length(ind1)
        v[ind2] = -1/length(ind2)
    }
    v = v/sqrt(sum(v*v))
    v = v*cp.sign
    return(v)
}

##' Construct union of halfspaces intersections given y and v.
##' @param y data vector.
##' @param v contrast vector.
##' @param cp changepoint location.
##' @param cp.sign changepoint sign.
get_union_of_pairs <- function(y, v, cp, cp.sign){

    n = length(y)

    ## Get all pairs, from 1st and 2nd step.
    first.pair = getvuplo(y, v, first_step(n, cp, cp.sign))
    second.pairs.list = lapply(second_step(n, cp, cp.sign), lapply, function(newrows){
        return(getvuplo(y, v, newrows))})

    ## Convert to matrices and take union
    second.pairs.mat = do.call(rbind, unlist(second.pairs.list, recursive=FALSE))
    first.pairs.mat = rbind(first.pair)
    all.pairs.mat = rbind(first.pairs.mat, second.pairs.mat)
    ## print(all.pairs.mat)
    all.union = get_union(all.pairs.mat)
    ## print(all.union)

    return(all.union)
}


##' Simulation driver for one-jump simulations.
##' @param nsim Number of simulations.
##' @param n Sample size.
##' @param l Size of window from which spike tests are to be form formed.
dosim <- function(nsim, lev, n=10, mc.cores=4, l=1, verbose=TRUE, sigma.add=0.2){

    start.time = Sys.time()
    pvs.list = mclapply(1:nsim, function(isim){
        if(verbose) printprogress(isim,nsim, start.time=start.time)
        mn = onejump(lev, n)
        y = mn + rnorm(n)
    
        ## Original contrasts
        out = bsfs(y, 2)
        vlist.orig = Map(function(cp,cp.sign){
            spike_contrast(n, cp, cp.sign, l=l)
        }, out$cp, out$cp.sign)
        names(vlist.orig) = out$cp * out$cp.sign

        ## Noisy contrasts
        out.noisy = bsfs(y, 2, sigma.add=sigma.add)
        vlist.noisy = Map(function(cp,cp.sign){
            spike_contrast(n, cp, cp.sign, l=l)
        }, out.noisy$cp, out.noisy$cp.sign)
        names(vlist.noisy) = out.noisy$cp * out.noisy$cp.sign
        
        ## Calculate TG p-values
        pvs.tg = addpv(out, sigma=1, vlist=vlist.orig)$pvs
        pvs.rand.tg = addpv(out.noisy, sigma=1, type="addnoise", sigma.add=sigma.add, vlist=vlist.noisy)$pvs
        pvs.tg.segment = addpv(out, sigma=1)$pvs
        pvs.rand.tg.segment = addpv(out.noisy, sigma=1, type="addnoise", sigma.add=sigma.add)$pvs
    
        ## Calculate TZ p-values
        pvs.tz = c()
        for(ii in 1:2){
            cp = out$cp[ii]
            cp.sign = out$cp.sign[ii]
            ## v = make_all_segment_contrasts_from_cp(cp=cp, cp.sign=cp.sign, n=n,
            ##                                        scaletype="unitnorm")[[1]]
            v = spike_contrast(n, cp, cp.sign, l=l)
            all.union = get_union_of_pairs(y, v, cp, cp.sign)
            pv = pv_from_endpoints(v, all.union, y)
            if(length(pv)==0) next
            pvs.tz[ii] = pv
        }
        names(pvs.tz) = names(vlist.orig)
        ## pvs.mat = rbind(pvs.orig, pvs.new)
        ## names(pvs.mat) = out$cp * out$cp.sign
        ## names(pvs.rand.orig) = out.noisy$cp * out.noisy$cp.sign
        ## return(list(orig=pvs.orig, orig.rand=pvs.rand.orig, new=pvs.new))
        return(list(tg=pvs.tg, tg.rand=pvs.rand.tg, tz=pvs.tz, tg.segment=pvs.tg.segment, tg.rand.segment=pvs.rand.tg.segment))
        ## return(pvs.mat)
    }, mc.cores=mc.cores)
    return(pvs.list)
}


##' Temporary simulation driver for n=10 simulations.
dosim_randomized <- function(nsim, lev, n=10, mc.cores=4){

    pvs.list = mclapply(1:nsim, function(isim){
        printprogress(isim,nsim)
        mn = onejump(lev, n)
        y = mn + rnorm(n)
        out = bsfs(y, 2, sigma.add=.5)
    
        ## Original
        pvs.rand.orig = addpv(out, sigma=1, type="addnoise", sigma.add=.5)$pvs
        names(pvs.rand.orig) = out$cp * out$cp.sign
        return(pvs.rand.orig)
    }, mc.cores=mc.cores)
    return(pvs.list)
}


