#' Creates halfspaces for the main maximize-max-cusums.
##' @param obj Result of running WBS.
max_of_max_obj <- function(obj){

    ## Basic check
    stopifnot(is_valid.wbsFt(obj))

    ## Apply max to all |semat|'s in |obj|
    polyobjs = lapply(obj$tree, max_of_max_semat, obj)

    ## Combine and return
    do.call(combine.polyhedra, polyobjs)
}


##' Creates halfspaces for the main maximize-max-cusums.  TODO: many of the rows
##' overlap; keep track of this and don't add duplicates. Currently, this is not
##' done carefully here; but rather, during aggregation of the polyhedra, the
##' duplicates are eliminated.
##' @param semat A matrix with results from a node of an obj$tree.
max_of_max_semat <- function(semat, obj){
    ## Basic checks
    ## stopifnot(is_valid.semat(semat))

    intervals = obj$intervals
    max.ind = semat[,"maxhere"]
    max.m = semat[max.ind,"m"]
  if(max.m==0){
      max.s = semat[max.ind, "s"]
      max.b = semat[max.ind, "b"]
      max.e = semat[max.ind, "e"]
  } else {
    max.s = intervals$se[[max.m]][1]
    max.b = semat[max.ind,"b"]
    max.e = intervals$se[[max.m]][2]
    }

    ## 1. First, characterize the sign of the max.cusum.contrast
    max.cusum.contrast = unsigned_contrast(max.s, max.b, max.e, y=obj$y)
    newrows1 = rbind(max.cusum.contrast)
    newu1 = rep(0,nrow(newrows1))


    ## 2. Second, Compare the max-cusum of the grand max
    M = semat[,"m"]
    if(length(M)==0)  return(NULL)
    newrows2 = do.call(rbind, lapply(M, function(m){

      if(m==0){
        se = c(max.s, max.e)
        b = max.b
      } else {
        se = intervals$se[[m]]
        b = semat[semat[,"m"]==m,"b"]
      }
        s.to.e = (se[1]:se[2])
        other.bs = s.to.e[-which(s.to.e == se[2])]
        if(m==max.m) other.bs = other.bs[other.bs!=max.b]
        if(length(other.bs)==0) return(rbind(rep(NA,length(obj$y)))[-1,])

        ## Subtract all other contrast from the maximum cusum contrast
        other.cusum.contrasts = do.call(rbind, lapply(other.bs, function(other.b){
            signed_contrast(se[1], other.b, se[2], y=obj$y)}))
        subtracted.contrasts = rbind(sweep(-rbind(other.cusum.contrasts), 2,
                                           max.cusum.contrast, "+" ),
                                     sweep(+rbind(other.cusum.contrasts), 2,
                                           max.cusum.contrast, "+" ))
        if(ncol(subtracted.contrasts)!=length(obj$y)) subtracted.contrasts = t(subtracted.contrasts)
        return(subtracted.contrasts)
    }))
    newu2 = rep(0,nrow(newrows2))


    ## 3. Check threshold exceedance:
    if(all(semat[,"passthreshold"]==FALSE)){
        newrows3 = -rbind(max.cusum.contrast)
        newu3 = -obj$thresh
    } else {
        newrows3 = rbind(max.cusum.contrast)
        newu3 = obj$thresh
    }


    newrows = rbind(newrows1, newrows2, newrows3)
    newu = c(newu1, newu2, newu3)

    ## If newrows is empty (no comparisons to be made), then don't do anything
    if(length(as.numeric(newrows))==0){ return(NULL)}

    ## Return it as a polyhedron
    return(polyhedra.matrix(obj = newrows, u = newu))
}


#' Generate polyhedra matrix from wbs output
#' Forms both Gamma matrix and u vector
#'
#' @param obj Output from wbs
#' @param ... not use now
#'
#' @return An object of class polyhedra
#' @export
polyhedra.wbsFt <- function(obj, ...){

    ## Basic checks
    stopifnot(is_valid.wbsFt(obj))

    ## Combine polyhedron
    return(max_of_max_obj(obj))
}

##' Temporary function to check if it is correct; NOT to be called at
##' runtime, but only at test time or internally.
##' @param poly An object of class polyhedra
##' @param y A data vector with the appropriate length, used when
##'     creating this polyhedron.
check_polyhedra <- function(poly, y){
    ## print(poly$gamma%*%y >= poly$u)
    stopifnot(all(poly$gamma %*% y >= poly$u))
}


