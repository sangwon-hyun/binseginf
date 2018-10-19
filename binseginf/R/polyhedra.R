#' Generate polyhedra
#'
#' @param obj numeric matrix to be the gamma matrix
#' @param u numeric vector
#' @param ... not used
#'
#' @return object of class polyhedra
#' @export
polyhedra.matrix <- function(obj, u, ...){
  structure(list(gamma = obj, u = u, ...), class = "polyhedra")
}

#' Check if polyhedra is valid
#'
#' @param obj object of class polyhedra
#'
#' @return TRUE if valid
#' @export
is_valid.polyhedra <- function(obj){
  if(!is.numeric(obj$gamma) | !is.matrix(obj$gamma)) {
    print(head(obj$gamma))
    stop("gamma is not a numeric matrix")
  }
  if(!is.numeric(obj$u)) stop("u is not a numeric vector")
  if(nrow(obj$gamma) != length(obj$u)) stop("nrow(gamma) does not match length(u)")

  TRUE
}


#' Generic for functions that combine same-type things.
#' @param obj \code{polyhedra} class object.
combine <- function(obj, ...) {UseMethod("combine")}

##' Combines several polyhedra to a single polyhedron. It is possible to add
##' \code{NULL} valued objects with other polyhedra (to no effect)
##' @param ... Polyhedra objects to add. It is permissible for some of these to
##'     be \code{NULL} valued.
combine.polyhedra <- function(...){

    polyobjs = list(...)
    polyobjs = polyobjs[which(!sapply(polyobjs, is.null))]

    ## NEW: Combine smartly. Handle better!!
    rownums <- sapply(polyobjs, function(a) nrow(a$gamma))
    nonemptyinds = which(rownums!=0)
    if(length(nonemptyinds)==0){
        return(NULL)
    } else {
        starts = cumsum(c(1,rownums[-length(rownums)]))[nonemptyinds]
        ends = cumsum(rownums)[nonemptyinds]
        rowindlist = Map(function(a,b)a:b, starts, ends)
        newgamma = matrix(NA,nrow=sum(rownums), ncol=ncol(polyobjs[[1]]$gamma))
        newu = rep(NA, sum(rownums))
        for(ii in 1:length(nonemptyinds)){
            newgamma[rowindlist[[ii]],] = polyobjs[[nonemptyinds[ii]]]$gamma
            newu[rowindlist[[ii]]] = polyobjs[[nonemptyinds[ii]]]$u
        }
        newpoly = structure(list(gamma = newgamma, u= newu), class = "polyhedra")
        stopifnot(is_valid.polyhedra(newpoly))
        return(newpoly)
    }

}



##' Prints polyhedra
##' @param obj \code{polyhedra} class object.
print.polyhedra <- function(mypoly){
    if(nrow(mypoly$gamma)==0 ){ print("Empty polyhedra object!")
    } else if (all(is.na(mypoly$gamma[1,])) & nrow(mypoly$gamma)==1){
        print("Empty polyhedra object!")
    } else {
        first.n = min(10, nrow(mypoly$gamma), ncol(mypoly$gamma)/2)
        print(paste("Gamma matrix (first", first.n, " rows&cols) looks like:"))
        print(signif((mypoly$gamma[1:first.n, (1:min(first.n*2, ncol(mypoly$gamma)))]),3))
        print(paste("u vector (first", first.n, "entries) looks like:"))
        print(signif(mypoly$u[1:first.n]))
    }
}

## #' Generic for functions that combine same-type things.
## #' @param obj object
## smartadd <- function(obj, ...) {UseMethod("smartadd")}



##' Check if y is in polyhedra (generic).
##' @param obj \code{polyhedra} class object.
##' @export
contained <- function(obj,...){UseMethod("contained")}

##' Check if y is in polyhedra.
##' @param obj \code{polyhedra} class object.
##' @export
contained.polyhedra <- function(obj, y){
    all(obj$gamma %*% y >= obj$u)
}

##' Make empty polyhedra object
## make_empty <- function(obj,...)
make_empty_polyhedra <- function(n){
            emptyrow = rbind(rep(NA,n))[-1,,drop=FALSE]
            return(polyhedra(obj=emptyrow, u=c()))
}


##' Obtains polyhedron from |path| object. A |path| object is created from
##' generalized lasso dual path algorithm in the |genlassoinf| R package.
##' @export
polyhedra.path <- function(obj, numSteps=obj$maxSteps){
    stopifnot(numSteps==obj$maxSteps)
    polyhedra(obj=obj$Gobj.naive$G, u=obj$Gobj.naive$u, nrow.by.step=obj$nkstep)
}


#' Function generic, for a function to get a snapshot of the polyhedron from a
#' particular algorithm step.
#' @param obj \code{polyhedra} class object.
#' @export
snapshot <- function(obj, ...) {UseMethod("snapshot")}

##' Gets a snapshot of the polyhedron from a particular algorithm step.
##' @param obj \code{polyhedra} class object.
##' @param numSteps Algorithm step you want a snapshot from
##' @export
snapshot.polyhedra <- function(obj, numSteps){
    assert_that("nrow.by.step" %in% names(obj))
    inds = (1:obj$nrow.by.step[numSteps])
    return(polyhedra(obj=obj$gamma[inds, ],
                     u=rep(0, length(inds)),
                     nrow.by.step=obj$nrow.by.step))
}


## addpv_fl = addpv.fl

## ##' Helper to harvest polyhedra from FL object.
## polyhedra.fl <- function(obj, numSteps=NULL, record.nrows=TRUE){
##     if(is.null(numSteps)) numSteps = obj$maxsteps
##     Gobj = genlassoinf::getGammat.naive(obj=obj, y=obj$y,
##                                         condition.step=numSteps)
    
##     ## Harvest number of rows per step
##     if(record.nrows){
##         nrow.by.step = obj$nkstep
##     } else {
##         nrow.by.step = NULL
##     }

##     poly = polyhedra(obj=Gobj$G, u=Gobj$u, nrow.by.step=nrow.by.step)
##     return(poly)
## }

## polyhedra_fl = polyhedra.fl


##' Proprietary print object for |path| class object. This is temporary, and
##' assumes that fused lasso (and not a different form of generalized lasso) is
##' run.
print.fl <- function(obj){
    cat("Detected changepoints using FL with", obj$numSteps, "steps is",
        obj$cp * obj$cp.sign, fill=TRUE)
    if(!is.null(obj$pvs)){
        cat("Pvalues of", names(obj$pvs), "are", obj$pvs, fill=TRUE)
    }
}
