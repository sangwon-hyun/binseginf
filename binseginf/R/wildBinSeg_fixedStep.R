##' Fixed step WBS function. There is a ``mimic'' option that helps you create
##' new polyhedron that says that the same |numSteps| winners win, but among a
##' new set of intervals.
##' @param intervals Optionally, add an |intervals| class object that contains
##'     intervals.
##' @param inference.type Specify whether you want to actually accumulate the
##'     halfspaces, or just precalculate G*v and G*y for direct calculation of
##'     p-values.
##' @param stop.time Manually inputted stop time. You can obtain a data
##'     dependent stop time using local IC rises, via \code{ic_wrapper()}.
##' @param ic.poly polyhedron for the ic stopping event
##' @param v contrast vector; used only for ic.poly
##' @return List of useful algorithm results, including the |gamma| matrix and
##'     |u|.
##' @export
wildBinSeg_fixedSteps <- function(y, numSteps, numIntervals=NULL,
                                  intervals=NULL, mimic=FALSE, wbs.obj=NULL,
                                  comprehensive=FALSE, cumsum.y=NULL, cumsum.v=NULL,
                                  inference.type =c("rows","pre-multiply", "none"),
                                  stop.time=numSteps,
                                  ic.poly = NULL,
                                  v=NULL){

    inference.type = match.arg(inference.type)

    ## Basic checks
    n = length(y)
    assert_that(!is.null(numIntervals) | !is.null(intervals),
                msg="Must provide either |numIntervals| or |intervals|.")
    if(mimic & is.null(wbs.obj))stop("When mimic=TRUE, provide a wbs object.")
    if(!mimic & !is.null(wbs.obj))stop("Even though mimic==FALSE, you're providing a wbs object.")

    ## Calculate all cusums for all intervals
    if(is.null(intervals)){
        intervals = intervals(numIntervals=numIntervals, n= n,
                              comprehensive=comprehensive)
    }
    intervals = addcusum(intervals=intervals, y=y)

    ## Make empty array to store selection events
    results = matrix(NA, nrow=numSteps, ncol = 6)
    colnames(results) = c("max.s", "max.b", "max.e", "max.sign", "scope.s", "scope.e")

    ## Initially, everything qualifies
    qual.inds = 1:nrow(intervals$cusummat)
    scope = rbind(c(1,n))
    colnames(scope)  = c("s", "e")
    new.rows = new.info = list()

    ## Pre-form IC halfspaces, if applicable
    if(!is.null(ic.poly)){
        Gy.ic = ic.poly$gamma %*% y
        Gv.ic = ic.poly$gamma %*% v
        u.ic = ic.poly$u
    } else {
        Gy.ic = Gv.ic = u.ic = c()
    }

    ## Iterate
    istep = 1
    time.to.stop = FALSE
    if(mimic) stopstep = (min(c(numSteps,stop.time, nrow(wbs.obj$results))))
    if(!mimic) stopstep = min(c(numSteps,stop.time))
    while(istep <= stopstep & !time.to.stop){

        ## Do maximization in the respective intervals.
        if(mimic){
            max.info = .get_max_info(wbs.obj, istep)
        } else {
            max.info = maximize.intervals(intervals=intervals, qual.inds=qual.inds)
        }
        if(inference.type=="rows"){
            ## Collect halfspaces for the model selection event
            new.rows[[istep]] = form_rows(intervals=intervals,
                                          max.s=max.info$max.s,
                                          max.b=max.info$max.b,
                                          max.e=max.info$max.e,
                                          max.sign=max.info$max.sign,
                                          qual.inds=qual.inds, n=n,
                                          y=(if(!mimic) y else NULL))
        } else if (inference.type=="pre-multiply"){

            ## New routine for new.rows, that collects Gy and Gv instead of
            ## actual rows
            new.info[[istep]] = form_info(intervals=intervals,
                                          max.s=max.info$max.s,
                                          max.b=max.info$max.b,
                                          max.e=max.info$max.e,
                                          max.sign=max.info$max.sign,
                                          qual.inds=qual.inds,
                                          cumsum.y=cumsum.y, cumsum.v=cumsum.v,
                                          Gy.ic=Gy.ic, Gv.ic=Gv.ic, u.ic)
        } else {
            new.info = c()
        }

        ## Calculate the maximizing scope
        max.scope.irow = which(apply(scope, 1,function(myscope){
            (myscope["s"]<= max.info$max.b) & (max.info$max.b <= myscope["e"])
        }))
        max.scope = scope[max.scope.irow,]

        ## Update qualifying indices in cusummat
        old.qual.inds = restrict(intervals=intervals,  s=max.scope["s"], e=max.scope["e"])
        new.qual.inds1 = restrict(intervals=intervals, s=max.scope["s"], e=max.info$max.b)
        new.qual.inds2 = restrict(intervals=intervals, s=max.info$max.b+1, e=max.scope["e"])
        qual.inds = qual.inds[!(qual.inds %in% old.qual.inds)]
        qual.inds = c(qual.inds, new.qual.inds1, new.qual.inds2)

        ## Update scopes
        scope = scope[-max.scope.irow,,drop=FALSE]
        scope = rbind(scope,
                      c(max.scope["s"], max.info$max.b),
                      c(max.info$max.b + 1, max.scope["e"]))
        too.short = which(scope[,"e"] - scope[,"s"] < 1)
        if(length(too.short)>=1){
            scope = scope[-too.short,,drop=FALSE]
        }

        ## Record some results
        results[istep,] <- c(max.info$max.s, max.info$max.b, max.info$max.e, max.info$max.sign,
                             max.scope["s"], max.scope["e"])

        ## Break if done.
        time.to.stop = (nrow(scope)==0 | length(qual.inds)==0)
        istep = istep+1
    }
    stoptime = istep-1
    ## if(stoptime != numSteps) print("Stopped early!")
    results = results[1:stoptime,,drop=FALSE]

    ## Aggregate things for inference
    if(inference.type=="rows"){
        gamma=do.call(rbind, new.rows)
        names(new.rows) = paste0("step-",1:stoptime)
        Gy=NULL
        Gv=NULL
        u = rep(0,nrow(gamma))
    } else if (inference.type == "pre-multiply"){
        Gy = do.call(c,lapply(new.info, function(a)a[["Gy"]]))## unlist()
        Gv = do.call(c,lapply(new.info, function(a)a[["Gv"]]))
        u = do.call(c,lapply(new.info, function(a)a[["u"]]))
        gamma = NULL
        ## u = rep(0, length(Gy))
    } else {
        Gy=Gv=u=gamma=NULL
    }

    return(structure(list(results=results, gamma=gamma, u=u,
                          Gy=Gy, Gv=Gv,
                          rows.list = new.rows,
                          info.list = new.info,
                          cp=results[,"max.b"], cp.sign=results[,"max.sign"],
                          numSteps=numSteps,
                          mimic=mimic,
                          intervals=intervals,
                          numIntervals=intervals$numIntervals,
                          inference.type=inference.type,
                          y=y), class="wbsFs"))
}



##' Print function for convenience, of |wbs| class object.
##' @export
print.wbsFs <- function(obj){
    if(obj$mimic) cat("Mimicked object!", fill=TRUE)
    cat("Detected changepoints using WBS with", obj$numSteps, "steps is", obj$cp * obj$cp.sign, fill=TRUE)
}

##' Checks if obj is a valid |wbs| class object
is_valid.wbsFs <- function(obj){
    return(all(names(obj) %in% c("results", "gamma", "u", "cp", "cp.sign", "y",
                                 "numSteps", "mimic", "rows.list","intervals",
                                 "numIntervals", "Gy", "Gv", "info.list","inference.type")))
}

.get_max_info <- function(wbs.obj,...){ UseMethod(".get_max_info")}
##' Helper to /manually/ obtain maximizing information from a |wbs| class object, at i'th step
.get_max_info.wbsFs <- function(wbs.obj, istep){
    assert_that(is_valid.wbsFs(wbs.obj))
    myrow  = wbs.obj$results[istep,]
    mylist = lapply(myrow, function(a) a)
    mylist = mylist[names(mylist) %in% c("max.s", "max.b", "max.e", "max.sign")]
    mylist[["max.i"]] = NA
    return(mylist)
}


##' Polyhedra export function.
##' @param obj Object of class |wbsFs|
##' @return object of class |polyhedra|.
##' @export
polyhedra.wbsFs <- function(obj){
    polyhedra(obj = obj$gamma, u = obj$u)
}
