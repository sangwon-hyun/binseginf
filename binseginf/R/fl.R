##' Wrapper for running FL.
##' @param sigma.add is important
##' @return list of information regarding the fitted algorithm. The list
##'     component \code{y} is the data used for actual fitting; \code{y.orig} is
##'     the pre-noise original data; \code{y.addnoise} (if not null) is the
##'     added noise.
##' @export
fl <- function(y, numSteps, sigma.add=NULL, ic.stop=FALSE, maxSteps=NULL, numIntervals=NULL){

    ## Basic checks
    if(ic.stop){
        ## warning("|numSteps| option will be ignored and replaced with maxSteps.")
        assert_that(!is.null(maxSteps))
        numSteps = maxSteps
    }    
    if(numSteps >= length(y)) stop("numSteps must be strictly smaller than the length of y")
    if(numSteps <= 0) step("numSteps must be at least 1.")
    if(!is.null(numIntervals)) warning("You provided |numIntervals| but this will not be used.")
    

    y.orig = y
    if(!is.null(sigma.add)){
        y.addnoise = rnorm(length(y), 0, sigma.add)
        y = y + y.addnoise
    }
    obj = genlassoinf::dualpathSvd2(y, maxsteps=numSteps,
                                    D=genlassoinf::makeDmat(length(y), ord=0))
    obj$numSteps = numSteps
    obj$y.orig = y.orig

    ## If applicable, collect IC stoppage
    obj$ic.stop = ic.stop
    if(ic.stop){

        ## Obtain IC information
        ic_obj = get_ic(obj$cp, obj$y, 2, sigma)
        obj$stoptime = ic_obj$stoptime
        obj$consec = consec
        obj$ic_poly = ic_obj$poly
        obj$ic_flag = ic_obj$flag
        obj$ic_obj = ic_obj

        ## Update changepoints with stopped model
        obj$cp.all = obj$cp
        obj$cp.sign.all = obj$cp.sign
        if(ic_obj$flag=="normal"){
            obj$cp = obj$cp[1:obj$stoptime]
            obj$cp.sign = obj$cp.sign[1:obj$stoptime]
        } else {
            obj$cp = obj$cp.sign = c()
        }
    }

    if(!is.null(sigma.add)){
        obj$noisy = TRUE
        obj$sigma.add = sigma.add
        obj$y.addnoise = y.addnoise
    }

    ## Make into FL class here.
    class(obj) = c("fl", "path")

    return(obj)
}



##' Class checker.
is.fl <- function(x) inherits(x, "fl")
