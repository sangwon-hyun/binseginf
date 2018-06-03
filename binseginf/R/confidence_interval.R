#' TG confidence interval from post-selection inference. Assumes that one-sided
#' test is in the right-sided direction i.e. H_0: v^T mu>0.
#'
#' @param y numeric vector
#' @param polyhedra polyhedra object
#' @param contrast contrast numeric vector
#' @param sigma numeric to denote the sd of the residuals
#' @param alpha numeric between 0 and 1 with default of 0.95. This is the
#'     significance level, guaranteeing that alpha percentage of the intervals
#'     will cover the true parameter.
#' @param gridsize numeric to denote how fine of a grid to invert the hypothesis
#'     test
#' @param alternative string of either "one.sided" or "two.sided" for the
#'     alternative. If one.sided, the alternative means the test statistic is
#'     positive.
#' @param precBits precision of Rmpfri
#' @param fac A positive factor that makes the search over \eqn{v^T\mu} to be
#'     over a range of \code{(-fac*(max(y)-min(y)), fac*(max(y)- min(y)))}.
#'     Defaults to \code{10}.
#'
#' @return a vector of two numbers, the lower and upper end of the confidence interval
#' @export
confidence_interval <- function(y, polyhedra, contrast, sigma = 1, alpha = 0.05,
                                gridsize = 250, alternative = c("two.sided", "one.sided"), precBits = NA,
                                fac=10){
    alternative <- match.arg(alternative, c("two.sided", "one.sided"))
    coverage = 1-alpha

    diff <- max(y) - min(y)
    seq.val <- seq(-fac*diff, fac*diff, length.out = gridsize)
    myfun <- get_nonzero_mean_pv_fun(y,polyhedra$gamma, polyhedra$u, contrast,sigma,alpha)
    pvalue <- sapply(seq.val, myfun)

    if(alternative == "two.sided"){

        idx <- c(.select_index(pvalue, alpha/2, T),
                 .select_index(pvalue, 1-alpha/2, F))
        c(seq.val[idx[1]], seq.val[idx[2]])

    } else if (alternative == "one.sided"){

        idx <- .select_index(pvalue, alpha, T)
        c(seq.val[idx], Inf)

    } else {
        stop("|alternative| option not known.")
    }
}


##' Selects the index of \code{vec} that is larger than \code{alpha}. Or
##' something like that.
.select_index <- function(vec, alpha, lower = T){
    idx <- ifelse(lower, min(which(vec >= alpha)), max(which(vec <= alpha)))
    if(length(idx) == 0 | is.na(idx) | is.infinite(idx)){
        warning("numeric precision suspected to be too low")
        if(lower) return(1) else return(length(vec))
    }

    if(lower & vec[idx] > alpha & idx > 1) idx <- idx - 1
    if(!lower & vec[idx] < alpha & idx < length(vec)) idx <- idx + 1

    idx
}



##' Produces a function whose single input is the null value of \eqn{v^T\mu},
##' and whose output is the TG p-value according to that nonzero null.
get_nonzero_mean_pv_fun <- function(y, G, u, v, sigma, alpha,
                                    gridrange=c(-100,100),gridpts=100, bits=NULL) {

    ## Added by justin; not really needed but making usre
    ## if(sum(v*v)!=1) stop("provide a standardized contrast vector")

    z = sum(v*y)
    vv = sum(v^2)
    sd = sigma*sqrt(vv)

    rho = G %*% v / vv
    vec = (u - G %*% y + rho*z) / rho
    vlo = suppressWarnings(max(vec[rho>0]))
    vup = suppressWarnings(min(vec[rho<0]))

    xg = seq(gridrange[1]*sd,gridrange[2]*sd,length=gridpts)
    fun = function(x) { tnorm.surv(z,x,sd,vlo,vup,bits) }

    return(fun)
}



##' A more barebones, readable version of a one-sided CI from a
##' right-sided-alternative TG test. Matches the output of
##' \code{confidence_interval()} when using the \code{alternative="one-sided"}
##' option.
my.one.sided.ci <- function(y,poly,contrast,sigma,alpha,fac=10,
                            gridsize=250, griddepth=2){
    ## Basic check
    if(sum(contrast*y)<0) warning("Is your contrast designed for a right-direction one-sided
test? I doubt it, since you currently have a negative sum(contrast*y)..")

    ## Make a range of v^T mu to scan over.
    diff <- max(y) - min(y)
    seq.val <- seq(-fac*diff, fac*diff, length.out = gridsize)

    ## For each value of hypothetical values of v^TY, calculate one-sided test
    ## result
    pvs = pvalue(y, p, contrast, sigma = sigma, null_mean = seq.val,
                 alternative = c("one.sided"), precBits = NA)

    ## Reject
    reject = (pvs < alpha)
    cutoff = seq.val[min(which(!reject))]

    ## Calculate the one sided confidence interval
    return(c(cutoff, +Inf))
}
