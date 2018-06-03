#' CpVector Class Constructor
#' 
#' Creates an instance of the CpVector, short for "Change Point Vector".
#' A CpVector objects has 3 elements, data (the realization from the
#' changepoint model), jump.height (the height of jumps in the true signal)
#' and jump.loc (the locations of the jumps in the true signal).
#' 
#' Here, CpVector is draws samples from a piecewise constant true signal from
#' points 1 to n. The jump.loc is a vector between 0 and 1 (inclusive and 
#' exclusive repsectively), and CpVector automatically scales jump.loc to
#' fit between 1 and n. Here, jump.loc must always be one less element than
#' jump.height.
#' 
#' Also, n must be large enough with respect to jump.loc that round(n*jump.loc)
#' must be distinct integers. Otherwise, an error will be thrown.
#' 
#' Here, func is a function that dictates the noise model. For example, the
#' default is rnorm.
#'
#' @param n the number of realizations
#' @param jump.height a numeric vector of the jump heights
#' @param jump.loc a numeric vector (between 0 (inclusive) and 1 (exclusive))
#' of the jump locations between 0 and 1.
#' @param func the function to dictate the noise model. The first parameter
#' of func must dictate how many samples are drawn.
#' @param ... additional parameters to pass into func.
#'
#' @return a CpVector instance
#' @export
CpVector <- function(n, jump.height = 0, jump.loc = NA, func = stats::rnorm, 
  ...){
  
  if(!is.numeric(n) | !is.numeric(jump.height))
    stop("n and jump.height must be numeric")
  if(!any(is.na(jump.loc))){
    if(!is.numeric(jump.loc)) stop("jump.loc must be numeric")
    if(length(jump.height) != length(jump.loc) + 1) 
      stop("jump.height must be one element more than jump.loc")
    if(min(jump.loc) < 0 | max(jump.loc) >= 1)
      stop("jump.loc must lie between 0 (inclusive) and 1 (exclusive)")
    if(!all(diff(jump.loc) > 0)) stop("jump.loc must be strictly increasing")
  } else {
    if(length(jump.height) != 1)
      stop("jump.height must be one element since jump.loc = NA")
    if(length(jump.loc) != 1)
      stop("jump.loc must be one element if it contains any NA")
  }
  
  if(n %% 1 != 0 | n < 0) stop("n must be a positive integer")
  
  jump.idx <- .computeJumpIdx(n, jump.loc)
  mean.vec <- .formMeanVec(n, jump.height, jump.idx)
  
  data <- mean.vec + func(n, ...)
  
  structure(list(data = data, jump.height = jump.height, jump.idx = jump.idx),
    class = "CpVector")
}

.computeJumpIdx <- function(n, jump.loc){
  if(any(is.na(jump.loc))) return(numeric(0))
  
  vec <- round(n*jump.loc)
  vec <- sapply(vec, function(x) {max(min(x,n),1)})
  
  if(length(vec) != length(unique(vec))) stop(paste("n is too small compared",
    "to jump.loc that changepoints are not unique"))
  
  vec
}

.formMeanVec <- function(n, jump.height, jump.idx){
  diff.vec <- diff(c(0,jump.idx,n))
  rep(jump.height, times = diff.vec)
}
