#' Computes the Hausdorff Distance between TWo Sets
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

