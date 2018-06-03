#' Fused lasso with fixed steps
#'
#' y must not have duplicated values. This is to avoid
#' degenerate behavior of binary segmentation. The reported \code{y.fit}
#' is set the lambda where the (k+1)th jump is about to appear.
#'
#' @param y  numeric vector to contain data
#' @param numSteps  numeric of number of steps
#' @param tol tolerance to handle divide by zero instances
#'
#' @return a flFs object
#' @export
fLasso_fixedSteps <- function(y, numSteps, tol = 1e-7){
  if(numSteps >= length(y)) stop("numSteps must be strictly smaller than the length of y")
  if(any(duplicated(y))) stop("y must contain all unique values")

  #initialization
  n <- length(y)
  model.mat <- as.data.frame(matrix(NA, numSteps + 1, 3))
  colnames(model.mat) <- c("Index", "Sign", "Lambda")
  D <- .form_Dmatrix(n)

  for(steps in 1:(numSteps + 1)){
    idx <- .select_nonactive(n, model.mat$Index)
    a.vec <- .compute_fused_numerator(D, idx, y)
    b.vec <- .compute_fused_denominator(D, idx, model.mat[1:(steps-1),,drop = F])
    
    pos.ratio <- a.vec/(1+b.vec)
    pos.ratio[abs(1 + b.vec) < tol] <- 0
    neg.ratio <- a.vec/(-1+b.vec)
    neg.ratio[abs(-1 + b.vec) < tol] <- 0 

    if(max(pos.ratio) > max(neg.ratio)){
      model.mat[steps,] <- c(idx[which.max(pos.ratio)], 1, max(pos.ratio))
    } else {
      model.mat[steps,] <- c(idx[which.max(neg.ratio)], -1, max(neg.ratio))
    }
  }
  
  y.fit <- .refit_flasso(y, model.mat)

  structure(list(model = model.mat[1:numSteps, ], y.fit = y.fit, 
    numSteps = numSteps), class = "flFs")
}

#' Get jumps from flFs objects
#' 
#' Enumerates the jumps. Sorted = F will return the jumps in order
#' of occurance in the binSeg algorithm. Sorted = T will list the jumps
#' in numeric order
#'
#' @param obj flFs object
#' @param sorted boolean
#' @param ... not used
#'
#' @return vector of jumps
#' @export
jumps.flFs <- function(obj, sorted = F, ...){
  idx <- obj$model$Index
  if(sorted) sort(idx) else idx
}

#' Get lambdas for flFs objects
#'
#' Enumerates the decreasing sequence of lambdas for the jumps.
#' Sorted = F will return the jumps in order
#' of occurance in the binSeg algorithm. Sorted = T will list the jumps
#' in numeric order
#'
#' @param obj flFs object
#' @param ... not used
#'
#' @return vector of decreasing lambdas
#' @export
jump_lambda.flFs <- function(obj, ...){
  obj$model$Lambda
}

#' Summary for flFs object
#'
#' @param object flFs object
#' @param ... not used
#'
#' @return summary of flFs (data frame)
#' @export
summary.flFs <- function(object, ...){
  object$model
}

.refit_flasso <- function(y, model.mat){
  n <- length(y)
  D <- .form_Dmatrix(n)
  u.vec <- rep(NA, n-1)
  lambda <- min(model.mat$Lambda)
  u.vec[model.mat$Index] <- lambda*model.mat$Sign
  
  idx <- .select_nonactive(n, model.mat$Index)
  a.vec <- .compute_fused_numerator(D, idx, y)
  b.vec <- .compute_fused_denominator(D, idx, model.mat)
  
  nonactive.idx <- c(1:(n-1))[-model.mat$Index]
  u.vec[nonactive.idx] <- a.vec - lambda*b.vec
  
  as.numeric(y - Matrix::crossprod(.form_Dmatrix(n)[-n,], u.vec))
}

#is square
.form_Dmatrix <- function(n){
  i <- rep(1:(n-1), each = 2)
  j <- c(1, rep(2:(n-1), each = 2), n)
  x <- rep(c(-1, 1), times = n-1)
  
  Matrix::sparseMatrix(i, j, x = x, triangular= TRUE)
}

.select_nonactive <- function(n, vec){
  val <- vec[!is.na(vec)]
  if(length(val) == 0) 1:(n-1) else c(1:(n-1))[-val]
}

.compute_fused_numerator <- function(D, idx, y){
  as.numeric(.compute_fused_numerator_polyhedra(D, idx) %*% y)
}

.compute_fused_denominator <- function(D, idx, model.mat){
  stopifnot(all(idx %% 1 == 0), !any(duplicated(idx)))
  stopifnot(min(idx) >= 1, max(idx) <= ncol(D) - 1)
  
  if(length(idx) == nrow(D) || length(model.mat) == 0 || any(is.na(model.mat$Index))) 
    return(rep(0, nrow(D) - 1))
  
  active.idx <- model.mat$Index; sign.vec <- model.mat$Sign
  DDT <- Matrix::tcrossprod(D[idx,,drop = F])
  DDTs <- Matrix::tcrossprod(D[idx,,drop = F], D[active.idx,,drop = F]) %*% sign.vec
  
  as.numeric(Matrix::solve(DDT, DDTs))
}