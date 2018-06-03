#' Generate polyhedra matrix from flFs object
#'
#' Forms both Gamma matrix and u vector
#' 
#' @param obj flFs object
#' @param ... not used
#'
#' @return An object of class polyhedra
#' @export
polyhedra.flFs <- function(obj, ...){
  k <- nrow(obj$model); n <- length(obj$y.fit)
  gamma.row.lis <- vector("list", k)
  D <- .form_Dmatrix(n)

  for(i in 1:k){
    gamma.row.lis[[i]] <- .gammaRows_from_flasso(n, D, obj$model[1:i,])
  }

  mat <- do.call(rbind, gamma.row.lis)
  polyhedra(obj = mat, u = rep(0, nrow(mat)))
}

.gammaRows_from_flasso <- function(n, D, model){
  sign.win <- model$Sign[nrow(model)]
  if(nrow(model) == 1){ 
    idx <- 1:(n-1) 
    denominator.vec <- .compute_fused_denominator(D, idx, model[-1,])
  } else {
    idx <- .select_nonactive(n, model$Index[1:(nrow(model) - 1)])
    denominator.vec <- .compute_fused_denominator(D, idx, model[1:(nrow(model) - 1),])
  }

  numerator.mat <- .compute_fused_numerator_polyhedra(D, idx)
  active.idx <- which(idx == model$Index[nrow(model)])

  contrasts <- .form_contrast_flasso(numerator.mat, denominator.vec, sign.win, active.idx)

  .vector_matrix_signedDiff(contrasts$win, contrasts$lose, 1, 
    rep(1, nrow(contrasts$lose)))
}

.compute_fused_numerator_polyhedra <- function(D, idx){
  stopifnot(all(idx %% 1 == 0), !any(duplicated(idx)))
  stopifnot(min(idx) >= 1, max(idx) <= nrow(D))
  
  DDT <- Matrix::tcrossprod(D[idx,,drop = F])
  Didx <- D[idx,,drop = F]
  
  Matrix::solve(DDT)%*%Didx
}

.form_contrast_flasso <- function(numerator.mat, denominator.vec,
  sign.win, active.idx, tol = 1e-7){
  stopifnot(active.idx <= nrow(numerator.mat))
  
  win <- numerator.mat[active.idx,]/(sign.win + denominator.vec[active.idx])

  denom <- denominator.vec[-active.idx]
  pos.denom <- denom; pos.denom[abs(denom + 1) < tol] <- 0
  neg.denom <- denom; neg.denom[abs(denom - 1) < tol] <- 0

  lose1 <- numerator.mat[-active.idx,]/(1 + pos.denom)
  lose2 <- numerator.mat[-active.idx,]/(-1 + neg.denom)

  list(win = win, lose = rbind(as.matrix(lose1), as.matrix(lose2)))
}