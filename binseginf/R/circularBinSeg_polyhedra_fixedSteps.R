#' Polyhedra for circular binary segmentation fixed steps
#'
#' The polyhedra has elements \code{gamma} (a matrix) and \code{u}
#' (a vector).
#'
#' @param obj cbsFs object
#' @param ... void, not used
#'
#' @return a polyhedra
#' @export
polyhedra.cbsFs <- function(obj, ...){
  n <- .get_startEnd(obj$tree$name)[2]
  numSteps <- obj$numSteps
  comp_lis <- .list_comparison(obj)
  sign_vec <- sign(jump_cusum(obj))
  gamma_row_lis <- vector("list", numSteps)

  for(i in 1:numSteps){
    gamma_row_lis[[i]] <- .gammaRows_from_comparisons_cbsfs(comp_lis[[i]]$winning,
                                                            comp_lis[[i]]$losing, sign_vec[i], n)
  }

  mat <- do.call(rbind, gamma_row_lis)
  polyhedra(obj = mat, u = rep(0, nrow(mat)))
}

#' Form contrast vector for circular binary segmentation
#'
#' The vector is made such that for a data vector \code{y}, if the
#' output of \code{.cusum_cbs_contrast_full} is \code{res}, then
#' \code{res \%*\% y} is the circular cusum statistic for a particular
#' set of indices.
#'
#' Here, \code{start}, \code{idx}, \code{end} are all integers between 1
#' and \code{n}
#'
#' @param start integer denoting the first (left) index of the left shoulder
#' @param idx a vector of 2 integers containing the first (left) index
#' of the hump and the last (right) index of the hump
#' @param end integer denoting the last (right) index of the right shoulder
#' @param n length of the contrast vector
#'
#' @return a vector
.cusum_cbs_contrast_full <- function(start, idx, end, n){
  stopifnot(length(idx) == 2, idx[1] >= 1, idx[2] <= n, idx[1] <= idx[2], start >= 1, end <= n)
  stopifnot(all(idx %% 1 == 0), start %% 1 == 0, end %% 1 == 0)
  stopifnot(!all(idx[1]==1, idx[2]==n))

  m <- idx[2]-idx[1]+1; n2 <- end-start+1
  const <- sqrt(1/(1/m + 1/(n2-m)))

  vec <- rep(0, n)
  if(idx[1] > start) vec[start:(idx[1]-1)] <- -1/(n2-m)
  vec[idx[1]:idx[2]] <- 1/m
  if(idx[2] < end) vec[(idx[2]+1):end] <- -1/(n2-m)

  const * vec
}

#' Form the lines of the gamma matrix
#'
#' \code{vec} is a vector of 4 positive integers less than \code{n}, where
#' (from left to right) denote the first (left) index of the left shoulder,
#' the first (left) index of the hump, the last (right) index of the hump,
#' and the last (right) index of the right shoulder. Similarly,
#' \code{mat} is a matrix with 4 columns.
#'
#' Typically, this function is used such that the circular cusum
#' statistic formed with the indices in \code{vec} is larger than the
#' circular cusum statistic formed with the indicies in any of the rows of
#' \code{mat}.
#'
#' This function returns a matrix with \code{n} columns, and \code{2*nrow(mat)}
#' rows.
#'
#' @param vec a vector of integers
#' @param mat a matrix of integers
#' @param sign_win 1 or -1, denoting the sign of the jump denoted in \code{vec}
#' @param n number of elements
#' @param add boolean
#'
#' @return a matrix
.gammaRows_from_comparisons_cbsfs <- function(vec, mat, sign_win, n, add = F){
  stopifnot(length(vec) == 4, ncol(mat) == 4)
  stopifnot(add | !any(is.na(mat)))

  if(add & any(is.na(mat))){
    win_contrast <- .cusum_cbs_contrast_full(vec[1], vec[2:3], vec[4], n)
    return(matrix(sign_win*win_contrast, 1, n))
  }

  lose_contrast <- t(apply(mat, 1, function(x){
    .cusum_cbs_contrast_full(x[1], x[2:3], x[4], n)
  }))

  if(!any(is.na(vec))){
    win_contrast <- .cusum_cbs_contrast_full(vec[1], vec[2:3], vec[4], n)

    # add inequalities to compare winning split to all other splits
    res <- .vector_matrix_signedDiff(win_contrast, lose_contrast, sign_win,
                                     rep(1, nrow(lose_contrast)))
    res2 <- .vector_matrix_signedDiff(win_contrast, lose_contrast, sign_win,
                                      -rep(1, nrow(lose_contrast)))

    res_mat <- rbind(res, res2)

    if(add) rbind(res_mat, sign_win*win_contrast) else res_mat
  } else {
    rbind(lose_contrast, -lose_contrast)
  }
}
