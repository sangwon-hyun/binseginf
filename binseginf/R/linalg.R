#' Projection of vector onto another vector
#'
#' Returns the component of \code{vec1} that is orthogonal to \code{vec2}
#'
#' @param vec1 vector
#' @param vec2 vector
#' @param tol small positive number
#'
#' @return vector
.projection <- function(vec1, vec2, tol = 1e-6){
  stopifnot(length(vec1) == length(vec2))

  d <- length(vec1)
  vec2 <- vec2/.l2norm(vec2)
  as.numeric((diag(d) - vec2%*%t(vec2))%*%vec1)
}

#' Projection of vector onto rows of a matrix
#'
#' Returns the component of \code{vec} that is orthogonal to all the rows
#' of \code{mat}
#'
#' @param vec vector
#' @param mat matrix
#'
#' @return vector
.projection_matrix <- function(vec, mat){
  stopifnot(length(vec) == ncol(mat), ncol(mat) >= nrow(mat))
  n <- length(vec)
  proj_mat <- diag(n) - t(mat) %*% solve(mat %*% t(mat)) %*% mat

  as.numeric(proj_mat %*% vec)
}

#' Sample unit vectors from the null space of a matrix
#'
#' This function assumes the row space of \code{mat} is the number of
#' rows in \code{mat}. If \code{num_vec} is \code{NA}, then the function
#' construct a basis for the entire null space of \code{mat}, which is
#' assumed to have dimension equal to the rank of \code{mat}.
#'
#' If \code{null} is \code{TRUE}, then this function samples from the nullspace.
#' Otherwise, sample from the rowspace.
#'
#' @param mat matrix
#' @param num_vec positive integer or \code{NA}
#' @param null boolean
#'
#' @return a matrix of vectors, with number of rows equal to \code{ncol(mat)}
#' and number of columns equal to \code{num_vec}
.sample_matrix_space <- function(mat, num_vec = NA, null = T){
  k <- Matrix::rankMatrix(mat)
  n <- ncol(mat)

  # form basis of either rowspace or nullspace
  if(null){
    res <- matrix(stats::rnorm(n*(n-k)), ncol = n-k)
    res <- apply(res, 2, .projection_matrix, mat = mat)
    res <- svd(t(res))$v
  } else {
    res <- svd(mat)$v
  }

  if(is.na(num_vec)) return(res)

  weights <- matrix(stats::rnorm(num_vec*ncol(res)), nrow = ncol(res))
  vec_mat <- res %*% weights

  #make the vectors orthogonal
  if(num_vec > 1){
    for(i in 2:num_vec){
      vec_mat[,i] <- .projection_matrix(vec_mat[,i], mat = t(vec_mat[,1:(i-1), drop = F]))
    }
  }

  sapply(1:ncol(vec_mat), function(x){vec_mat[,x]/.l2norm(vec_mat[,x])})
}

.compute_nullspace <- function(mat){
  k <- Matrix::rankMatrix(mat)
  n <- ncol(mat)

  res <- matrix(stats::rnorm(n*(n-k)), ncol = n-k)
  res <- apply(res, 2, .projection_matrix, mat = mat)
  svd(t(res))$v
}

.sample_two_vectors <- function(mat){
  num_vec <- 2
  weights <- matrix(stats::rnorm(num_vec*ncol(mat)), nrow = ncol(mat))
  vec_mat <- mat %*% weights

  vec_mat[,2] <- .projection(vec_mat[,2], vec_mat[,1])

  sapply(1:ncol(vec_mat), function(x){vec_mat[,x]/.l2norm(vec_mat[,x])})
}
