#' Construct a hyperplane
#'
#' Hyperplane is of form \code{a%*%x = b}.
#' If \code{a} is a matrix, this object represents an intersection of hyperplanes.
#'
#' @param a matrix
#' @param b vector
#'
#' @return a \code{plane} object
.plane <- function(a, b = 0){
  if(!is.matrix(a)){a <- matrix(a, nrow = 1)}
  stopifnot(length(b) == nrow(a))
  stopifnot(length(which(a != 0)) > 0)

  l2_vec <- apply(a, 1, .l2norm)
  b <- as.numeric(b/l2_vec)
  if(length(l2_vec) > 1){
    a <- diag(1/l2_vec)%*%a
  } else {
    a <- a/l2_vec
  }

  structure(list(a = a, b = b), class = "plane")
}
