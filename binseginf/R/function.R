#' L2 norm
#'
#' @param x vector
#'
#' @return scalar
#' @export
l2norm <- function(x){
  sqrt(sum(x^2))
}