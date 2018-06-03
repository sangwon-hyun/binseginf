#' Generate contrast vector
#'
#' @param obj object
#' @param jump.idx index among the list of jumps to jump at
#' @param sorted boolean on whether or not jumps should be sorted
#' @param type segment or spike
#' @param ... not used
#'
#' @return numeric vector
#' @export
contrast_vector <- function(obj, jump.idx, sorted = F,
  type = c("segment", "spike"), ...){
  stopifnot(class(obj) %in% c("bsFs", "flFs","wbsFs"))
  if(!is.character(type[1])) stop("type must be character")
  if(!type[1] %in% c("segment", "spike")) stop("type must be either segment or spike")
  type <- type[1]

  jump.vec <- jumps(obj, sorted)
  jump <- jump.vec[jump.idx]
  n <- .get_length(obj)

  jump_mat <- jump_sign(obj)
  jump_sign <- jump_mat[which(jump_mat[,"Jump"] == jump), "Sign"]

  if(type == "segment") {
    v <- .contrast_vector_segment(obj, jump, n)
  } else {
    v <- .contrast_vector_spike(obj, jump, n)
  }

  v * jump_sign
}

.get_length.bsFs <- function(obj){
  .get_startEnd(obj$tree$name)[2]
}

.get_length.flFs <- function(obj){
  length(obj$y.fit)
}

.contrast_vector_segment <- function(obj, jump, n){
  jumpSorted.vec <- c(1, jumps(obj, T), n)

  idx <- which(jumpSorted.vec == jump)[1]
  if(idx == 1) idx <- 2

  start <- jumpSorted.vec[idx-1]; split <- jump; end <- jumpSorted.vec[idx+1]
  res <- sign(.cusum_contrast_full(start, split, end, n))
  res[res>0] <- 1/sum(res > 0)
  res[res<0] <- -1/sum(res < 0)

  res
}

.contrast_vector_spike <- function(obj, jump, n){
  res <- rep(0, n)
  res[c(jump, jump+1)] <- c(-1, 1)

  res
}
