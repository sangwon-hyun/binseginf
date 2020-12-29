#' Form comparisons for circular binary segmentation
#' 
#' For a particular node in \code{nodeName} (for example, \code{"1-10"}) and
#' a desired \code{breakpoint} (for example \code{c(3,5)}), form matrices
#' with 4 columns on all the possible comparisons made.
#' 
#' The returned object has a \code{winning} element, a matrix with 1 row and
#' 4 columns. It denotes where the split dictated by \code{breakpoint} occurs,
#' in this example, \code{matrix(c(1,3,5,10), nrow = 1)}. 
#' 
#' The returned object also has a \code{losing} element, a matrix with 4 columns.
#' Each row denotes all the other splits that could have occurred that is not
#' in \code{winning}.
#' 
#' \code{excluding} is a boolean. If \code{TRUE}, then the comparisons will
#' include other nodes in \code{losing}. Otherwise, it will not. Also, \code{breakpoint}
#' can be \code{NA}.
#'
#' @param tree tree object
#' @param nodeName node name to split on right now
#' @param breakpoint vector of 2 positive integers where the split will occur
#' @param excluding boolean
#'
#' @return a list
.form_comparison_cbs <- function(tree, nodeName, breakpoint, excluding = T){
  stopifnot(length(breakpoint) == 2 | is.na(breakpoint))
  
  winning <- matrix(NA, ncol = 4, nrow = 1)
  nodeNumeric <- .get_startEnd(nodeName)
  winning[c(1,4)] <- nodeNumeric; winning[2:3] <- breakpoint
  
  if(excluding) leaves.mat <- .get_leaves_matrix_excluding(tree, nodeName)
  
  losing <- .fourColumnMatrix_from_nodeVec(nodeNumeric, breakpoint)
  if(excluding && !any(is.na(leaves.mat))) {
    losing <- rbind(losing, do.call(rbind, .fourColumnMatrix_from_nodeMatrix(leaves.mat)))
  }
  
  list(winning = winning, losing = losing)
}

#' Form comparisons from a matrix
#' 
#' Form the comparisons in a 4-column matrix 
#'
#' @param mat a matrix with 2 columns, where each row denotes the left (first)
#' index of the left shoulder of some segment, and the right (last) index of the
#' right shoulder of the same segment.
#'
#' @return a matrix
.fourColumnMatrix_from_nodeMatrix <- function(mat){
  plyr::alply(mat, 2, .fourColumnMatrix_from_nodeVec)
}

#' Form comparisons from a vector
#' 
#' This function returns a 4-column matrix of positive integers.
#' 
#' \code{vec} is two positive integers in increasing order. This
#' returns all the possible circular binary segmentation comparisons
#' where \code{vec[1]} is the left (first) index of the left shoulder, and \code{vec[2]} 
#' is the right (last) index of the right shoulder.
#' 
#' \code{exclude} is an optional argument. If included (i.e., not \code{NA}),
#' it should be a vector of 2 positive integers and
#' the returned matrix will omit the row including \code{exclude}.
#'
#' @param vec vector of 2 integers
#' @param exclude \code{NA} or a vector of 2 integers
#'
#' @return a matrix
.fourColumnMatrix_from_nodeVec <- function(vec, exclude = NA){
  stopifnot(length(vec) == 2, is.numeric(vec), vec[2] > vec[1])
  stopifnot(any(is.na(exclude)) || (is.numeric(exclude) & length(exclude) == 2) & 
              exclude[1] >= vec[1] & exclude[2] <= vec[2])
  
  if(any(is.na(exclude))){
    cbind(vec[1], .enumerate_breakpoints_cbs(vec[2]-vec[1]+1, vec[1]), vec[2])
  } else if(vec[2]-vec[1] == 1) {
    matrix(NA, 1, 4)
  } else {
    mid.vec <- .enumerate_breakpoints_cbs(vec[2]-vec[1]+1, vec[1])
    mid.vec <- mid.vec[-which(apply(cbind(exclude[1] == mid.vec[,1], 
                                         exclude[2] == mid.vec[,2]), 1, all)),]
    if(length(mid.vec) == 0) return(numeric(0))
    cbind(vec[1], mid.vec, vec[2])
  }
}