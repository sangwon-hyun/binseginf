## functions to help convert a tree model into binSeg comparisons

.comparison_mat <- function(winning, losing){
  if(length(winning) != 3 | !is.numeric(winning)) stop("winning must be numeric of length 3")
  if(!is.matrix(losing) | !is.numeric(losing) | ncol(losing) != 3)
    stop("losing must be a numeric matrix with 3 columns")
  
  if(any(losing[,1] > losing[,2]) | any(losing[,2] >= losing[,3]))
    stop("losing must have second column >= first column and < third column")
  
  list(winning = winning, losing = losing)
}


.form_comparison <- function(tree, nodeName, breakpoint){
  winning <- matrix(NA, ncol = 3, nrow = 1)
  nodeNumeric <- .get_startEnd(nodeName)
  winning[c(1,3)] <- nodeNumeric; winning[2] <- breakpoint

  leaves.mat <- .get_leaves_matrix_excluding(tree, nodeName)

  losing <- .threeColumnMatrix_from_nodeVec(nodeNumeric, breakpoint)
  if(!any(is.na(leaves.mat))) losing <- rbind(losing, do.call(rbind, 
      .threeColumnMatrix_from_nodeMatrix(leaves.mat)))
  
  .comparison_mat(winning, losing)
}

.threeColumnMatrix_from_nodeMatrix <- function(mat){
  plyr::alply(mat, 2, .threeColumnMatrix_from_nodeVec)
}

.threeColumnMatrix_from_nodeVec <- function(vec, exclude = NA){
  stopifnot(length(vec) == 2, is.numeric(vec))
  stopifnot(is.na(exclude) || (is.numeric(exclude) & length(exclude) == 1) & 
      exclude >= vec[1] & exclude < vec[2])
  
  if(is.na(exclude)){
    cbind(vec[1], vec[1]:(vec[2]-1), vec[2])
  } else{
    mid.vec <- vec[1]:(vec[2]-1)
    mid.vec <- mid.vec[mid.vec != exclude]
    if(length(mid.vec) == 0) return(numeric(0))
    cbind(vec[1], mid.vec, vec[2])
  }
}

.get_leaves_matrix_excluding <- function(tree, nodeName){
  leaves.names <- .get_leaves_names(tree)
  stopifnot(nodeName %in% leaves.names)
  leaves.names <- leaves.names[leaves.names != nodeName]
  
  #remove leaves if start = end
  idx <- sapply(leaves.names, function(x){
    tmp <- strsplit(x, split = "-")
    if(tmp[[1]][1] == tmp[[1]][2]) return(FALSE) else return(TRUE)
  })
  if(length(idx) > 0) leaves.names <- leaves.names[which(idx)]
  
  if(length(leaves.names) == 0){
    return(NA)
  } else {
    sapply(leaves.names, .get_startEnd)
  }
}

.get_startEnd <- function(nodeName){
  stopifnot(length(nodeName) == 1, is.character(nodeName), grep("-", nodeName) == 1)
  
  res <- as.numeric(strsplit(nodeName, split = "-")[[1]])
  stopifnot(length(res) == 2, is.numeric(res))
  
  res
}

.get_comparisonSigns <- function(y, mat){
  stopifnot(is.matrix(mat), ncol(mat) == 3)
  
  values <- apply(mat, 1, function(x){
    .cusum(y, x[1], x[2], x[3])
  })
  
  sign(values)
}