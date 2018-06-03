.create_node <- function(start, end, breakpoint = NA, cusum = NA, active = NA){
  node <- data.tree::Node$new(paste0(start, "-", end))

  node$start <- start
  node$end <- end
  node$breakpoint <- breakpoint
  node$cusum <- cusum
  node$active <- active

  is_valid(node)

  node
}

#' Check whether tree is valid or not
#'
#' @param obj The tree (of class Node)
#'
#' @return TRUE if valid
#' @export
is_valid.Node <- function(obj){
  if(obj$start > obj$end) stop("the start must be less or equal to end")
  if(!is.na(obj$breakpoint) & (obj$start > obj$breakpoint & obj$end < obj$breakpoint))
    stop("breakpoint must be between start and end (inclusive)")

  TRUE
}

#' Get jumps from a tree
#'
#' Enumerates the jumps. Sorted = F will return the jumps in order
#' of occurance in the binSeg algorithm. Sorted = T will list the jumps
#' in numeric order
#'
#' @param obj object of class Node
#' @param sorted boolean
#' @param ... not used
#'
#' @return vector of jumps
#' @export
jumps.Node <- function(obj, sorted = F, ...){
  leaves <- .enumerate_splits(obj)
  res <- as.numeric(sapply(leaves, function(x){data.tree::FindNode(obj, x)$breakpoint}))
  ## res <- sapply(leaves, function(x){data.tree::FindNode(obj, x)$breakpoint}) ## delete later./>>>

  if(sorted) sort(res) else res
}

#' Get the cusum for jumps for Node objects
#'
#' Enumerates the cusum for each jump. Sorted = F will return the jumps in order
#' of occurance in the binSeg algorithm. Sorted = T will list the jumps
#' in numeric order
#'
#' @param obj  Node object
#' @param sorted  boolean
#' @param ... not used
#'
#' @return vector of cusum numerics
#' @export
jump_cusum.Node <- function(obj, sorted = F, ...){
  leaves <- .enumerate_splits(obj)

  res <- sapply(leaves, function(x){data.tree::FindNode(obj, x)$cusum})

  if(sorted){
    jumps <- jumps(obj, sorted = T)
    idx <- order(jumps)
    res[idx]
  } else {
    res
  }
}

#' Summary of Node object
#'
#' @param object Node object
#' @param ... not used
#'
#' @return matrix of summary statistics
#' @export
summary.Node <- function(object, ...){
  leaves <- .enumerate_splits(object)
  jumps <- jumps(object)
  cusum <- jump_cusum(object)

  mat <- t(sapply(leaves, .get_startEnd))
  dat <- data.frame("Split_Number" = 1:length(leaves), "Start" = mat[,1],
    "End" = mat[,2], "Breakpoint" = jumps, "Cusum" = cusum)
  rownames(dat) <- NULL

  dat
}

.get_leaves_names <- function(tree){
  leaves <- tree$leaves
  vec <- sapply(leaves, function(x){x$name})
  names(vec) <- NULL
  sort(vec)
}

.find_leadingBreakpoint <- function(tree){
  leaves.names <- .get_leaves_names(tree)
  cusum.vec <- sapply(leaves.names, function(x){
    data.tree::FindNode(tree, x)$cusum
  })

  leaves.names[which.max(abs(cusum.vec))]
}

.split_node <- function(node){
  if(is.na(node$breakpoint)) stop("node does not have a set breakpoint yet")
  if(node$breakpoint >= node$end) stop("node breakpoint must be less than end")

  left <- .create_node(node$start, node$breakpoint)
  right <- .create_node(node$breakpoint + 1, node$end)

  list(left = left, right = right)
}

.enumerate_splits <- function(tree){
  names(sort(tree$Get("active")))
}

