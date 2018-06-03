#' Binary segmentation with fixed steps

#' y must not have duplicated values. This is to avoid
#' degenerate behavior of binary segmentation
#'
#' @param y numeric vector to contain data
#' @param numSteps numeric of number of steps
#'
#' @return a bsFs object
#' @export
binSeg_fixedSteps <- function(y, numSteps){
  if(numSteps >= length(y)) stop("numSteps must be strictly smaller than the length of y")
  if(numSteps <= 0) step("numSteps must be at least 1.")

  #initialization
  n <- length(y); tree <- .create_node(1, n)
  cp <- c()

  for(steps in 1:numSteps){
    leaves.names <- .get_leaves_names(tree)
    for(i in 1:length(leaves.names)){
      leaf <- data.tree::FindNode(tree, leaves.names[i])

      res <- .find_breakpoint(y, leaf$start, leaf$end)

      leaf$breakpoint <- res$breakpoint; leaf$cusum <- res$cusum
    }

    node.name <- .find_leadingBreakpoint(tree)
    node.selected <- data.tree::FindNode(tree, node.name)
    node.selected$active <- steps
    node.pairs <- .split_node(node.selected)
    node.selected$AddChildNode(node.pairs$left)
    node.selected$AddChildNode(node.pairs$right)
  }

  y.fit <- .refit_binseg(y, jumps(tree))
  obj <- structure(list(tree = tree, y.fit = y.fit, numSteps = numSteps), class = "bsFs")
  cp <- jumps(obj)
  leaves <- .enumerate_splits(tree)
  cp.sign <- sign(as.numeric(sapply(leaves, function(x){
      data.tree::FindNode(tree, x)$cusum})))
  obj <- structure(list(tree = tree, y.fit = y.fit, numSteps = numSteps, cp = cp,
                        cp.sign=cp.sign, y=y), class = "bsFs")

}

#' is_valid for bsFs
#'
#' @param obj bsFs object
#'
#' @return TRUE if valid
#' @export
is_valid.bsFs <- function(obj){
  if(class(obj$tree)[1] != "Node") stop("obj$tree must a Node")
  if(!is.numeric(obj$numSteps)) stop("obj$numSteps must be a numeric")
  if(length(.enumerate_splits(obj$tree)) != obj$numSteps)
    stop("obj$tree and obj$numSteps disagree")

  TRUE
}

#' Get jumps from bsFs objects
#'
#' Enumerates the jumps. Sorted = F will return the jumps in order
#' of occurance in the binSeg algorithm. Sorted = T will list the jumps
#' in numeric order
#'
#' @param obj bsFs object
#' @param sorted boolean
#' @param ... not used
#'
#' @return vector of jumps
#' @export
jumps.bsFs <- function(obj, sorted = F, ...){
  jumps(obj$tree, sorted)
}

#' Get the cusum for jumps for bsFs objects
#'
#' Enumerates the cusum for each jump. Sorted = F will return the jumps in order
#' of occurance in the binSeg algorithm. Sorted = T will list the jumps
#' in numeric order
#'
#' @param obj  bsFs object
#' @param sorted  boolean
#' @param ... not use
#'
#' @return vector of cusum numerics
#' @export
jump_cusum.bsFs <- function(obj, sorted = F, ...){
  jump_cusum(obj$tree, sorted)
}

##' Summary of bsFs object
##'
##' @param object  bsFs object
##' @param ... not used
##'
##' @return matrix of summary statistics
##' @export
summary.bsFs <- function(object, ...){
  summary(object$tree)
}

.refit_binseg <- function(y, jumps){
  stopifnot(max(jumps) < length(y), min(jumps) > 0)
  stopifnot(all(jumps %% 1 == 0), length(jumps) == length(unique(jumps)))

  n <- length(y); y.fit <- numeric(n)
  jumps <- c(0, sort(jumps), n)
  for(i in 2:length(jumps)){
    y.fit[(jumps[i-1]+1):jumps[i]] <- mean(y[(jumps[i-1]+1):jumps[i]])
  }

  y.fit
}

.find_breakpoint <- function(y, start, end){
  ## stopifnot(!any(duplicated(y)))
  if(start > end) stop("start must be smaller than or equal to end")
  if(start == end) return(list(breakpoint = start, cusum = 0))

  breakpoint <- seq(from = start, to = end - 1, by = 1)
  cusum.vec <- sapply(breakpoint, .cusum, y = y, start = start, end = end)

  idx <- which.max(abs(cusum.vec))
  list(breakpoint = breakpoint[idx], cusum = cusum.vec[idx])
}

.cusum <- function(y, start, idx, end){
  v <- .cusum_contrast(start, idx, end)
  as.numeric(v %*% y[start:end])
}

#n1 is denoted as start to idx (inclusive)
.cusum_contrast <- function(start, idx, end){
  if(start > idx) stop("start must be smaller or equal than idx")
  if(idx >= end) stop("idx must be smaller to end")

  n1 <- idx - start + 1
  n2 <- end - idx

  c(rep(-1/n1, n1), rep(1/n2, n2)) * sqrt(1/((1/n1) + (1/n2)))
}

.cusum_contrast_full <- function(start, idx, end, n){
  res <- rep(0, n)
  res[start:end] <- .cusum_contrast(start, idx, end)

  res
}


##' Print function for convenience, of |wbs| class object.
print.bsFs <- function(obj){
    cat("Detected changepoints using WBS with", obj$numSteps, "steps is", obj$cp * obj$cp.sign, fill=TRUE)
}
