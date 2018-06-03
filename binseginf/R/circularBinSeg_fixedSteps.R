#' Circular binary segmentation - fixed steps
#'
#' Returns a list containing a tree denoting how the input \code{y} is
#' segmented and the input \code{numSteps}. If the output is stored
#' in a variable \code{res}, use \code{jumps(res)} to enumerate all the
#' jump locations.
#'
#' @param y data, numeric vector
#' @param numSteps integer for number of jumps
#'
#' @return cbs object
#' @export
circularBinSeg_fixedSteps <- function(y, numSteps){
  if(numSteps >= length(y)) stop("numSteps must be strictly smaller than the length of y")

  #initialization
  n <- length(y); tree <- .create_node(1, n); vec <- cumsum(y)

  for(steps in 1:numSteps){
    leaves.names <- .get_leaves_names(tree)
    for(i in 1:length(leaves.names)){
      leaf <- data.tree::FindNode(tree, leaves.names[i])

      res <- .find_breakpoint_cbs(y[leaf$start:leaf$end])

      leaf$breakpoint <- res$breakpoint+leaf$start-1; leaf$cusum <- res$cusum
    }

    node_name <- .find_leadingBreakpoint(tree)
    node_selected <- data.tree::FindNode(tree, node_name)
    node_selected$active <- steps
    node_pairs <- .split_node_cbs(node_selected)
    if(!any(is.na(node_pairs$left))) node_selected$AddChildNode(node_pairs$left)
    node_selected$AddChildNode(node_pairs$middle)
    if(!any(is.na(node_pairs$right))) node_selected$AddChildNode(node_pairs$right)
  }

  obj <- structure(list(tree = tree, numSteps = numSteps), class = "cbsFs")

  ## Extract cp and cp.sign
  cp <- jumps(obj, unique=FALSE, sorted=FALSE)
  leaves <- .enumerate_splits(tree)
  hump.signs <- sign(as.numeric(sapply(leaves, function(x){
      data.tree::FindNode(tree, x)$cusum})))
  all.signs = lapply(hump.signs, function(hump.sign){
      cp.sign.add = c(1,-1)*hump.sign
  })
  cp.sign = unlist(all.signs)

  ## Handling when cp is NA, from jump.cbsFs
  cp.sign = cp.sign[which(!is.na(cp))]
  cp = cp[which(!is.na(cp))]

  ## Return with cp and cp.sign
  obj <- structure(list(y=y, tree = tree, numSteps = numSteps,
                        cp = cp, cp.sign=cp.sign), class = "cbsFs")

}

#' Get jumps from cbsFs objects
#'
#' Enumerates the jumps for circular binary segmentation. If \code{sorted}
#' is true, then only the unique changepoints are reported, and no
#' \code{NA} is recorded. If \code{sorted} is false, then each pair of
#' elements in the result denote the jump locations, and \code{NA} means
#' the left (first) element of the hump is the first index.
#'
#' @param obj cbs object
#' @param sorted boolean
#' @param unique prunes the changepoints to be unique. Not to be used for CBS
#' @param ... not used
#'
#' @return vector of all the jump locations
#' @export
jumps.cbsFs <- function(obj, sorted = T, unique=TRUE,cbs=FALSE, remove.na=FALSE,...){
  leaves <- .enumerate_splits(obj$tree)
  if(length(leaves) == 0) return(NA)

  res <- sapply(leaves, function(x){data.tree::FindNode(obj$tree, x)$breakpoint})
  res[1,] <- sapply(res[1,], function(x){ifelse(x>1, x-1, NA)})
  res <- as.numeric(res)
  if(remove.na){
      res <- res[!is.na(res)];
  }
  if(sorted) {
      sort(res)
  }
  if(unique){
      res = unique(res)
  }
  return(res)
}

#' Find breakpoints for circular binary segmentation
#'
#' Given a vector \code{y} (representating a segment of the data vector),
#' find the best place to split according to circular binary segmentation.
#' It returns a list with \code{breakpoint}, a vector of length 2 (for the
#' start and end of the hump), and \code{cusum} for the circular
#' cusum statistic.
#'
#' Here, the \code{breakpoint} is encoded so \code{breakpoint[1]-1}
#' is the last index of left shoulder (if larger than 0), and
#' \code{breakpoint[2]+1} is the first index of the right shoulder.
#'
#' @param y numeric vector
#'
#' @return a list
.find_breakpoint_cbs <- function(y){
  if(length(y) == 1) return(list(breakpoint = NA, cusum = 0))
  n <- length(y)

  vec <- cumsum(y)
  breakpoint <- .enumerate_breakpoints_cbs(n)
  cusum_vec <- apply(breakpoint, 1, .cusum_cbs, vec = vec)

  idx <- which.max(abs(cusum_vec))
  list(breakpoint = breakpoint[idx,], cusum = cusum_vec[idx])
}

#' Enumerate all the possible jump locations for circular binary segmentation
#'
#' Enumerates all the possible places the start and end of the hump
#' can occur for a segment of length \code{n}. This is essentially the same
#' as \code{t(cbind(combn(1:(n-1), 2), rbind(2:(n-1), rep(n, n-2)), rbind(1:n, 1:n)))}.
#' The first part, \code{combn(1:(n-1), 2)}, accounts for any humps with
#' a right shoulder and length larger than 1, possibly no left shoulder.
#' The second part, \code{rbind(2:n, rep(n, n-1)))} accounts for any humps with no
#' right shoulder and length larger than 1, but has a left shoulder.
#' The third part, \code{rbind(1:n, 1:n)}
#' accounts for humps of length 1.
#'
#' \code{start} shifts all the indices appropriately to the right.
#'
#' @param n integer
#' @param start integer
#'
#' @return a matrix
.enumerate_breakpoints_cbs <- function(n, start = 1){
  if(n == 1) return(matrix(c(1,1), nrow = 1))

  res <- cbind(rep(1:n, times = c((n-1), (n-1):1)),
        unlist(lapply(1:n, function(x){
          if(x == 1) x:(n-1) else x:n
        })))

  res + start - 1
}

#' Compute the circular cusum statistcs
#'
#' \code{x} are the breakpoints (vector of length 2), and \code{vec}
#' denotes the vector of cumulative sums. Output the circular cusum
#' statistic for this proposed split.
#'
#' \code{x} must be 2 integers where 1) the first element must be
#' larger than or equal to 1, 2) the second element must be smaller than
#' or equal to \code{length(vec)}, 3) the first element must be
#' smaller than or equal to the second element, and 4) the first element
#' cannot be 1 when the last element is \code{length(vec)}.
#'
#' @param x vector of 2 integers.
#' @param vec numeric vector
#'
#' @return a value
.cusum_cbs <- function(x, vec){
  stopifnot(length(x) == 2, x[1] >= 1, x[2] <= length(vec), x[1] <= x[2])
  stopifnot(all(x %% 1 == 0))
  stopifnot(!all(x[1]==1, x[2]==length(vec)) | length(vec) == 1)

  if(length(vec) == 1) return(0)

  n <- length(vec); m <- x[2]-x[1]+1
  sum1 <- vec[x[2]] - ifelse(x[1] > 1, vec[x[1]-1], 0)
  sum2 <- (vec[n] - vec[x[2]]) + ifelse(x[1] > 1, vec[x[1]-1], 0)

  const <- sqrt(1/(1/m + 1/(n-m)))
  const*(sum1/m - sum2/(n-m))
}

#' Create the children nodes for circular binary segmentation
#'
#' For a \code{node} (of class \code{Node}), where \code{node$breakpoint} is
#' a valid vector of breakpoints (of length 2), create up to three
#' children nodes according to the split. The first child, \code{left},
#' represents the left shoulder if it exists, \code{NA} otherwise. The
#' second child, \code{middle}, represents the hump. The third child,
#' \code{right}, represents the right shoulder if it exists, \code{NA} otherwise.
#' All three children are returned as a list.
#'
#' @param node \code{Node} object
#'
#' @return a list
.split_node_cbs <- function(node){
  if(any(is.na(node$breakpoint))) stop("node does not have a set breakpoint yet")
  stopifnot(length(node$breakpoint) == 2, node$breakpoint[1] >= node$start,
            node$breakpoint[2] <= node$end)

  if(node$breakpoint[1] > node$start){left <- .create_node(node$start, node$breakpoint[1] - 1)} else{left <- NA}
  middle <- .create_node(node$breakpoint[1], node$breakpoint[2])
  if(node$breakpoint[2] < node$end){right <- .create_node(node$breakpoint[2] + 1, node$end)} else{right <- NA}

  stopifnot(!any(is.na(left)) | !any(is.na(right)))

  list(left = left, middle = middle, right = right)
}

#' Placeholder function
#'
#' Included for technical coding reasons.
#'
#' @param x \code{Node} object
#'
#' @return FALSE
#' @export
is.na.Node <- function(x){FALSE}


jump_cusum.cbsFs <- function(obj, ...){
  jump_cusum(obj$tree, F)
}

