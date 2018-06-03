.list_comparison.bsFs <- function(obj){
  is_valid(obj)

  numSteps <- obj$numSteps
  comp_lis <- vector("list", numSteps)
  active_vec <- .enumerate_splits(obj$tree)
  tree2 <- .create_node(1, obj$tree$end)

  for(i in 1:numSteps){
    breakpoint <- data.tree::FindNode(obj$tree, active_vec[i])$breakpoint
    comp_lis[[i]] <- .form_comparison(tree2, active_vec[i], breakpoint)

    node_selected <- data.tree::FindNode(tree2, active_vec[i])
    node_selected$breakpoint <- breakpoint
    
    node_pairs <- .split_node(node_selected)
    node_selected$AddChildNode(node_pairs$left)
    node_selected$AddChildNode(node_pairs$right)
  }

  comp_lis
}

#' List comparisons for circular binary segmentation fixed steps object
#' 
#' Returns a list with number of elements equal to the number of splits
#' in \code{obj}. Each element of the list is a list with elements
#' \code{winning} and \code{losing}, each a matrix with 4 columns.
#' \code{winning} only has 1 row, denoting the comparison that was selected
#' at a particular step of ciruclar binary segmentation.
#' \code{losing} contains all the other comparisons not chosen.
#' 
#' The 4 columns represent (from left to right) the first (left) index of the left shoulder,
#' the first (left) index of the hump, the last (right) index of the hump,
#' and the last (right) index of the right shoulder. 
#'
#' @param obj 
#'
#' @return a list of lists of matrices
.list_comparison.cbsFs <- function(obj){
  numSteps <- obj$numSteps
  comp_lis <- vector("list", numSteps)
  active_vec <- .enumerate_splits(obj$tree)
  tree2 <- .create_node(1, obj$tree$end)
  
  for(i in 1:numSteps){
    breakpoint <- data.tree::FindNode(obj$tree, active_vec[i])$breakpoint
    comp_lis[[i]] <- .form_comparison_cbs(tree2, active_vec[i], breakpoint)
    
    node_selected <- data.tree::FindNode(tree2, active_vec[i])
    node_selected$breakpoint <- breakpoint
    
    node_pairs <- .split_node_cbs(node_selected)
    if(!any(is.na(node_pairs$left))) node_selected$AddChildNode(node_pairs$left)
    node_selected$AddChildNode(node_pairs$middle)
    if(!any(is.na(node_pairs$right))) node_selected$AddChildNode(node_pairs$right)
  }
  
  comp_lis
}




