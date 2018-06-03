#' List comparisons for circular binary segmentation fixed threshold object
#' 
#' Returns a list with number of elements equal to the number of nodes
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
#' If an element in the returned list has \code{winning} set to 4 \code{NA}'s,
#' it means for that particular split, all possible future splits resulted in
#' a cusum statistic lower than the desired threshold.
#'
#' @param obj cbsFt object
#'
#' @return a list of lists of matrices
.list_comparison.cbsFt <- function(obj){
  nodes <- obj$tree$Get("active")
  active_vec <- names(sort(nodes))
  inactive_vec <- names(nodes)[!names(nodes) %in% active_vec]
  
  active_vec <- .remove_singleton_nodes(active_vec)
  inactive_vec <- .remove_singleton_nodes(inactive_vec)
  num_nodes <- length(active_vec) + length(inactive_vec)
  
  comp_lis <- vector("list", num_nodes)
  u_lis <- vector("list", num_nodes)
  tree2 <- .create_node(1, obj$tree$end)
  
  if(length(active_vec) > 0){
    for(i in 1:length(active_vec)){
      breakpoint <- data.tree::FindNode(obj$tree, active_vec[i])$breakpoint
      comp_lis[[i]] <- .form_comparison_cbs(tree2, active_vec[i], breakpoint, excluding = F)
      
      node_selected <- data.tree::FindNode(tree2, active_vec[i])
      node_selected$breakpoint <- breakpoint
      
      node_pairs <- .split_node_cbs(node_selected)
      if(!any(is.na(node_pairs$left))) node_selected$AddChildNode(node_pairs$left)
      node_selected$AddChildNode(node_pairs$middle)
      if(!any(is.na(node_pairs$right))) node_selected$AddChildNode(node_pairs$right)
    }
  }
 
  if(length(inactive_vec) > 0){
    for(i in 1:length(inactive_vec)){
      comp_lis[[i+length(active_vec)]] <- .form_comparison_cbs(tree2, inactive_vec[i], 
                                                               NA, excluding = F)
      comp_lis[[i+length(active_vec)]]$winning[1:4] <- NA
    }
  }

  comp_lis
}

#' Remove singleton nodes
#' 
#' This a specialized function to remove nodes (based on their names)
#' that start and end with the same index. This is to prevent degenercy 
#' in the methods.
#'
#' @param node_vec Vector of node names
#'
#' @return vector of node names
.remove_singleton_nodes <- function(node_vec){
  if(length(node_vec) == 0) return(node_vec)
  vec <- sapply(node_vec, .get_startEnd)
  idx <- which(apply(vec, 2, function(x){ifelse(x[1]==x[2], FALSE, TRUE)}))
  node_vec[idx]
}