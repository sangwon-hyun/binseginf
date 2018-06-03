#' Polyhedra for circular binary segmentation fixed threshold
#' 
#' The polyhedra has elements \code{gamma} (a matrix) and \code{u}
#' (a vector).
#'
#' @param obj cbsFt object
#' @param ... void, not used
#'
#' @return a polyhedra 
#' @export
polyhedra.cbsFt <- function(obj, ...){
  n <- .get_startEnd(obj$tree$name)[2] 
  comp_lis <- .list_comparison(obj)
  numNodes <- length(comp_lis)
  gamma_row_lis <- vector("list", numNodes)
  u_lis <- vector("list", numNodes)
  
  sign_vec <- .get_signs_cbsFt(obj)
  
  for(i in 1:numNodes){
    add <- !all(is.na(comp_lis[[i]]$winning))
    gamma_row_lis[[i]] <- .gammaRows_from_comparisons_cbsfs(comp_lis[[i]]$winning,
                                                            comp_lis[[i]]$losing, sign_vec[i], n,
                                                            add)
    if(all(is.na(comp_lis[[i]]$winning))){
      stopifnot(nrow(gamma_row_lis[[i]]) %% 2 == 0)
      u_lis[[i]] <- c(rep(-obj$thres, nrow(gamma_row_lis[[i]])/2),
                      rep(-obj$thres, nrow(gamma_row_lis[[i]])/2))
    } else {
      row_minus_one <- ifelse(length(nrow(gamma_row_lis[[i]])) > 0, nrow(gamma_row_lis[[i]]) - 1, 0)
      u_lis[[i]] <- c(rep(0, row_minus_one), obj$thres)
    }
  }
  
  mat <- do.call(rbind, gamma_row_lis)
  u <- do.call(c, u_lis)
  polyhedra(obj = mat, u = u)
}

#' Get signs of a circular binary segmentation fixed threshold
#' 
#' The reason this is a private function is because it's somewhat non-standard.
#' It's unclear what "order" the jumps of a fixed threshold algorithm should
#' be. In this procedure, we order the jumps by their "rounds" (in 
#' \code{obj$tree$Get(\"active\")}), but it's unclear how the jumps are ordered
#' within the same round. 
#' 
#' The returned vector has length equal to the number of nodes in the 
#' circular binary segmentation tree, and all the active nodes are in
#' the front of the vector. All the inactive nodes are towards the back of
#' the returned vector with element 0.
#'
#' @param obj cbsFt object
#'
#' @return a vector containing (1,-1,0)
.get_signs_cbsFt <- function(obj){
  nodes <- obj$tree$Get("active")
  active_vec <- names(sort(nodes))
  
  sign_vec <- rep(0, length(nodes))
  if(length(active_vec) == 0) return(sign_vec)
  
  for(i in 1:length(active_vec)){
    node_selected <- data.tree::FindNode(obj$tree, active_vec[i])
    sign_vec[i] <- sign(node_selected$cusum)
  }
  
  sign_vec
}