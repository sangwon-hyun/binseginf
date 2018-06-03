circularBinSeg_fixedThres<- function(y, thres){
  stopifnot(thres > 0)
  
  #initialization
  n <- length(y); tree <- .create_node(1, n); vec <- cumsum(y)
  
  bool <- TRUE; round <- 1
  while(bool){
    bool <- FALSE
    
    leaves.names <- .get_leaves_names(tree)
    for(i in 1:length(leaves.names)){
      leaf <- data.tree::FindNode(tree, leaves.names[i])
      
      if(is.na(leaf$cusum) & leaf$start != leaf$end){
        res <- .find_breakpoint_cbs(y[leaf$start:leaf$end])
        
        leaf$breakpoint <- res$breakpoint+leaf$start-1; leaf$cusum <- res$cusum
        
        if(abs(leaf$cusum) > thres){
          bool <- TRUE
          leaf$active <- round
          node.pairs <- .split_node_cbs(leaf)
          if(!any(is.na(node.pairs$left))) leaf$AddChildNode(node.pairs$left)
          leaf$AddChildNode(node.pairs$middle)
          if(!any(is.na(node.pairs$right))) leaf$AddChildNode(node.pairs$right)
        }
      }
    }
    
    round <- round+1
  }
  
  obj <- structure(list(tree = tree, thres = thres), class = "cbsFt")
}

#' Get jumps from cbsFt objects
#' 
#' Enumerates the jumps for circular binary segmentation. 
#'
#' @param obj cbs object
#' @param ... not used
#'
#' @return vector of all the jump locations
#' @export
jumps.cbsFt <- function(obj, ...){
  jumps.cbsFs(obj, ...)
}