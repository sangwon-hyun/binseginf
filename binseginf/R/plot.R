#' Plot function for CpVector
#'
#' @param x the CpVector instance
#' @param col.dots the color of the data points (x$data)
#' @param col.line the color of the true signal (x$jump.height and x$jump.idx)
#' @param ylim the graphical parameter ylim
#' @param ... additional graphical parameters such as pch, par, etc.
#'
#' @return void
#' @export
#'
#' @examples
#' set.seed(10)
#' res <- CpVector(n = 100, jump.height = 1:4, jump.loc = c(.25,.5,.75))
#' defaultPlotDefaults()
#' plot(res, ylab = "Value", xlab = "n")
plot.CpVector <- function(x, col.dots = defaultColors(-1, 0.75), 
  col.line = defaultColors(1), ylim = NA, ...){
  
  if(is.na(ylim)) ylim <- c(min(x$data, x$jump.height), 
    max(x$data, x$jump.height))
  
  graphics::plot(x$data, col = col.dots, ylim = ylim, ...)
  
  lis <- .splitChangepoints(length(x$data), x$jump.height, x$jump.idx)
  
  for(i in 1:length(lis)){
    graphics::lines(lis[[i]]$x, lis[[i]]$y, col = col.line, ...)
  }
  
  invisible()
}

.splitChangepoints <- function(n, jump.height, jump.idx){
  jump.idx <- c(0, jump.idx, n)
  
  lis <- vector("list", length(jump.idx) - 2)
  for(i in 2:(length(jump.idx))){
    x <- (jump.idx[i-1]+1):min(jump.idx[i]+1,n)
    y <- rep(jump.height[i-1], length(x))
    lis[[i-1]] <- list(x = x, y = y)
  }
  
  lis
}

