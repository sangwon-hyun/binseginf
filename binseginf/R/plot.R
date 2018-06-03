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

#' Default Colors
#'
#' @param vec a vector of integers ranging between 
#' .defaultColorRange() (inclusive)
#' @param transparency a numeric between 0 and 1 (inclusive) to set the
#' transparency of the colors
#'
#' @return void
#' @export
defaultColors <- function(vec = NA, transparency = 1){
  if(!is.na(vec)) stopifnot(is.numeric(vec), all(vec %% 1 == 0))
  stopifnot(is.numeric(transparency), length(transparency) == 1,
    transparency >= 0, transparency <= 1)
  
  if(all(is.na(vec))) vec <- .defaultColorRange()
  trans <- round(transparency*255)
  
  color.numeric <- .defaultColorRange()
  color.rgb <- c(
    grDevices::rgb(150, 150, 150, trans, max = 255), #gray
    grDevices::rgb(0, 0, 0, trans, max = 255), #black
    grDevices::rgb(255, 0, 0, trans, max = 255), #red
    grDevices::rgb(0, 205, 0, trans, max = 255), #green3
    grDevices::rgb(0, 0, 255, trans, max = 255)) #blue
  
  plyr::mapvalues(vec, from = color.numeric, to = color.rgb, warn_missing = F)
}

.defaultColorRange <- function(){-1:3}

.demoColors <- function(){
  col.vec <- defaultColors()
  defaultPlotDefaults()
  graphics::plot(x = 1:length(col.vec), y = rep(1, length(col.vec)), 
    col = col.vec, cex = 4, xlab = "", ylab = "", yaxt = "n")
}

.defaultFont <- function(){"sans"}

.defaultCex <- function(){1.5}

.defaultPch <- function(){16}

.defaultLwd <- function(){3}


#' Default Plot Settings
#' 
#' Removes all the plots and currently rendering pdf figures and sets
#' all the default graphical parameters. Useful for standardizing the figures.
#' 
#' This sets the graphical parameters, "pch", "cex", "family" and "mar".
#'
#' @return void
#' @export
defaultPlotDefaults <- function(){
  res <- names(grDevices::dev.cur())
  if(res != "null device" & res != "pdf") grDevices::graphics.off()
  graphics::par(pch = .defaultPch(), cex = .defaultCex(), 
    lwd = .defaultLwd(), family = .defaultFont(), mar = c(4,4,1,1))
  
  invisible()
}
