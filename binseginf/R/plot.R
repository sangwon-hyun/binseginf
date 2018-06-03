#' CpVector Class Constructor
#' 
#' Creates an instance of the CpVector, short for "Change Point Vector".
#' A CpVector objects has 3 elements, data (the realization from the
#' changepoint model), jump.height (the height of jumps in the true signal)
#' and jump.loc (the locations of the jumps in the true signal).
#' 
#' Here, CpVector is draws samples from a piecewise constant true signal from
#' points 1 to n. The jump.loc is a vector between 0 and 1 (inclusive and 
#' exclusive repsectively), and CpVector automatically scales jump.loc to
#' fit between 1 and n. Here, jump.loc must always be one less element than
#' jump.height.
#' 
#' Also, n must be large enough with respect to jump.loc that round(n*jump.loc)
#' must be distinct integers. Otherwise, an error will be thrown.
#' 
#' Here, func is a function that dictates the noise model. For example, the
#' default is rnorm.
#'
#' @param n the number of realizations
#' @param jump.height a numeric vector of the jump heights
#' @param jump.loc a numeric vector (between 0 (inclusive) and 1 (exclusive))
#' of the jump locations between 0 and 1.
#' @param func the function to dictate the noise model. The first parameter
#' of func must dictate how many samples are drawn.
#' @param ... additional parameters to pass into func.
#'
#' @return a CpVector instance
#' @export
CpVector <- function(n, jump.height = 0, jump.loc = NA, func = stats::rnorm, 
  ...){
  
  if(!is.numeric(n) | !is.numeric(jump.height))
    stop("n and jump.height must be numeric")
  if(!any(is.na(jump.loc))){
    if(!is.numeric(jump.loc)) stop("jump.loc must be numeric")
    if(length(jump.height) != length(jump.loc) + 1) 
      stop("jump.height must be one element more than jump.loc")
    if(min(jump.loc) < 0 | max(jump.loc) >= 1)
      stop("jump.loc must lie between 0 (inclusive) and 1 (exclusive)")
    if(!all(diff(jump.loc) > 0)) stop("jump.loc must be strictly increasing")
  } else {
    if(length(jump.height) != 1)
      stop("jump.height must be one element since jump.loc = NA")
    if(length(jump.loc) != 1)
      stop("jump.loc must be one element if it contains any NA")
  }
  
  if(n %% 1 != 0 | n < 0) stop("n must be a positive integer")
  
  jump.idx <- .computeJumpIdx(n, jump.loc)
  mean.vec <- .formMeanVec(n, jump.height, jump.idx)
  
  data <- mean.vec + func(n, ...)
  
  structure(list(data = data, jump.height = jump.height, jump.idx = jump.idx),
    class = "CpVector")
}

.computeJumpIdx <- function(n, jump.loc){
  if(any(is.na(jump.loc))) return(numeric(0))
  
  vec <- round(n*jump.loc)
  vec <- sapply(vec, function(x) {max(min(x,n),1)})
  
  if(length(vec) != length(unique(vec))) stop(paste("n is too small compared",
    "to jump.loc that changepoints are not unique"))
  
  vec
}

.formMeanVec <- function(n, jump.height, jump.idx){
  diff.vec <- diff(c(0,jump.idx,n))
  rep(jump.height, times = diff.vec)
}


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
