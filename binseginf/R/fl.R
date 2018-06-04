##' Wrapper for running FL.
##' @export
fl <- function(y, numSteps, sigma.add=NULL){
  if(numSteps >= length(y)) stop("numSteps must be strictly smaller than the length of y")
  if(numSteps <= 0) step("numSteps must be at least 1.")

  if(!is.null(sigma.add)){
      y.addnoise = rnorm(length(y), 0, sigma.add)
      y = y + y.addnoise
  }

  obj = genlassoinf::dualpathSvd2(y, maxsteps=numSteps,
                                  D=genlassoinf::makeDmat(length(y),ord=0))

  if(!is.null(sigma.add)){
      obj$noisy = TRUE
      obj$sigma.add = sigma.add
      obj$y.addnoise = y.addnoise
  }
  return(obj)
}
