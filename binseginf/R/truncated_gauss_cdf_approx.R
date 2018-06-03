# The below functions are approximations for the truncated Gaussian cdf

# Returns Prob(Z>z | Z in [a,b]), where mean can be a vector, based on
# A UNIFORM APPROXIMATION TO THE RIGHT NORMAL TAIL INTEGRAL, W Bryc
# Applied Mathematics and Computation
# Volume 127, Issues 23, 15 April 2002, Pages 365--374
# https://math.uc.edu/~brycw/preprint/z-tail/z-tail.pdf

.truncated_gauss_cdf_approx <- function(z, mean=0, sd=1, a, b) {
  z <- (z-mean)/sd
  a <- (a-mean)/sd
  b <- (b-mean)/sd
  n <- length(z)
  
  term1 <- exp(z*z)
  o <- a > -Inf
  term1[o] <- ff(a[o])*exp(-(a[o]^2 - z[o]^2)/2)
  term2 <- rep(0,n)
  oo <- b < Inf
  term2[oo] <- ff(b[oo])*exp(-(b[oo]^2 - z[oo]^2)/2)
  p <- (ff(z)-term2)/(term1-term2)
  
  pmin(1, pmax(0,p))
}

ff <- function(z) {
  return((z^2+5.575192695*z+12.7743632)/
           (z^3*sqrt(2*pi)+14.38718147*z^2+31.53531977*z+2*12.77436324))
}