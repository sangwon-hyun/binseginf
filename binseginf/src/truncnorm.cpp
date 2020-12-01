#include <Rcpp.h>
using namespace Rcpp;

// code modified from https://github.com/selective-inference/R-software/blob/master/forLater/maxZ/funs.constraints.R
// [[Rcpp::export]]
void c_gibbs_step(NumericVector & state, const NumericVector & direction,
                NumericVector & U, const NumericVector & alpha){
  int nstate = state.length();
  int nconstraint = U.length();

  /* Compute V=\eta^Ty */

  double value = sum(direction * state);

  /* Compute upper and lower bounds */

  double lower_bound = -1e12;
  double upper_bound = 1e12;
  double bound_val = 0;
  double tol = 1.e-7;

  for (int iconstraint = 0; iconstraint < nconstraint; iconstraint++) {

    bound_val = -U[iconstraint] / alpha[iconstraint] + value;

    if ((alpha[iconstraint] > tol) && (bound_val < upper_bound)) {
      upper_bound = bound_val;
    }
    else if ((alpha[iconstraint] < -tol) && (bound_val > lower_bound)) {
      lower_bound = bound_val;
    }
  }

  /* Ensure constraints are satisfied */

  if (lower_bound > value) {
    lower_bound = value - tol;
  }
  else if (upper_bound < value) {
    upper_bound = value + tol;
  }

  /* Now, take a step */

  double tnorm; /* the 1D gaussian variable */
  double cdfU, cdfL, unif; /* temp variables */



  if (upper_bound < -10) {

    /* use Exp approximation */
    /* the approximation is that */
    /* Z | lower_bound < Z < upper_bound */
    /* is fabs(upper_bound) * (upper_bound - Z) = E approx Exp(1) */
    /* so Z = upper_bound - E / fabs(upper_bound) */
    /* and the truncation of the exponential is */
    /* E < fabs(upper_bound - lower_bound) * fabs(upper_bound) = D */

    /* this has distribution function (1 - exp(-x)) / (1 - exp(-D)) */
    /* so to draw from this distribution */
    /* we set E = - log(1 - U * (1 - exp(-D))) where U is Unif(0,1) */
    /* and Z (= tnorm below) is as stated */

    unif = R::runif(0, 1) * (1 - exp(-fabs((lower_bound - upper_bound) * upper_bound)));
    tnorm = (upper_bound + log(1 - unif) / fabs(upper_bound));
  }
  else if (lower_bound > 10) {

    /* here Z = lower_bound + E / fabs(lower_bound) (though lower_bound is positive) */
    /* and D = fabs((upper_bound - lower_bound) * lower_bound) */

    unif = R::runif(0., 1.) * (1 - exp(-fabs((upper_bound - lower_bound) * lower_bound)));
    tnorm = (lower_bound - log(1 - unif) / lower_bound);
  }
  else if (lower_bound < 0) {
    cdfL = R::pnorm(lower_bound, 0., 1., 1, 0);
    cdfU = R::pnorm(upper_bound, 0., 1., 1, 0);
    unif = R::runif(0., 1.) * (cdfU - cdfL) + cdfL;
    if (unif < 0.5) {
      tnorm = R::qnorm(unif, 0., 1., 1, 0);
    }
    else {
      tnorm = -R::qnorm(1-unif, 0., 1., 1, 0);
    }
  }
  else {
    cdfL = R::pnorm(-lower_bound, 0., 1., 1, 0);
    cdfU = R::pnorm(-upper_bound, 0., 1., 1, 0);
    unif = R::runif(0., 1.) * (cdfL - cdfU) + cdfU;
    if (unif < 0.5) {
      tnorm = -R::qnorm(unif, 0., 1., 1, 0);
    }
    else {
      tnorm = R::qnorm(1-unif, 0., 1., 1, 0);
    }
  }

  /* Now update the state and U */

  double delta = tnorm - value;

  for (int istate = 0; istate < nstate; istate++) {
    state[istate] += delta * direction[istate];
  }
  for (int iconstraint = 0; iconstraint < nconstraint; iconstraint++) {
    U[iconstraint] += delta * alpha[iconstraint] ;
  }
}

// [[Rcpp::export]]
NumericMatrix c_sample_truncnorm_white(NumericVector & state, NumericVector & U,
                            const NumericMatrix & directions, const NumericMatrix & alphas,
                            const int & burnin, const int & ndraw){
  int which_direction;
  int ndirection = directions.ncol();
  int nstate = state.length();

  NumericVector direction;
  NumericVector alpha;
  NumericMatrix output(nstate, ndraw);

  for (int iter_count = 0; iter_count < burnin + ndraw; iter_count++) {

    which_direction = (int) floor(R::runif(0, 1) * ndirection);
    direction = directions(_ , which_direction);
    alpha = alphas(_, which_direction);

    /* take a step, which implicitly updates `state` and `U` */

    c_gibbs_step(state, direction, U, alpha);

    if (iter_count >= burnin) {
      output(_, iter_count-burnin) = state;
    }
  }

  return output;
}
