#ifndef _INTERVAL_H
#define _INTERVAL_H

#include <Rcpp.h>
#include <math.h>
#include "plane.h"

double c_euclidean_to_radian(const Circle & circle,
                             const Rcpp::NumericVector & point);
double c_initial_theta(const Rcpp::NumericVector & y,
                       const Rcpp::NumericVector & v,
                       const Rcpp::NumericVector & w);
Rcpp::NumericVector c_basic_interval(const Rcpp::NumericVector & endpoints,
                                     const double & theta);
Rcpp::NumericMatrix c_partition_interval(const Rcpp::NumericVector & interval);
Rcpp::NumericMatrix c_interval(const Rcpp::NumericVector & endpoints,
                               const double & initial_theta);
Rcpp::NumericMatrix c_form_interval(const Rcpp::NumericVector & a,
                                    const Rcpp::NumericVector & b,
                                    const Rcpp::NumericVector & y,
                                    const Rcpp::NumericVector & v,
                                    const Rcpp::NumericVector & w);

#endif
