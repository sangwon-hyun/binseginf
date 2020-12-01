#ifndef _RADIANS_H
#define _RADIANS_H

#include <Rcpp.h>

Rcpp::NumericVector c_unique_sort_native(const Rcpp::NumericVector & x);
Rcpp::NumericVector c_construct_midpoints(const Rcpp::NumericVector & x);
Rcpp::IntegerVector c_which_native(const Rcpp::LogicalVector & x);
Rcpp::IntegerMatrix c_consecutive_true(const Rcpp::LogicalVector & vec);
Rcpp::NumericVector c_unlist_native(const Rcpp::List & list);
Rcpp::LogicalVector c_theta_in_matrix(const double & x, const Rcpp::NumericMatrix & mat);
Rcpp::LogicalVector c_theta_in_all_matrix(const double & x, const Rcpp::List & list);
Rcpp::NumericMatrix c_intersect_intervals(const Rcpp::List & list);

#endif
