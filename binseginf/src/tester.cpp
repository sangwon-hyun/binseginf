#include "interval.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix c_intersect_circle_tester(Rcpp::NumericVector a,
                                              Rcpp::NumericVector b,
                                              Rcpp::NumericVector center,
                                              double radius){
  Plane plane(a, b);
  Circle circle(center, radius);
  Rcpp::NumericMatrix mat = plane.c_intersect_circle(circle);
  return mat;
}

// [[Rcpp::export]]
Rcpp::List c_intersect_basis_tester(Rcpp::NumericVector a,
                                              Rcpp::NumericVector b,
                                              Rcpp::NumericVector y,
                                              Rcpp::NumericVector v,
                                              Rcpp::NumericVector w){
  Plane plane(a, b);
  plane.c_intersect_basis(y, v, w);

  return Rcpp::List::create(Rcpp::Named("a") = plane.a,
                            Rcpp::Named("b") = plane.b);
}

// [[Rcpp::export]]
double c_euclidean_to_radian_tester(Rcpp::NumericVector center,
                                              double radius,
                                              Rcpp::NumericVector point){
  Circle circle(center, radius);
  double val = c_euclidean_to_radian(circle, point);
  return val;
}

// c_euclidean_to_radian_tester(c(1,1), sqrt(2), c(1, 1+sqrt(2)))
