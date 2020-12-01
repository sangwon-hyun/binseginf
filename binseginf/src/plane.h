#ifndef _PLANE_H
#define _PLANE_H

#include <Rcpp.h>
#include "radians.h"

class Circle{
public:
  Circle(Rcpp::NumericVector center_ = 0, double radius_ = 0);
  Rcpp::NumericVector center;
  double radius;
};

class Plane{
public:
  Plane(Rcpp::NumericVector a_ = 0, Rcpp::NumericVector b_ = 0);
  Rcpp::NumericVector a;
  Rcpp::NumericVector b;
  void c_normalize();
  void c_intersect_basis(const Rcpp::NumericVector &, const Rcpp::NumericVector &,
                         const Rcpp::NumericVector &);
  Rcpp::NumericVector c_point_on_plane();
  double c_distance_point_to_plane(const Rcpp::NumericVector &);
  Rcpp::NumericMatrix c_intersect_circle(const Circle &);
  void print();
};

double c_l2norm(const Rcpp::NumericVector & vec);
Rcpp::NumericVector c_quadratic(const double & a,
                                const double & b,
                                const double & c);
#endif
