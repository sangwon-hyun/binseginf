#include "interval.h"

double c_pi() { return std::atan(1)*4; }

// from https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T> int c_sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

// from https://stackoverflow.com/questions/21609934/ordering-permutation-in-rcpp-i-e-baseorder
Rcpp::IntegerVector c_order(Rcpp::NumericVector x) {
  Rcpp::NumericVector sorted = clone(x).sort();
  return match(sorted, x);
}

double c_euclidean_to_radian(const Circle & circle,
                             const Rcpp::NumericVector & point){
  Rcpp::NumericVector tmp = Rcpp::no_init(1);
  tmp = point[0];
  double point1 = Rcpp::as<double>(tmp);
  tmp = point[1];
  double point2 = Rcpp::as<double>(tmp);
  tmp = circle.center[0];
  double circle1 = Rcpp::as<double>(tmp);
  tmp = circle.center[1];
  double circle2 = Rcpp::as<double>(tmp);

  if(pow(point1 - circle1, 2) + pow(point2 - circle2, 2) -
     pow(circle.radius, 2) > 1e-3) {
    throw std::runtime_error("point not on circle");
  }
  double val;

  if(point2 != 0){
    val = atan(point1/point2);
  } else {
    val = atan(-circle2/circle1);
  }

  return val;
}

// [[Rcpp::export]]
double c_initial_theta(const Rcpp::NumericVector & y,
                       const Rcpp::NumericVector & v,
                       const Rcpp::NumericVector & w){
  double top = 0;
  double bottom = 0;
  int len = y.length();

  for(int i = 0; i < len; i++){
    top += y[i] * w[i];
    bottom += y[i] * v[i];
  }

  double val = atan(-top/bottom);
  return val;
}

// [[Rcpp::export]]
Rcpp::NumericVector c_basic_interval(const Rcpp::NumericVector & endpoints,
                                     const double & theta){
  if(endpoints.length() != 2){
    throw std::runtime_error("endpoint is not a vector of 2");
  }

  Rcpp::NumericVector endpoints_copy = Rcpp::no_init(2);
  std::copy(endpoints.begin(), endpoints.end(), endpoints_copy.begin());

  Rcpp::NumericVector tmp1 = Rcpp::no_init(1);
  Rcpp::NumericVector tmp2 = Rcpp::no_init(1);
  tmp1 = endpoints[0]; tmp2 = endpoints[1];
  if(fabs(Rcpp::as<double>(tmp1)) > c_pi()/2 |
     fabs(Rcpp::as<double>(tmp2)) > c_pi()/2){
    throw std::runtime_error("endpoint out of range");
  }

  endpoints_copy.sort(false);

  tmp1 = endpoints_copy[0];
  tmp2 = endpoints_copy[1];

  if(Rcpp::as<double>(tmp1) > theta | theta > Rcpp::as<double>(tmp2)){
    endpoints_copy[0] = Rcpp::as<double>(tmp2);
    endpoints_copy[1] = Rcpp::as<double>(tmp1) + c_pi();
  }
  return endpoints_copy;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix c_partition_interval(const Rcpp::NumericVector & interval){
  Rcpp::NumericVector tmp = Rcpp::no_init(1);
  tmp = interval[0];
  double a = Rcpp::as<double>(tmp);
  tmp = interval[1];
  double b = Rcpp::as<double>(tmp);

  if(a > b){
    throw std::runtime_error("interval not sorted");
  }

  Rcpp::NumericVector vec = Rcpp::NumericVector::create(-2*c_pi(), -3*c_pi()/2, -c_pi(), -c_pi()/2, 0,
                                                        c_pi()/2, c_pi(), 3*c_pi()/2, 2*c_pi());

  Rcpp::IntegerVector idx1 = c_which_native(vec >= a);
  Rcpp::IntegerVector idx2 = c_which_native(vec <= b);
  Rcpp::IntegerVector idx = Rcpp::intersect(idx1, idx2);

  Rcpp::NumericVector vec2 = Rcpp::no_init(idx.length() + 2);
  vec2[0] = a; vec2[vec2.length()-1] = b;
  if(idx.length() > 0){
    for(int i = 1; i < vec2.length()-1; i++){
      vec2[i] = vec[idx[i - 1] - 1];
    }
  }

  Rcpp::NumericMatrix mat = Rcpp::no_init(vec2.length()-1, 2);
  for(int i = 0; i < mat.nrow(); i++){
    mat(i,0) = vec2[i];
    mat(i,1) = vec2[i+1];
  }

  Rcpp::NumericVector idx_nv = Rcpp::no_init(mat.nrow());
  Rcpp::NumericVector tmp_nv = Rcpp::no_init(1);
  double tmp_d;
  for(int i = 0; i < idx_nv.length(); i++){
    tmp_nv = (mat(i,0) + mat(i,1))/2;
    tmp_d = Rcpp::as<double>(tmp_nv);
    idx_nv[i] = c_sgn(tmp_d) * std::ceil(fabs(tmp_d)/(c_pi()/2));
  }

  for(int i = 0; i < mat.nrow(); i++){
    tmp_nv = idx_nv[i];
    tmp_d = Rcpp::as<double>(tmp_nv);
    if(fabs(tmp_d) > 1){
      mat(i,0) = mat(i,0) - c_sgn(tmp_d)*c_pi();
      mat(i,1) = mat(i,1) - c_sgn(tmp_d)*c_pi();
    }
  }

  Rcpp::NumericMatrix mat2;
  if(mat.nrow() == 1){
    return mat;
  } else {
    idx = c_order(mat(Rcpp::_, 0));
    mat2 = Rcpp::no_init(mat.nrow(), mat.ncol());
    for(int i = 0; i < mat.nrow(); i++){
      mat2(i,0) = mat(idx[i]-1,0);
      mat2(i,1) = mat(idx[i]-1,1);
    }
    mat = mat2;
  }

  tmp_nv = mat(Rcpp::_,1) - mat(Rcpp::_,0);
  idx = c_which_native(tmp_nv > 1e-6);
  mat2 = Rcpp::no_init(idx.length(), mat.ncol());
  for(int i = 0; i < idx.length(); i++){
    mat2(i,0) = mat(idx[i]-1,0);
    mat2(i,1) = mat(idx[i]-1,1);
  }

  return mat2;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix c_interval(const Rcpp::NumericVector & endpoints,
                               const double & initial_theta){
  Rcpp::NumericMatrix mat;

  if(is_true(any(is_na(endpoints)))){
    mat = Rcpp::no_init(1,2);
    mat(0,0) = -c_pi()/2;
    mat(0,1) = c_pi()/2;
  } else {
    Rcpp::NumericVector interval = c_basic_interval(endpoints, initial_theta);
    mat = c_partition_interval(interval);
  }

  return mat;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix c_form_interval(const Rcpp::NumericVector & a,
                                    const Rcpp::NumericVector & b,
                                    const Rcpp::NumericVector & y,
                                    const Rcpp::NumericVector & v,
                                    const Rcpp::NumericVector & w){

  Plane plane = Plane(a, b);
  plane.c_intersect_basis(y, v, w);

  Rcpp::NumericMatrix mat;
  if(is_true(any(is_na(plane.a)))){
    mat = Rcpp::no_init(1,2);
    mat(0,0) = -c_pi()/2;
    mat(0,1) = c_pi()/2;
    return mat;
  }

  Rcpp::NumericVector center = Rcpp::no_init(2);
  std::fill(center.begin(), center.end(), 0);
  int n = y.length();
  for(int i = 0; i < n; i++){
    center[0] -= y[i]*v[i];
    center[1] -= y[i]*w[i];
  }
  double radius = c_l2norm(center);
  Circle circle = Circle(center, radius);
  double dis = plane.c_distance_point_to_plane(center);

  if(dis >= radius){
    mat = Rcpp::no_init(1,2);
    mat(0,0) = -c_pi()/2;
    mat(0,1) = c_pi()/2;
  } else {
    Rcpp::NumericMatrix mat2 = plane.c_intersect_circle(circle);
    Rcpp::NumericVector endpoints = Rcpp::no_init(2);
    endpoints[0] = c_euclidean_to_radian(circle, mat2(Rcpp::_,0));
    endpoints[1] = c_euclidean_to_radian(circle, mat2(Rcpp::_,1));
    double init_theta = c_initial_theta(y, v, w);
    mat = c_interval(endpoints, init_theta);
  }

  return mat;
}


// c_partition_interval(c(pi/3, 3*pi/4))



