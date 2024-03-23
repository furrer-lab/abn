#include <Rcpp.h>
#include <math.h>

//' @title logit functions
//'
//' @description transform \code{x} either via the logit, or expit.
//'
//'
//' @param x a numeric vector
//' @returns a numeric vector
//' @export
//' @rdname logit_cpp
// [[Rcpp::export]]
Rcpp::NumericVector logit_cpp(Rcpp::NumericVector x) {
  int n = x.size();
  Rcpp::NumericVector result(n);

  for(int i = 0; i < n; ++i) {
    result[i] = std::log(static_cast<double>( x[i] / (1.0 - x[i]) ));
  }

  return result;
}

//' @title expit function
//'
//' @description transform \code{x} either via the logit, or expit.
//'
//'
//' @param x a numeric vector
//' @returns a numeric vector
//' @export
//' @rdname expit_cpp
// [[Rcpp::export]]
Rcpp::NumericVector expit_cpp(Rcpp::NumericVector x) {
  int n = x.size();
  Rcpp::NumericVector result(n);

  for (int i=0; i < n; ++i) {
    result[i] = 1.0 / (1.0 + std::exp(static_cast<double>( (-1.0 * x[i]))));
  }
  return result;
}
