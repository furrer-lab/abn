#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//' @title Rank of a matrix
//' @description similar to \code{base::rank}
//' @keywords internal
//' @returns an integer
//' @export
// [[Rcpp::export]]

int rank_cpp(arma::mat A)
{
  int out;
  out = rank(A);
  //return
  return out;

}
